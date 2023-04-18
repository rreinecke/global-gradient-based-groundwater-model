#ifndef TESTING_SIMPLEDATAREADER_HPP
#define TESTING_SIMPLEDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"

namespace GlobalFlow {
namespace DataProcessing {


class SimpleDataReader : public DataReader {
    public:
        SimpleDataReader(int step) { stepMod = step; }

        virtual void readData(Simulation::Options op) {
            LOG(userinfo) << "Building the initial model layer";
            std::vector<std::vector<int>> grid;
            grid = readGrid(nodes,
                            buildDir(op.getNodesDir()),
                            op.getNumberOfNodes(),
                            op.getNumberOfRows(),
                            op.getNumberOfCols(),
                            op.getInitialK(),
                            op.getInitialHead(),
                            op.getAquiferDepth()[0],
                            op.getAnisotropy(),
                            op.getSpecificYield(),
                            op.getSpecificStorage(),
                            op.getEdgeLengthLeftRight(),
                            op.getEdgeLengthFrontBack(),
                            op.isConfined(0),
                            op.isDensityVariable(),
                            op.getInitialZetas(),
                            op.getMaxTipSlope(),
                            op.getMaxToeSlope(),
                            op.getMinDepthFactor(),
                            op.getSlopeAdjFactor(),
                            op.getVDFLock(),
                            op.getDensityZones());

            LOG(userinfo) << "Building the bottom layers";
            DataProcessing::buildBottomLayers(nodes,
                                              op.getNumberOfLayers(),
                                              op.getConfinements(),
                                              op.getAquiferDepth());

            LOG(userinfo) << "Reading hydraulic parameters";
            readConduct(buildDir(op.getLithology()));
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Initializing head";
            readInitialHeads((buildDir(op.getInitialHeadsDir())));

            LOG(userinfo) << "Defining rivers";
            readRiver(buildDir(op.getKRiver()));

            LOG(userinfo) << "Connecting the layers";
            DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodes(), op.getNumberOfLayers(),
                                        op.getGHBConduct(),
                                        op.getBoundaryCondition());
        }

    private:
        template<class T>
        using Matrix = std::vector<std::vector<T>>;

        Matrix<int>
        readGrid(NodeVector nodes,
                 std::string path,
                 int numberOfNodes,
                 int numberOfRows,
                 int numberOfCols,
                 double defaultK,
                 double initialHead,
                 double aquiferDepth,
                 double anisotropy,
                 double specificYield,
                 double specificStorage,
                 double edgeLengthLeftRight,
                 double edgeLengthFrontBack,
                 bool confined,
                 bool isDensityVariable,
                 vector<double> initialZetas,
                 double maxTipSlope,
                 double maxToeSlope,
                 double minDepthFactor,
                 double slopeAdjFactor,
                 double vdfLock,
                 vector<double> densityZones) {
            Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

            io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "X", "Y", "cell_area", "col", "row");

            double x{0};
            double y{0};
            double area{0};
            int globid{0};
            int i{0};
            int row{0};
            int col{0};
            lookupSpatIDtoID.reserve(numberOfNodes);
            vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
            vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);
            vector<Model::quantity<Model::Meter>> initialZetasDim;
            for (int i = 0; i < initialZetas.size(); i++) {
                initialZetasDim.push_back(initialZetas[i] * Model::si::meter);
            }
            while (in.read_row(globid, x, y, area, col, row)) {
                out[col][row] = i;
                //area is in km needs to be in m
                nodes->emplace_back(new Model::StandardNode(nodes,
                                                            x,
                                                            y,
                                                            area * Model::si::square_meter,
                                                            edgeLengthLeftRight * Model::si::meter,
                                                            edgeLengthFrontBack * Model::si::meter,
                                                            (unsigned long) globid,
                                                            i,
                                                            defaultK * (Model::si::meter / Model::day),
                                                            initialHead * Model::si::meter,
                                                            stepMod,
                                                            aquiferDepth,
                                                            anisotropy,
                                                            specificYield,
                                                            specificStorage,
                                                            confined,
                                                            isDensityVariable,
                                                            initialZetasDim,
                                                            delnus,
                                                            nusInZones,
                                                            maxTipSlope,
                                                            maxToeSlope,
                                                            minDepthFactor,
                                                            slopeAdjFactor,
                                                            vdfLock * Model::si::meter));
                lookupSpatIDtoID[globid] = i;
                i++;
            }

            return out;
        }

        void readConduct(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                if (data > 10) {
                    LOG(debug) << "Very high conductance value at spatID "<<pos<<". Possible Data Error";
                }
                nodes->at(pos)->setK(data * (Model::si::meter / Model::day));
            });
        };

        void
        readElevation(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->setElevation(data * Model::si::meter);
            });
        };

        void readRiver(std::string path) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "Head", "Bottom", "Conduct");
            int spatID{0};
            double head{0};
            double conduct{0};
            double bottom{0};

            while (in.read_row(spatID, head, bottom, conduct)) {
                int i = 0;
                try {
                    i = lookupSpatIDtoID.at(spatID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(i)->addExternalFlow(Model::RIVER, head * Model::si::meter, conduct,
                                              bottom * Model::si::meter);

            }
        }

        void readGWRecharge(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->addExternalFlow(Model::RECHARGE,
                                                0 * Model::si::meter,
                                                data * nodes->at(pos)->getProperties().get < Model::quantity <
                                                Model::SquareMeter > ,
                                                Model::Area > ().value(),
                                                0 * Model::si::meter);
            });
        }

    void readInitialHeads(std::string path) {
        readTwoColumns(path, [this](double data, int pos) {
            nodes->at(pos)->setHead_direct(data);
        });
    }

    vector<Model::quantity<Model::Dimensionless>> calcNusInZones(vector<double> densityZones){
        double densityFresh = 1000.0;
        vector<Model::quantity<Model::Dimensionless>> nusInZones;
        for (int id = 0; id < densityZones.size(); id++) {
            // nus of zones is equal to nus of zeta surface below
            nusInZones.push_back(((densityZones[id] - densityFresh) / densityFresh) * Model::si::si_dimensionless);
        }
        return nusInZones;
    }

    vector<Model::quantity<Model::Dimensionless>> calcDelnus(vector<double> densityZones) {
        double densityFresh = 1000.0;
        vector<Model::quantity<Model::Dimensionless>> nusInZones;
        vector<Model::quantity<Model::Dimensionless>> delnus;
        for (int id = 0; id < densityZones.size(); id++) {
            // nus of zones is equal to nus of zeta surface below
            nusInZones.push_back(((densityZones[id] - densityFresh) / densityFresh) * Model::si::si_dimensionless);
            if (id == 0) {
                delnus.push_back(nusInZones[id]); // density difference in top zone
            } else {
                delnus.push_back((nusInZones[id] - nusInZones[id - 1]));
            }
        }
        return delnus;
    }
};
}
}
#endif //TESTING_SIMPLEDATAREADER_HPP
