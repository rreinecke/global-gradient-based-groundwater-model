#ifndef TESTING_SIMPLEDATAREADER_HPP
#define TESTING_SIMPLEDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"

namespace GlobalFlow {
namespace DataProcessing {


class SimpleDataReader : public DataReader {
    public:
        SimpleDataReader() { }

        virtual void readData(Simulation::Options op) {
            LOG(userinfo) << "Building the initial model layer";
            LOG(userinfo) << "Building the initial model layer";
            readLandMask(nodes,
                         buildDir(op.getNodesDir()),
                         op.getNumberOfNodesPerLayer(),
                         op.getNumberOfLayers(),
                         op.getInitialK()[0],
                         op.getInitialHead(),
                         op.getAquiferDepth()[0],
                         op.getAnisotropy()[0],
                         op.getSpecificYield(),
                         op.getSpecificStorage(),
                         op.useEfolding(),
                        op.isConfined(0),
                         op.isDensityVariable(),
                         op.getEffectivePorosity(),
                         op.getMaxTipSlope(),
                         op.getMaxToeSlope(),
                         op.getMinDepthFactor(),
                         op.getSlopeAdjFactor(),
                         op.getVDFLock(),
                        op.getDensityZones());

            LOG(userinfo) << "Building grid by spatial ID";
            DataProcessing::buildBySpatID(nodes,
                                          this->getMappingSpatIDtoNodeIDs(),
                                          1, // resolution = 0.0833 <- input for global models
                                          10, // lonRange = 360
                                          10, // latRange = 180
                                          false, // isGlobal = true
                                          op.getNumberOfNodesPerLayer(),
                                          op.getGHBConduct(),
                                          op.getBoundaryCondition());

            LOG(userinfo) << "Building the bottom layers";
            DataProcessing::buildBottomLayers(nodes,
                                              op.getNumberOfLayers(),
                                              op.getConfinements(),
                                              op.getAquiferDepth(),
                                              op.getInitialK(),
                                              op.getAnisotropy());

            LOG(userinfo) << "Reading hydraulic parameters";
            readConduct(buildDir(op.getLithology()));
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Initializing head";
            readInitialHeads((buildDir(op.getInitialHeadsDir())));

            LOG(userinfo) << "Defining rivers";
            readRiverConductance(buildDir(op.getKRiver()));

            //LOG(userinfo) << "Building grid by rows and columns (boundaries need to be specified in with a file)";
            //DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodesPerLayer(), op.getNumberOfLayers());
        }

    private:
        template<class T>
        using Matrix = std::vector<std::vector<T>>;

        /*Matrix<int>
        readGrid(NodeVector nodes,
                 std::string path,
                 int numberOfNodesPerLayer,
                 int numberOfLayers,
                 int numberOfRows,
                 int numberOfCols,
                 double defaultK,
                 double initialHead,
                 double aquiferDepth,
                 double anisotropy,
                 double specificYield,
                 double specificStorage,
                 bool useEfolding,
                 double edgeLengthLeftRight,
                 double edgeLengthFrontBack,
                 bool confined,
                 bool isDensityVariable,
                 double effPorosity,
                 double maxTipSlope,
                 double maxToeSlope,
                 double minDepthFactor,
                 double slopeAdjFactor,
                 double vdfLock,
                 std::vector<double> densityZones) {
            Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

            io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "X", "Y", "cell_area", "col", "row");

            double x{0};
            double y{0};
            double area{0};
            int spatID{0};
            int nodeID{0};
            int row{0};
            int col{0};
            int refID{0};
            lookupSpatIDtoNodeIDs.reserve(numberOfNodesPerLayer);
            std::vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
            std::vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

            while (in.read_row(spatID, x, y, area, col, row)) {
                out[col][row] = nodeID;
                nodes->emplace_back(new Model::StandardNode(nodes,
                                                            y,
                                                            x,
                                                            area * Model::si::square_meter,
                                                            edgeLengthLeftRight * Model::si::meter,
                                                            edgeLengthFrontBack * Model::si::meter,
                                                            (large_num) spatID,
                                                            (large_num) nodeID,
                                                            defaultK * (Model::si::meter / Model::day),
                                                            initialHead * Model::si::meter,
                                                            aquiferDepth,
                                                            anisotropy,
                                                            specificYield,
                                                            specificStorage,
                                                            useEfolding,
                                                            confined,
                                                            refID,
                                                            isDensityVariable,
                                                            delnus,
                                                            nusInZones,
                                                            effPorosity,
                                                            maxTipSlope,
                                                            maxToeSlope,
                                                            minDepthFactor,
                                                            slopeAdjFactor,
                                                            vdfLock * Model::si::meter));
                for (int layer = 0; layer < numberOfLayers; layer++) { // todo: not ideal. move to neighbouring?
                    lookupSpatIDtoNodeIDs[spatID][layer][refID] = nodeID + (numberOfNodesPerLayer * layer);
                }
                nodeID++;
            }

            return out;
        }*/

    /**
         * @brief Initial reading of node definitions - without col and row
         * @note Without col and row
         * Reads a csv file with x and y coordinates for predefined grid of cells
         * @param nodes Vector of nodes
         * @param path Path to read definitions from
         * @param numberOfNodesPerLayer The number of expected computation nodes
         * @param defaultK The default conductivity
         * @param aquiferDepth The default depth per cell
         * @param anisotropy The default relation of vertical and horizontal conductivity
         * @param specificYield The default specific yield
         * @param specificStorage The default specific storage
         * @param confined If node is part of a confined layer?
         * @return number of total top nodes
         */
    int
    readLandMask(NodeVector nodes,
                        std::string path,
                        large_num numberOfNodesPerLayer,
                        int numberOfLayers,
                        double defaultK,
                        double initialHead,
                        double aquiferDepth,
                        double anisotropy,
                        double specificYield,
                        double specificStorage,
                        bool useEfolding,
                        bool confined,
                        bool isDensityVariable,
                        double effPorosity,
                        double maxTipSlope,
                        double maxToeSlope,
                        double minDepthFactor,
                        double slopeAdjFactor,
                        double vdfLock,
                        std::vector<double> densityZones) {
        io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
        in.read_header(io::ignore_no_column, "spatID", "X", "Y", "area"); // todo use new refIDs (1, 11, 111 for deeper levels) only for reading
        double X{0};
        double Y{0};
        double area{0};
        int refID{0};
        large_num nodeID{0};
        large_num spatID{0};

        lookupSpatIDtoNodeIDs.reserve(numberOfNodesPerLayer);
        std::vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
        std::vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

        while (in.read_row(spatID, X, Y, area)) {
            //area is in km needs to be in m
            //TODO implement a container for these parameters
            nodes->emplace_back(new Model::StandardNode(nodes,
                                                        Y,
                                                        X,
                                                        area * Model::si::square_meter,
                                                        std::sqrt(area)*Model::si::meter,
                                                        std::sqrt(area)*Model::si::meter,
                                                        spatID, // spatID
                                                        nodeID,
                                                        defaultK * (Model::si::meter / Model::day),
                                                        initialHead * Model::si::meter,
                                                        aquiferDepth,
                                                        anisotropy,
                                                        specificYield,
                                                        specificStorage,
                                                        useEfolding,
                                                        confined,
                                                        refID,
                                                        isDensityVariable,
                                                        delnus,
                                                        nusInZones,
                                                        effPorosity,
                                                        maxTipSlope,
                                                        maxToeSlope,
                                                        minDepthFactor,
                                                        slopeAdjFactor,
                                                        vdfLock * Model::si::meter));
            for (int layer = 0; layer < numberOfLayers; layer++) { // todo: not ideal. move to neighbouring?
                lookupSpatIDtoNodeIDs[spatID][layer][refID] = nodeID + (numberOfNodesPerLayer * layer);
                //LOG(debug) << "spatID: " << spatID << ", nodeID: " << nodeID << "refID = " << refID;
            }
            nodeID++;
        }
        //Return number of total top nodes
        return nodeID - 1;
    };

        void readConduct(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                if (data > 10) {
                    LOG(debug) << "Very high conductance value at spatID "<<pos<<". Possible Data Error";
                }
                nodes->at(pos)->setK(data * (Model::si::meter / Model::day));
            });
        };

        void readElevation(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->setElevation(data * Model::si::meter);
            });
        };

        void readRiverConductance(std::string path) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "Head", "Bottom", "Conduct");
            large_num spatID{0};
            double head{0};
            double conduct{0};
            double bottom{0};
            large_num nodeID{0};

            while (in.read_row(spatID, head, bottom, conduct)) {
                try {
                    nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(0).at(0); //layer = 0, refID = 0
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->addExternalFlow(Model::RIVER, head * Model::si::meter, conduct,
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

    std::vector<Model::quantity<Model::Dimensionless>> calcNusInZones(std::vector<double> densityZones){
        double densityFresh = 1000.0;
        std::vector<Model::quantity<Model::Dimensionless>> nusInZones;
        for (int id = 0; id < densityZones.size(); id++) {
            // nus of zones is equal to nus of zeta surface below
            nusInZones.push_back(((densityZones[id] - densityFresh) / densityFresh) * Model::si::si_dimensionless);
        }
        return nusInZones;
    }

    std::vector<Model::quantity<Model::Dimensionless>> calcDelnus(std::vector<double> densityZones) {
        double densityFresh = 1000.0;
        std::vector<Model::quantity<Model::Dimensionless>> nusInZones;
        std::vector<Model::quantity<Model::Dimensionless>> delnus;
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
