#ifndef TESTING_SIMPLEVDFDATAREADER_HPP
#define TESTING_SIMPLEVDFDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"

namespace GlobalFlow {
namespace DataProcessing {


class SimpleVDFDataReader : public DataReader {
    public:
        SimpleVDFDataReader(int step) { stepMod = step; }

        virtual void readData(Simulation::Options op) {
            LOG(userinfo) << "Building the initial model layer";
            std::vector<std::vector<int>> grid;
            grid = readGrid(nodes,
                            buildDir(op.getNodesDir()),
                            op.getNumberOfNodes(),
                            op.getInitialK(),
                            op.getAquiferDepth()[0],
                            op.getAnisotropy(),
                            op.getSpecificYield(),
                            op.getSpecificStorage(),
                            op.isConfined(0));

            LOG(userinfo) << "Reading hydraulic parameters";
            readConduct(buildDir(op.getLithology()));
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Connecting the model cells";
            DataProcessing::buildByGrid(nodes, grid, op.getNumberOfLayers(), op.getOceanConduct(),
                                        op.getBoundaryCondition());


        }

    private:
        template<class T>
        using Matrix = std::vector<std::vector<T>>;

        Matrix<int>
        readGrid(NodeVector nodes, std::string path, int numberOfNodes, double defaultK, double aquiferDepth,
                 double anisotropy,
                 double specificYield,
                 double specificStorage, bool confined) {
            Matrix<int> out = Matrix<int>(numberOfNodes, std::vector<int>(numberOfNodes));

            io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "global_ID", "X", "Y", "cell_area", "row", "col");

            double x{0};
            double y{0};
            double area{0};
            int globid{0};
            int i{0};
            int row{0};
            int col{0};
            lookupglobIDtoID.reserve(numberOfNodes);

            while (in.read_row(globid, x, y, area, row, col)) {
                out[row][col] = i;
                //area is in km needs to be in m
                nodes->emplace_back(new Model::StandardNode(nodes,
                                                            x,
                                                            y,
                                                            area * Model::si::square_meter,
                                                            (unsigned long) globid,
                                                            i,
                                                            defaultK * (Model::si::meter / Model::day),
                                                            stepMod,
                                                            aquiferDepth,
                                                            anisotropy,
                                                            specificYield,
                                                            specificStorage, confined));
                lookupglobIDtoID[globid] = i;
                i++;
            }

            return out;
        }

        void readConduct(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                if (data > 10) {
                    LOG(debug) << "Very high conductance value at global_ID "<<pos<<". Possible Data Error";
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
            in.read_header(io::ignore_no_column, "global_ID", "Head", "Bottom", "Conduct");
            int arcid{0};
            double head{0};
            double conduct{0};
            double bottom{0};

            while (in.read_row(arcid, head, bottom, conduct)) {
                int i = 0;
                try {
                    i = lookupglobIDtoID.at(arcid);
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

        void readHeadBoundary(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                data * Model::si::meter,
                                                data * nodes->at(pos)->getProperties().get < Model::quantity <
                                                        Model::SquareMeter > ,
                                                        Model::Area > ().value(),
                                                0 * Model::si::meter);
            });
        }
};
}
}
#endif //TESTING_SIMPLEVDFDATAREADER_HPP

