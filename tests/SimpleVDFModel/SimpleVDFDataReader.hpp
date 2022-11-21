#ifndef TESTING_SIMPLEVDFDATAREADER_HPP
#define TESTING_SIMPLEVDFDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"
#include "../../src/Model/Node.hpp"
#include "../../src/Model/Units.hpp"

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
                                op.getNumberOfRows(),
                                op.getNumberOfCols(),
                                op.getInitialK(),
                                op.getAquiferDepth()[0],
                                op.getAnisotropy(),
                                op.getSpecificYield(),
                                op.getSpecificStorage(),
                                op.getEdgeLengthLeftRight(),
                                op.getEdgeLengthFrontBack(),
                                op.isConfined(0),
                                op.isDensityVariable(),
                                op.getDensityZones(),
                                op.getMaxToeSlope(),
                                op.getMaxTipSlope());

                LOG(userinfo) << "Reading hydraulic parameters";
                readConduct(buildDir(op.getLithology()));
                readElevation(buildDir(op.getElevation()));

                LOG(userinfo) << "Reading the groundwater recharge";
                readGWRecharge(buildDir(op.getRecharge()));

                LOG(userinfo) << "Reading the boundary condition";
                readHeadBoundary(buildDir(op.getKGHBDir()));

                LOG(userinfo) << "Reading parameters for variable density flow";
                readInitialZetas(buildDir(op.getInitialZetasDir()));
                readEffectivePorosity(buildDir(op.getEffectivePorosity()));
                readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinks()), op.getDensityZones());

                LOG(userinfo) << "Initializing head";
                readInitialHeads((buildDir(op.getInitialHeadsDir())));

                LOG(userinfo) << "Connecting the model cells";
                DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodes(), op.getNumberOfLayers(), op.getGHBConduct(),
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
                     double aquiferDepth,
                     double anisotropy,
                     double specificYield,
                     double specificStorage,
                     double edgeLengthLeftRight,
                     double edgeLengthFrontBack,
                     bool confined,
                     bool densityVariable,
                     vector<double> densityZones,
                     double maxToeSlope,
                     double maxTipSlope) {
                Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

                io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "X", "Y", "cell_area", "col", "row");

                double x{0};
                double y{0};
                double area{0};
                int spatID{0};
                int i{0};
                int row{0};
                int col{0};
                lookupGlobalIDtoID.reserve(numberOfNodes);
                Model::DensityProperties densityProperties =
                        Model::DensityProperties::setDensityProperties(densityVariable,
                                                                       densityZones,maxToeSlope, maxTipSlope);

                while (in.read_row(spatID, x, y, area, col, row)) {
                    out[col][row] = i;
                    nodes->emplace_back(new Model::StandardNode(nodes,
                                                                x,
                                                                y,
                                                                area * Model::si::square_meter,
                                                                edgeLengthLeftRight * Model::si::meter,
                                                                edgeLengthFrontBack * Model::si::meter,
                                                                (unsigned long) spatID,
                                                                i,
                                                                defaultK * (Model::si::meter / Model::day),
                                                                stepMod,
                                                                aquiferDepth,
                                                                anisotropy,
                                                                specificYield,
                                                                specificStorage,
                                                                confined,
                                                                densityProperties
                                                                ));
                    lookupGlobalIDtoID[spatID] = i;
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

            void readHeadBoundary(std::string path) {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "elevation", "conduct");
                int spatID{0};
                double elevation{0};
                double conduct{0};

                while (in.read_row(spatID, elevation, conduct)) {
                    int pos = 0;
                    try {
                        pos = lookupGlobalIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                    elevation * Model::si::meter,
                                                    conduct,
                                                    elevation * Model::si::meter);
                }
            }

            void readGWRecharge(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    nodes->at(pos)->addExternalFlow(Model::RECHARGE,
                                                    0 * Model::si::meter,
                                                    data * nodes->at(pos)->getProperties().get<Model::quantity<Model::SquareMeter>,Model::Area>().value(),
                                                    0 * Model::si::meter);
                });
            }

            void readZonesSourcesSinks(std::string path, vector<double> densityZones) {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "zoneOfSinks", "zoneOfSources");
                int spatID{0};
                double zoneOfSinks{0};
                double zoneOfSources{0};

                while (in.read_row(spatID, zoneOfSinks, zoneOfSources)) {
                    int pos = 0;
                    try {
                        pos = lookupGlobalIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(pos)->setZoneOfSinksAndSources(zoneOfSinks, zoneOfSources, densityZones.size());
                }
            }

            void readInitialZetas(std::string pathZetas) {
                int spatID{0};
                double height{0};

                // read initial data for density surfaces
                io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(pathZetas);
                inZetas.read_header(io::ignore_no_column, "spatID", "height");
                while (inZetas.read_row(spatID, height)) {
                    int pos = 0;
                    try {
                        pos = lookupGlobalIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(pos)->addInitialZeta(height * Model::si::meter);
                }
            }

            void readEffectivePorosity(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    nodes->at(pos)->setEffectivePorosity(data * Model::si::si_dimensionless);
                });
            }

            void readInitialHeads(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    nodes->at(pos)->setHead_direct(data);
                });
            }



        };
    }
}
#endif //TESTING_SIMPLEVDFDATAREADER_HPP