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
                                op.isDensityStratified(),
                                op.getDensityFresh(),
                                op.getDensityZones(),
                                op.getNumberOfDensityZones(),
                                op.getMaxToeSlope(),
                                op.getMaxTipSlope());

                LOG(userinfo) << "Reading hydraulic parameters";
                readConduct(buildDir(op.getLithology()));
                readElevation(buildDir(op.getElevation()));

                LOG(userinfo) << "Reading the groundwater recharge";
                readGWRecharge(buildDir(op.getRecharge()));

                LOG(userinfo) << "Reading the zones of sources and sinks";
                readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinks()));

                LOG(userinfo) << "Reading the boundary condition";
                readHeadBoundary(buildDir(op.getKGHBDir()));

                LOG(userinfo) << "Reading parameters for variable density flow";
                readInitialZetas(op.isDensityVariable(), op.getDensityFresh(),buildDir(op.getInitialZetasDir()), op.getAquiferDepth()[0]);
                readEffectivePorosity(buildDir(op.getEffectivePorosity()));

                LOG(userinfo) << "Initializing head";
                readInitialHeads((buildDir(op.getInitialHeadsDir())));

                LOG(userinfo) << "Connecting the model cells";
                DataProcessing::buildByGrid(nodes, grid, op.getNumberOfLayers(), op.getGHBConduct(),
                                            op.getBoundaryCondition());
            }

        private:
            template<class T>
            using Matrix = std::vector<std::vector<T>>;

            Matrix<int>
            readGrid(NodeVector nodes, std::string path, int numberOfNodes, int numberOfRows, int numberOfCols,
                     double defaultK,
                     double aquiferDepth,
                     double anisotropy,
                     double specificYield,
                     double specificStorage,
                     double edgeLengthLeftRight,
                     double edgeLengthFrontBack,
                     bool confined,
                     bool densityVariable, // todo: not needed?
                     bool densityStratified,
                     double densityFresh,
                     vector<double> densityZones,
                     int numberOfDensityZones,
                     double maxToeSlope,
                     double maxTipSlope) {
                Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

                io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "arcID", "X", "Y", "cell_area", "row", "col");

                double x{0};
                double y{0};
                double area{0};
                int arcID{0};
                int i{0};
                int row{0};
                int col{0};
                lookuparcIDtoID.reserve(numberOfNodes);
                LOG(debug) << "numberOfDensityZones: " << numberOfDensityZones << std::endl;
                Model::DensityProperties densityProperties =
                        Model::DensityProperties::setDensityProperties(densityVariable,
                                                                       densityStratified,
                                                                       densityFresh,
                                                                       densityZones,
                                                                       numberOfDensityZones, maxToeSlope, maxTipSlope);

                while (in.read_row(arcID, x, y, area, row, col)) {
                    out[row][col] = i;
                    nodes->emplace_back(new Model::StandardNode(nodes,
                                                                x,
                                                                y,
                                                                area * Model::si::square_meter,
                                                                edgeLengthLeftRight * Model::si::meter,
                                                                edgeLengthFrontBack * Model::si::meter,
                                                                (unsigned long) arcID,
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
                    lookuparcIDtoID[arcID] = i;
                    i++;
                }

                return out;
            }

            void readConduct(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    if (data > 10) {
                        LOG(debug) << "Very high conductance value at arcID "<<pos<<". Possible Data Error";
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
                in.read_header(io::ignore_no_column, "arcID", "elevation", "conduct");
                int arcid{0};
                double elevation{0};
                double conduct{0};

                while (in.read_row(arcid, elevation, conduct)) {
                    int pos = 0;
                    try {
                        pos = lookuparcIDtoID.at(arcid);
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
                                                    data * nodes->at(pos)->getProperties().get<Model::quantity<
                                                            Model::SquareMeter>,
                                                            Model::Area>().value(),
                                                    0 * Model::si::meter);
                });
            }

            void readZonesSourcesSinks(std::string path) {
                readTwoColumns(path, [this](int data, int pos) {
                    nodes->at(pos)->setZoneOfSourcesAndSinks(data);
                });
            }

            void readInitialZetas(bool densityVariable, double densityFresh, std::string pathZetas, double aquiferDepth) {
                if (densityVariable){
                    int arcid{0};
                    double density{0};
                    double height{0};
                    unsigned long int globalZetaID{0};

                    // read initial data for density surfaces
                    io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(pathZetas);
                    inZetas.read_header(io::ignore_no_column, "arcID", "density", "height");
                    while (inZetas.read_row(arcid, density, height)) {
                        int pos = 0;
                        try {
                            pos = lookuparcIDtoID.at(arcid);
                        }
                        catch (const std::out_of_range &ex) {
                            //if Node does not exist ignore entry
                            continue;
                        }
                        nodes->at(pos)->addInitialZeta(height * Model::si::meter, globalZetaID);
                        globalZetaID++;
                    }
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