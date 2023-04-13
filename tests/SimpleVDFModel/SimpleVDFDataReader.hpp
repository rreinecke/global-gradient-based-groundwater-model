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
                                op.getMaxTipToeSlope(),
                                op.getMinDepthFactor(),
                                op.getSlopeAdjFactor(),
                                op.getVDFLock(),
                                op.getDensityZones());

                LOG(userinfo) << "Reading elevation";
                readElevation(buildDir(op.getElevation()));

                if (op.isKFromFile()) {
                    LOG(userinfo) << "Reading hydraulic conductivity";
                    readConduct(buildDir(op.getLithology()));
                }

                LOG(userinfo) << "Reading the groundwater recharge";
                readGWRecharge(buildDir(op.getRecharge()));

                LOG(userinfo) << "Reading the boundary condition";
                readHeadBoundary(buildDir(op.getKGHBDir()));

                if (op.isDensityVariable()) {
                    LOG(userinfo) << "Reading input for variable density flow";
                    if (op.isInitialZetasFromFile()){
                        readInitialZetas(buildDir(op.getInitialZetasDir()));
                    }
                    if (op.isEffectivePorosityFromFile()){
                        readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
                    }
                    if (op.isZonesSourcesSinksFromFile()){
                        readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinksDir()), op.getDensityZones());
                    }
                }

                if (op.isInitialHeadFromFile()){
                    LOG(userinfo) << "Initializing head";
                    readInitialHeads((buildDir(op.getInitialHeadsDir())));
                }

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
                     double maxTipToeSlope,
                     double minDepthFactor,
                     double slopeAdjFactor,
                     double vdfLock,
                     vector<double> densityZones) {
                Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

                io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "X", "Y", "area","col","row");

                double x{0};
                double y{0};
                double area{0};
                int spatID{0};
                int nodeID{0};
                int row{0};
                int col{0};
                lookupSpatIDtoID.reserve(numberOfNodes);
                vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
                vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);
                vector<Model::quantity<Model::Meter>> initialZetasDim;
                for (int i = 0; i < initialZetas.size(); i++) {
                    initialZetasDim.push_back(initialZetas[i] * Model::si::meter);
                }
                while (in.read_row(spatID, x, y, area, col, row)) {
                    out[col][row] = nodeID;
                    //area is in km needs to be in m
                    nodes->emplace_back(new Model::StandardNode(nodes,
                                                                x,
                                                                y,
                                                                area * Model::si::square_meter,
                                                                edgeLengthLeftRight * Model::si::meter,
                                                                edgeLengthFrontBack * Model::si::meter,
                                                                (unsigned long) spatID,
                                                                nodeID,
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
                                                                maxTipToeSlope,
                                                                minDepthFactor,
                                                                slopeAdjFactor,
                                                                vdfLock * Model::si::meter));
                    lookupSpatIDtoID[spatID] = nodeID;
                    nodeID++;
                }

                return out;
            };

            void readConduct(std::string path) {

                readTwoColumns(path, [this](double data, int nodeID) {
                    if (data > 10) {
                        LOG(debug) << "Very high conductance value at spatID "<<nodeID<<". Possible Data Error";
                    }
                    nodes->at(nodeID)->setK(data * (Model::si::meter / Model::day));
                });
            };

            void
            readElevation(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setElevation(data * Model::si::meter);
                });
            };

            void readHeadBoundary(std::string path) {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "elevation", "conduct");
                int spatID{0};
                double elevation{0};
                double conduct{0};

                while (in.read_row(spatID, elevation, conduct)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                    elevation * Model::si::meter,
                                                    conduct,
                                                    elevation * Model::si::meter);
                }
            }

            void readGWRecharge(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                    0 * Model::si::meter,
                                                    data * nodes->at(nodeID)->getProperties().get<Model::quantity<Model::SquareMeter>,Model::Area>().value(),
                                                    0 * Model::si::meter);
                });
            }

            void readInitialHeads(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setHead_direct(data);
                });
            }

            void readZonesSourcesSinks(std::string path, vector<double> densityZones) {
                /**
                 * Here we use zoneOfSinks and zoneOfSources (containing values between 0 and number of density zones).
                 * Thus, sources and sinks are associated to the respective zone. Rule: zoneOfSinks <= zoneOfSources
                 * For simulation of submarine groundwater discharge:
                 * - zoneOfSources: an integer between 1 and the number of zones (brackish/saline water)
                 * - zoneOfSinks: 0 (fresh water)
                 */

                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "zoneOfSinks", "zoneOfSources");
                int spatID{0};
                double zoneOfSinks{0};
                double zoneOfSources{0};

                while (in.read_row(spatID, zoneOfSinks, zoneOfSources)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->setZoneOfSinksAndSources(zoneOfSinks, zoneOfSources, densityZones.size());
                }
            }

            void readInitialZetas(std::string pathZetas) {
                int spatID{0};
                double zeta0{0};
                double zeta1{0};
                double zeta2{0};

                // read initial data for density surfaces
                io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(pathZetas);
                inZetas.read_header(io::ignore_no_column, "spatID", "zeta0", "zeta1", "zeta2");
                while (inZetas.read_row(spatID, zeta0, zeta1, zeta2)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) { // if node does not exist ignore entry
                        continue;
                    }

                    vector<Model::quantity<Model::Meter>> initialZetas{zeta0 * Model::si::meter,
                                                                       zeta1 * Model::si::meter,
                                                                       zeta2 * Model::si::meter};

                    nodes->at(nodeID)->setInitialZetas(initialZetas);
                }
            }

            void readEffectivePorosity(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setEffectivePorosity(data * Model::si::si_dimensionless);
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
                    nusInZones.push_back(
                            ((densityZones[id] - densityFresh) / densityFresh) * Model::si::si_dimensionless);
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
#endif //TESTING_SIMPLEVDFDATAREADER_HPP
