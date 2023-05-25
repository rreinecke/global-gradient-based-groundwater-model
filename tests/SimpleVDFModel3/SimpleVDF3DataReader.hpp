#ifndef TESTING_SIMPLEVDF3DATAREADER_HPP
#define TESTING_SIMPLEVDF3DATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"
#include "../../src/Model/Node.hpp"
#include "../../src/Model/Units.hpp"

namespace GlobalFlow {
    namespace DataProcessing {

        class SimpleVDF3DataReader : public DataReader {
        public:
            SimpleVDF3DataReader() { }

            virtual void readData(Simulation::Options op) {
                LOG(userinfo) << "Building the initial model layer";
                std::vector<std::vector<int>> grid;
                grid = readGrid(nodes,
                                buildDir(op.getNodesDir()),
                                op.getNumberOfNodesPerLayer(),
                                op.getNumberOfLayers(),
                                op.getNumberOfRows(),
                                op.getNumberOfCols(),
                                op.getInitialK()[0],
                                op.getInitialHead(),
                                op.getAquiferDepth()[0],
                                op.getAnisotropy()[0],
                                op.getSpecificYield(),
                                op.getSpecificStorage(),
                                op.getEdgeLengthLeftRight(),
                                op.getEdgeLengthFrontBack(),
                                op.isConfined(0),
                                op.isDensityVariable(),
                                op.getEffectivePorosity(),
                                op.getMaxTipSlope(),
                                op.getMaxToeSlope(),
                                op.getMinDepthFactor(),
                                op.getSlopeAdjFactor(),
                                op.getVDFLock(),
                                op.getDensityZones());

                LOG(userinfo) << "Building the model layer(s) below";
                DataProcessing::buildBottomLayers(nodes,
                                                  op.getNumberOfLayers(),
                                                  op.getConfinements(),
                                                  op.getAquiferDepth(),
                                                  op.getInitialK(),
                                                  op.getAnisotropy());

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
                    LOG(userinfo) << "Reading zetas";
                    readInitialZetas(op.getNumberOfNodesPerLayer(), op.getNumberOfLayers(),
                                     buildDir(op.getInitialZetasDir()));
                    if (op.isEffectivePorosityFromFile()){
                        LOG(userinfo) << "Reading effective porosity";
                        readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
                    }
                    if (op.isZonesSourcesSinksFromFile()){
                        LOG(userinfo) << "Reading zones of sinks and sources";
                        readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinksDir()), op.getDensityZones());
                    }
                }

                if (op.isInitialHeadFromFile()){
                    LOG(userinfo) << "Initializing head";
                    readInitialHeads((buildDir(op.getInitialHeadsDir())));
                }

                LOG(userinfo) << "Building grid by rows and columns (boundaries need to be specified in with a file)";
                DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodesPerLayer(), op.getNumberOfLayers());
            }

        private:
            template<class T>
            using Matrix = std::vector<std::vector<T>>;

            Matrix<int>
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
                lookupSpatIDtoNodeIDs.reserve(numberOfNodesPerLayer);
                vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
                vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

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
                                                                aquiferDepth,
                                                                anisotropy,
                                                                specificYield,
                                                                specificStorage,
                                                                confined,
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
                        lookupSpatIDtoNodeIDs[spatID].push_back(nodeID + (numberOfNodesPerLayer * layer));
                    }
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

            void readElevation(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setElevation(data * Model::si::meter);
                });
            };

            /**
             * @brief Read in a custom definition for the general head boundary
             * @param path Where to read from
             */
            void readHeadBoundary(std::string path) {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "head", "conduct");
                int spatID{0};
                double head{0};
                double conduct{0};

                while (in.read_row(spatID, head, conduct)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoNodeIDs[spatID][0]; // only at layer 0
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                       head * Model::si::meter,
                                                       conduct,
                                                       head * Model::si::meter);
                }
            }

            void readGWRecharge(std::string path) {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "layer", "recharge");
                int spatID{0};
                int layer{0};
                double recharge{0};

                while (in.read_row(spatID, layer, recharge)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoNodeIDs[spatID][layer]; // at layer
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                       0 * Model::si::meter,
                                                       recharge * nodes->at(nodeID)->getProperties().get<Model::quantity<Model::SquareMeter>,Model::Area>().value(),
                                                       0 * Model::si::meter);
                }
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
                    vector<int> nodeIDs;
                    try {
                        nodeIDs = lookupSpatIDtoNodeIDs[spatID];
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    for (int nodeID : nodeIDs) { // apply to all layers at this spatID
                        nodes->at(nodeID)->setZoneOfSinksAndSources(zoneOfSinks, zoneOfSources, densityZones.size());
                    }
                }
            }

            void readInitialZetas(int numberOfNodesPerLayer, int numberOfLayers, std::string pathZetas) {
                double topOfNode;
                double bottomOfNode;
                int spatID{0};
                double localZetaID{0};
                double zeta{0};

                // add zeta surfaces to top and bottom of each node
                int numberOfNodes = numberOfNodesPerLayer * numberOfLayers;
                for (int nodeIter = 0; nodeIter < numberOfNodes; ++nodeIter){
                    topOfNode = nodes->at(nodeIter)->getProperties().get<Model::quantity<Model::Meter>,Model::Elevation>().value();
                    bottomOfNode = topOfNode - nodes->at(nodeIter)->getProperties().get<Model::quantity<Model::Meter>,Model::VerticalSize>().value();

                    nodes->at(nodeIter)->addZeta(0, topOfNode * Model::si::meter);
                    nodes->at(nodeIter)->addZeta(1, bottomOfNode * Model::si::meter);
                }

                // read initial data for density surfaces
                int nodeID = 0;
                for (int layer = 0; layer < numberOfLayers; layer++) {
                    io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(pathZetas);
                    inZetas.read_header(io::ignore_no_column, "spatID", "localZetaID", "zeta");
                    while (inZetas.read_row(spatID, localZetaID, zeta)) {
                        try {
                            nodeID = lookupSpatIDtoNodeIDs[spatID][layer];
                        }
                        catch (const std::out_of_range &ex) { // if node does not exist ignore entry
                            continue;
                        }
                        nodes->at(nodeID)->addZeta(localZetaID, zeta * Model::si::meter);
                    }
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
#endif //TESTING_SIMPLEVDF3DATAREADER_HPP
