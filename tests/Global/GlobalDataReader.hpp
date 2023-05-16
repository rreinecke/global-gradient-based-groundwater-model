/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "../../src/DataProcessing/DataReader.hpp"
#include "../../src/Model/Node.hpp"
#include "../../src/Model/Units.hpp"

namespace GlobalFlow {
    namespace DataProcessing {
        /**
         * @class GlobalDataReader
         * @implements DataReader
         * @brief This class provides methods for loading large input data
         * The paths are specified in the json file in in the data folder
         */
        class GlobalDataReader : public DataReader {
        public:
            /**
             * @brief Constructor
             * @param step Day, Month, ...
             */
            GlobalDataReader(int step) { stepMod = step; }

            virtual void readData(Simulation::Options op) {
                LOG(userinfo) << "Building the top model layer";
                std::vector<std::vector<int>> grid;
                if (op.isRowCol()) {
                    LOG(userinfo) << "- reading grid by rows and columns";
                    grid = readGrid(nodes,
                                    buildDir(op.getNodesDir()),
                                    op.getNumberOfNodesPerLayer(),
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
                                    op.getEffectivePorosity(),
                                    op.getMaxTipSlope(),
                                    op.getMaxToeSlope(),
                                    op.getMinDepthFactor(),
                                    op.getSlopeAdjFactor(),
                                    op.getVDFLock(),
                                    op.getDensityZones());
                } else {
                    LOG(userinfo) << "- reading land mask";
                    readLandMask(nodes, buildDir(op.getNodesDir()), op.getNumberOfNodesPerLayer(),
                                 op.getInitialK(), op.getInitialHead(),op.getAquiferDepth()[0],
                                 op.getAnisotropy(), op.getSpecificYield(), op.getSpecificStorage(),
                                 op.isConfined(0), op.isDensityVariable(),
                                 op.getEffectivePorosity(), op.getMaxToeSlope(), op.getMaxToeSlope(),
                                 op.getMinDepthFactor(), op.getSlopeAdjFactor(), op.getVDFLock(), op.getDensityZones());
                    LOG(userinfo) << "- reading mapping";
                    readSpatIDtoArcID(buildDir(op.getMapping()));
                }

                LOG(userinfo) << "Building the model layer(s) below";
                DataProcessing::buildBottomLayers(nodes,
                                                  op.getNumberOfLayers(),
                                                  op.getConfinements(),
                                                  op.getAquiferDepth());

                LOG(userinfo) << "Reading groundwater recharge";
                //readGWRecharge(buildDir(op.getRecharge()));
                readGWRechargeMapping(buildDir(op.getRecharge()),
                                      [](const double &recharge, const double &area) {
                                          return (((recharge / 1000) * area) / 365);});

                if (op.isKFromFile()) {
                    LOG(userinfo) << "Reading hydraulic conductivity";
                    readConduct(buildDir(op.getLithology()));
                }

                LOG(userinfo) << "Reading elevation";
                readElevation(buildDir(op.getElevation()));

                LOG(userinfo) << "Reading e-folding";
                readEfold(buildDir(op.getEfolding()), op.getEfolding_a());

                LOG(userinfo) << "Reading slope";
                readSlope(buildDir(op.getSlope()), op.getSlope_a());

                // read either initial head (default) or equal water table depth from file, if available
                if (op.isInitialHeadFromFile()){
                    LOG(userinfo) << "Reading initial head";
                    readInitialHeads((buildDir(op.getInitialHeadsDir())));
                } else if (op.isEqWTDFromFile()){
                    LOG(userinfo) << "Reading equal water table depth";
                    readEqWTD(buildDir(op.getEqWTD())); // requires elevation to be set
                }


                if (op.isKRiverFromFile()) {
                    LOG(userinfo) << "Reading river conductance";
                    readRiverConductance(buildDir(op.getKRiver()));
                } else {
                    LOG(userinfo) << "Reading blue cells";
                    readBlueCells(buildDir(op.getSurfaceWaterElevation()),
                                  calculateRiverStage(buildDir(op.getRiverExtent())));
                }

                LOG(userinfo) << "Reading lakes and wetlands";
                readLakesAndWetlands(buildDir(op.getLocalWetlands()),
                                     buildDir(op.getGlobalWetlands()),
                                     buildDir(op.getGlobalLakes()),
                                     buildDir(op.getGlobalLakes()));

                if (op.isDensityVariable()) {
                    LOG(userinfo) << "Reading initial zeta heights";
                    readInitialZetas(op.getNumberOfNodesPerLayer() * op.getNumberOfLayers(),
                                     buildDir(op.getInitialZetasDir())); // requires elevation to be set

                    if (op.isEffectivePorosityFromFile()) {
                        LOG(userinfo) << "Reading effective porosity";
                        readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
                    }
                    if (op.isZonesSourcesSinksFromFile()) {
                        LOG(userinfo) << "Reading zones of sources and sinks";
                        readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinksDir()), op.getDensityZones().size());
                    }
                }

                if (op.isRowCol()) {
                    LOG(userinfo) << "Building grid by rows and columns (boundaries need to be specified in with a file)";
                    DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodesPerLayer(), op.getNumberOfLayers());
                } else {
                    LOG(userinfo) << "Building grid by spatial ID";
                    //DataProcessing::buildNeighbourMap(nodes, i, op.getNumberOfLayers(), op.getOceanConduct(), op.getBoundaryCondition());
                    DataProcessing::buildBySpatID(nodes, this->getMappingSpatIDtoID(), 60*5 , op.getNumberOfLayers(),
                                                  op.getGHBConduct(), op.getBoundaryCondition());
                }

                if(op.isKGHBFile()) {
                    LOG(userinfo) << "Reading the boundary condition";
                    readHeadBoundary(buildDir(op.getKGHBDir()));
                }

                if (op.isDensityVariable()) {
                    LOG(userinfo) << "Setting the variable density conditions at general head boundaries";
                    // Needs to be called after GHB was set (after buildByGrid/buildBySpatID and readHeadBoundary)
                    setVariableDensityConditionsAtBoundary(op.getNumberOfLayers(),
                                                           op.getDensityZones().size(),
                                                           op.getAquiferDepth()[0]);
                }

                LOG(userinfo) << "Copying neighbours to bottom layer(s)";
                DataProcessing::copyNeighboursToBottomLayers(nodes,op.getNumberOfLayers());
            }

            template<class T>
            using Matrix = std::vector<std::vector<T>>;

            /**
             * @brief Method for already gridded definitions - that is structured in row and column
             * @note Structured in row, col
             * @param nodes Vector of nodes
             * @param path Path to read definitions from
             * @param numberOfNodesPerLayer The number of expected computation nodes
             * @param defaultK The default conductivity
             * @param aquiferDepth The default depth per cell
             * @param anisotropy The default relation of vertical and horizontal conductivity
             * @param specificYield The default specific yield
             * @param specificStorage The default specific storage
             * @param confined If node is part of a confined layer?
             * @return A Matrix of computational nodes
             */
            Matrix<int>
            readGrid(NodeVector nodes,
                     std::string path,
                     int numberOfNodesPerLayer,
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
                in.read_header(io::ignore_no_column, "spatID", "lon", "lat", "area","col","row");

                double lat{0};
                double lon{0};
                double area{0};
                int spatID{0};
                int nodeID{0};
                int row{0};
                int col{0};
                lookupSpatIDtoNodeID.reserve(numberOfNodesPerLayer);
                vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
                vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

                while (in.read_row(spatID, lon, lat, area, col, row)) {
                    out[col][row] = nodeID;
                    //area is in km needs to be in m
                    nodes->emplace_back(new Model::StandardNode(nodes,
                                                                lat,
                                                                lon,
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
                                                                delnus,
                                                                nusInZones,
                                                                effPorosity,
                                                                maxTipSlope,
                                                                maxToeSlope,
                                                                minDepthFactor,
                                                                slopeAdjFactor,
                                                                vdfLock * Model::si::meter));
                    lookupSpatIDtoNodeID[spatID] = nodeID;
                    nodeID++;
                }

                return out;
            };

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
                         int numberOfNodesPerLayer,
                         double defaultK,
                         double initialHead,
                         double aquiferDepth,
                         double anisotropy,
                         double specificYield,
                         double specificStorage,
                         bool confined,
                         bool isDensityVariable,
                         double effPorosity,
                         double maxTipSlope,
                         double maxToeSlope,
                         double minDepthFactor,
                         double slopeAdjFactor,
                         double vdfLock,
                         vector<double> densityZones) {
                io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "lon", "lat", "area");
                double lon{0};
                double lat{0};
                double area{0};
                large_num spatID{0};
                int nodeID{0};
                lookupSpatIDtoNodeID.reserve(numberOfNodesPerLayer);
                vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
                vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

                while (in.read_row(spatID, lon, lat, area)) {
                    //area is in km needs to be in m
                    //TODO implement a container for these parameters
                    nodes->emplace_back(new Model::StandardNode(nodes,
                                                                lat,
                                                                lon,
                                                                1e+6 * area * Model::si::square_meter,
                                                                std::sqrt(1e+6 * area)*Model::si::meter,
                                                                std::sqrt(1e+6 * area)*Model::si::meter,
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
                                                                delnus,
                                                                nusInZones,
                                                                effPorosity,
                                                                maxTipSlope,
                                                                maxToeSlope,
                                                                minDepthFactor,
                                                                slopeAdjFactor,
                                                                vdfLock * Model::si::meter));
                    lookupSpatIDtoNodeID[spatID] = nodeID;
                    nodeID++;
                }
                //Return number of total top nodes
                return nodeID - 1;
            };

            /**
             * @brief Read in a custom definition for the general head boundary
             * @param path Where to read from
             */
            void readHeadBoundary(std::string path) {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "elevation", "conduct");
                int spatID{0};
                double elevation{0};
                double conduct{0};

                while (in.read_row(spatID, elevation, conduct)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoNodeID.at(spatID);
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

            void setVariableDensityConditionsAtBoundary(int layers, int numZetas, double aquiferDepth){
                vector<Model::quantity<Model::Meter>> zetas;
                Model::quantity<Model::Meter> zeta;
                large_num nodes_per_layer = nodes->size() / layers;
                LOG(userinfo) << "numZetas: " << numZetas;
                for (int nodeID = 0; nodeID < nodes_per_layer; ++nodeID) {
                    if (nodes->at(nodeID)->hasGHB()){
                        for (int zetaID = 0; zetaID < numZetas; ++zetaID){
                            double elevation = nodes->at(nodeID)->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>().value();

                            if (zetaID == numZetas-1){
                                zeta = (elevation - aquiferDepth) * Model::si::meter;
                            } else {
                                zeta = elevation * Model::si::meter;
                            }
                            nodes->at(nodeID)->setZeta(zetaID, zeta);
                        }
                        nodes->at(nodeID)->setZoneOfSinksAndSources(0, numZetas-2, numZetas-1);
                    }
                }
            }

            /**
             * @brief Read in a custom definition file for initial heads
             * @param path Where to read the file from
             */
            void readInitialHeads(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setHead_direct(data);
                });
            }

            /**
             * @brief Read in a custom river definition file
             * Structured as: spatID, Head, Bottom, Conduct
             * @param path Where to read the file from
             */
            void readRiverConductance(std::string path) {
                io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "Head", "Bottom", "Conduct");
                int spatID{0};
                double head{0};
                double conduct{0};
                double bottom{0};

                while (in.read_row(spatID, head, bottom, conduct)) {
                    int i = 0;
                    try {
                        i = lookupSpatIDtoNodeID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(i)->addExternalFlow(Model::RIVER, head * Model::si::meter, conduct, bottom * Model::si::meter);
                }
            }

            /**
             * @brief Read e-folding data from a specified path
             * @param path Where to read the file from
             * @param files If different files for different regions are given
             */
            void readEfold(std::string path, std::vector<std::string> files) {
                loopFiles(path, files, [this](std::string path) {
                    readTwoColumns(path, [this](double data, int nodeID) {
                        nodes->at(nodeID)->setEfold(data);
                    });
                });
            };

            /**
             * @brief Read slope data from a specified path
             * @param path Where to read the file from
             * @param files If different files for different regions are given
             */
            void readSlope(std::string path, std::vector<std::string> files) {
                loopFiles(path, files, [this](std::string path) {
                    readTwoColumns(path, [this](double data, int nodeID) {
                        nodes->at(nodeID)->setSlope(data);
                    });
                });
            };

            /**
             * @brief Read elevation data from a specified path
             * Used for multiple files
             * @note !Uses setElevation() function. Should only be called after all layers are build
             * as it affects layers below
             * @param path Where to read the file from
             * @param files If different files for different regions are given
             */
            void readElevation(std::string path, std::vector<std::string> files) {
                loopFiles(path, files, [this](std::string path) {
                    readTwoColumns(path, [this](double data, int nodeID) {
                        nodes->at(nodeID)->setElevation(data * Model::si::meter);
                    });
                });
            };

            /**
             * @brief Read elevation data from a specified path
             * @note !Uses setElevation() function. Should only be called after all layers are build
             * as it affects layers below
             * @param path Where to read the file from
             */
            void readElevation(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setElevation(data * Model::si::meter);
                });
            };

            /**
             * @brief Read equilibrium water-table information used for the dynamic river computation
             * @param path Where to read the file from
             */
            void readEqWTD(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setEqHead(data * Model::si::meter);
                });
            };

            /**
             * @brief Read diffuse gw-recharge
             * @param path Where to read the file from
             */
            void readGWRecharge(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                       0 * Model::si::meter,
                                                       data *
                                                       nodes->at(nodeID)->getProperties().get<Model::quantity<Model::SquareMeter>,
                                                               Model::Area>().value(),
                                                       0 * Model::si::meter);
                });
            }

            /**
            * @brief Read diffuse groundwater recharge from a file and map the value using a conversion function
            * @tparam ConversionFunction Allows the dynamic recalculation of recharge based on cell area
            * @param path Where to read the file from
            * @param convertToRate The conversion function
            */
            template<typename ConversionFunction>
            void readGWRechargeMapping(std::string path, ConversionFunction convertToRate) {
                io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "arcID", "data");

                int arcID = -1;
                std::vector<int> tmp;

                //Count missing mappings due to landmask diff
                //int missingMapping = 0;

                double recharge = 0;

                while (in.read_row(arcID, recharge)) {
                    //lookup nodes which get the special_flow
                    try {
                        tmp = this->lookupArcIDtoSpatIDs.at(arcID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        //missingMapping++;
                        continue;
                    }

                    for (int spatID : tmp) {
                        int i = 0;
                        try {
                            i = this->lookupSpatIDtoNodeID.at(spatID);
                        }
                        catch (const std::out_of_range &ex) {
                            //if Node does not exist ignore entry
                            continue;
                        }
                        /**
                        * GW Recharge is head independent
                        * -> p = 0
                        * Bottom is undefined and c is same as q for recharge flows
                        */
                        //double K = nodes->at(i)->getProperties().get<quantity<Model::Velocity>, Model::K>().value();
                        //double size = nodes->at(
                        //        i)->getProperties().get<quantity<Model::SquareMeter>, Model::SideSurface>().value();
                        //double QMax = K * size;
                        NANChecker(recharge, "Broken recharge value");
                        double q_rech = convertToRate(recharge,
                                                      nodes->at(i)->getProperties().get<Model::quantity<Model::SquareMeter>, Model::Area>().value());
                        //NANChecker(QMax, "QMax Problem");
                        NANChecker(q_rech, "Recharge-init Problem");
                        //LOG(debug) << QMax << "," << q_rech << "q,r";

                        nodes->at(i)->addExternalFlow(Model::RECHARGE,
                                                      0 * Model::si::meter,
                                                      q_rech,
                                                      0 * Model::si::meter);

                    }
                }
            }

            /**
             * @brief Read cell conductance definition
             * @note currently does check if val > 10 m/day
             * @param path Where to read the file from
             */
            void readConduct(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    if (data != 0) {
                        //0 is possible data error, known to occur with Gleeson based map
                        nodes->at(nodeID)->setK(data * (Model::si::meter / Model::day));
                    }
                });
            }

            /**
             * @brief Helper function for @see readBlueCells
             * @param path Where to read the file from
             * @return a map of bankfull depth, stream width, and length
             */
            std::unordered_map<int, std::array<double, 3>> calculateRiverStage(std::string path) {
                io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "POINTID", "length", "bankfull", "spatID", "width");

                int id = 0;
                double lenght = 0;
                double bankfull = 0;
                double width = 0;
                int spatID = 0;

                std::unordered_map<int, std::array<double, 3>> out;

                while (in.read_row(id, lenght, bankfull, spatID, width)) {
                    int i = 0;
                    try {
                        i = lookupSpatIDtoNodeID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    double bankfull_depth = 0.349 * std::pow(bankfull, 0.341);
                    out[i] = {{bankfull_depth, width, lenght * 1000}};
                }
                return out;
            }

            /**
             * @brief Reads in river definitions based on a specific elevation data-set
             * @param file to read from
             * @param bankfull_depth A map with addition information @see calculateRiverStage
             */
            void readBlueCells(std::string file,
                               std::unordered_map<int, std::array<double, 3>> bankfull_depth) {
                io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(file);
                in.read_header(io::ignore_no_column, "spatID", "data");
                double spatID = 0;
                double elevation = 0;

                while (in.read_row(spatID, elevation)) {
                    int i = 0;
                    try {
                        i = lookupSpatIDtoNodeID.at((int) spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        //cout << "ID in elevation that has no corresponding node";
                        continue;
                    }
                    if (nodes->at(i)->getProperties().get<large_num, Model::SpatID>() != spatID) {
                        throw "Error in reading spatID";
                    }

                    double flowHead = elevation;
                    double depth = bankfull_depth[i][0];
                    if (depth <= 1) {
                        depth = 1;
                    }
                    double riverBottom = elevation - depth;
                    double K = nodes->at(i)->getProperties().get<Model::quantity<Model::Velocity>, Model::K>().value();
                    double riverWidthFactor =
                            K * bankfull_depth[i][2] * bankfull_depth[i][1] / (flowHead - riverBottom);
                    if (riverWidthFactor <= 0) { riverWidthFactor = 1; }

                    nodes->at(i)->addExternalFlow(Model::RIVER_MM,
                                                  flowHead * Model::si::meter,
                                                  riverWidthFactor,
                                                  riverBottom * Model::si::meter);
                }
            };

            /**
             * @brief Reads in lakes and wetlands definitions based on @cite Lehner and Döll
             * @param pathGlobalLakes
             * @param pathGlobalWetlands
             * @param pathLocalLakes
             * @param pathLocalWetlands
             */
            void readLakesAndWetlands(std::string pathGlobalLakes,
                                      std::string pathGlobalWetlands,
                                      std::string pathLocalLakes,
                                      std::string pathLocalWetlands) {

                std::vector<std::string> paths = {pathGlobalLakes, pathGlobalWetlands, pathLocalLakes,
                                                  pathLocalWetlands};

                int itter = 0;
                for (std::string path : paths) {
                    io::CSVReader<2, io::trim_chars<'"', '\t'>, io::no_quote_escape<','>> in(path);
                    in.read_header(io::ignore_no_column, "spatID", "data");

                    double percentage = 0;
                    float spatID = 0;
                    while (in.read_row(spatID, percentage)) {
                        int i = 0;
                        try {
                            i = lookupSpatIDtoNodeID.at((int) spatID);
                        }
                        catch (const std::out_of_range &ex) {
                            //if Node does not exist ignore entry
                            continue;
                        }
                        if (nodes->at(i)->getProperties().get<large_num, Model::SpatID>() != (int) spatID) {
                            throw "Error in reading spatID";
                        }

                        if (percentage == 0) {
                            continue;
                        }

                        percentage = percentage / 100;

                        double elevation = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>().value();
                        try {
                            elevation = nodes->at(i)->getExternalFlowByName(Model::RIVER_MM).getFlowHead().value();
                        } catch (const std::out_of_range &ex) {}

                        double flowHead = elevation;
                        double bottom = elevation;
                        double K_s = nodes->at(i)->getProperties().get<Model::quantity<Model::Velocity>, Model::K>().value();
                        if (itter == 1) {
                            //global wetlands
                            percentage = percentage * 0.8;
                        }
                        double A_s =
                                nodes->at(i)->getProperties().get<Model::quantity<Model::SquareMeter>, Model::Area>().value() *
                                percentage;
                        double M = 5;
                        //Simple_ Criv=KLW/M, M is the thickness of the riverbed and K is the hydraulic conductivity of the riverbed
                        double conduct = (K_s * A_s) / M;
                        //if(conduct > 1e6){
                        //    std::cout << "To high K:" << K_s << "area: " << A_s << "at: " << spatID <<"\n";
                        //}

                        if (itter == 0) {
                            //nodes->at(i)->removeExternalFlow(Model::RIVER_MM);
                            //Global LAKE
                            //flowHead -= 10;
                            bottom -= 100;
                            nodes->at(i)->addExternalFlow(Model::LAKE,
                                                          flowHead * Model::si::meter,
                                                          conduct,
                                                          bottom * Model::si::meter);
                        } else if (itter == 1) {
                            //Global WETLANDS
                            bottom -= 2;
                            nodes->at(i)->addExternalFlow(Model::GLOBAL_WETLAND,
                                                          flowHead * Model::si::meter,
                                                          conduct,
                                                          bottom * Model::si::meter);
                        } else if (itter == 2) {
                            //Local LAKE
                            //flowHead -= 10;
                            bottom -= 10;
                            nodes->at(i)->addExternalFlow(Model::LAKE,
                                                          flowHead * Model::si::meter,
                                                          conduct,
                                                          bottom * Model::si::meter);
                        } else {
                            //Local WETLANDS
                            bottom -= 2;
                            nodes->at(i)->addExternalFlow(Model::WETLAND,
                                                          flowHead * Model::si::meter,
                                                          conduct,
                                                          bottom * Model::si::meter);
                        }
                    }
                    itter++;
                }
            }

            /**
             * @brief Adds an evapotranspiration module to cells
             * @deprecated
             */
            void addEvapo() {
                for (const std::unique_ptr<Model::NodeInterface> &node : *nodes.get()) {
                    if (node->getProperties().get<int, Model::Layer>() == 0 and
                        node->getProperties().get<Model::quantity<Model::Meter>,
                                Model::Elevation>().value() <= 600) {
                        node->addExternalFlow(Model::EVAPOTRANSPIRATION,
                                              node->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>(),
                                              0.001 *
                                              node->getProperties().get<Model::quantity<Model::SquareMeter>, Model::Area>().value(),
                                              0.5 * Model::si::meter);
                    }
                }
            }

            void readInitialZetas(int numberOfNodes, std::string pathZetas) {
                double topOfNode;
                double bottomOfNode;
                int spatID{0};
                double localZetaID{0};
                double zeta{0};

                // add zeta surfaces to top and bottom of each node
                for (int nodeIter = 0; nodeIter < numberOfNodes; ++nodeIter){
                    topOfNode = nodes->at(nodeIter)->getProperties().get<Model::quantity<Model::Meter>,Model::Elevation>().value();
                    bottomOfNode = topOfNode - nodes->at(nodeIter)->getProperties().get<Model::quantity<Model::Meter>,Model::VerticalSize>().value();

                    nodes->at(nodeIter)->addZeta(0, topOfNode * Model::si::meter);
                    nodes->at(nodeIter)->addZeta(1, bottomOfNode * Model::si::meter);
                }

                // read initial data for density surfaces
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(pathZetas);
                inZetas.read_header(io::ignore_no_column, "spatID", "localZetaID", "zeta");
                while (inZetas.read_row(spatID, localZetaID, zeta)) {
                    int nodeID = 0;
                    try {
                        nodeID = lookupSpatIDtoNodeID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) { // if node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->addZeta(localZetaID, zeta * Model::si::meter);
                }
            }

            void readEffectivePorosity(std::string path) {
                readTwoColumns(path, [this](double data, int nodeID) {
                    nodes->at(nodeID)->setEffectivePorosity(data * Model::si::si_dimensionless);
                });
            }

            void readZonesSourcesSinks(std::string path, int numDensityZones) {
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
                        nodeID = lookupSpatIDtoNodeID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->setZoneOfSinksAndSources(zoneOfSinks, zoneOfSources, numDensityZones);
                }
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

            vector<Model::quantity<Model::Dimensionless>> calcDelnus(vector<double> densityZones){
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