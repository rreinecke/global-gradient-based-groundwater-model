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
                                op.isConfined(0));

                LOG(userinfo) << "Reading hydraulic parameters";
                readConduct(buildDir(op.getLithology()));
                readElevation(buildDir(op.getElevation()));

                LOG(userinfo) << "Reading the groundwater recharge";
                readGWRecharge(buildDir(op.getRecharge()));

                //LOG(userinfo) << "Reading the boundary condition";
                //readHeadBoundary(buildDir(op.getKGHBDir()));

                //LOG(userinfo) << "Reading parameters for variable density flow";
                //readInitialZetas(buildDir(op.getInitialZetasDir()));
                //readEffectivePorosity(buildDir(op.getEffectivePorosity()));
                //readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinks()), op.getDensityZones());

                //LOG(userinfo) << "Initializing head";
                //readInitialHeads((buildDir(op.getInitialHeadsDir())));

                //LOG(userinfo) << "Defining rivers";
                //readRiver(buildDir(op.getKRiverDir()));

                //LOG(userinfo) << "Reading efolding";
                //readEfold(buildDir(op.getEfolding()));

                LOG(userinfo) << "Reading equal water table depth";
                readEqWTD(buildDir(op.getEqWTD()));

                LOG(userinfo) << "Connecting the model cells";
                DataProcessing::buildByGrid(nodes, grid, op.getNumberOfLayers(), op.getGHBConduct(),
                                            op.getBoundaryCondition());
            }

            template<class T>
            using Matrix = std::vector<std::vector<T>>;

            /**
             * @brief Method for already gridded defintions - that is structered in row and column
             * @note Structured in row, col
             * @param nodes Vector of nodes
             * @param path Path to read definitions from
             * @param numberOfNodes The number of expected computation nodes
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
                     bool confined) {
                Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

                io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);

                in.read_header(io::ignore_no_column, "spatID", "X", "Y", "area","col","row");

                double x{0};
                double y{0};
                double area{0};
                int spatID{0};
                int i{0};
                int row{0};
                int col{0};
                lookupGlobalIDtoID.reserve(numberOfNodes);

                Model::DensityProperties densityProperties = Model::DensityProperties::setDensityProperties();
                LOG(debug) << "path: " << path << std::endl;

                while (in.read_row(spatID, x, y, area, col, row)) {
                    out[col][row] = i;
                    //area is in km needs to be in m
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
                                                                densityProperties));
                    lookupGlobalIDtoID[spatID] = i;
                    i++;
                }

                return out;
            }

            /**
             * @brief Initial reading of node definitions - without col and row
             * @note Without col and row
             * Reads a csv file with x and y coordinates for predefined grid of cells
             * @param nodes Vector of nodes
             * @param path Path to read definitions from
             * @param numberOfNodes The number of expected computation nodes
             * @param defaultK The default conductivity
             * @param aquiferDepth The default depth per cell
             * @param anisotropy The default relation of vertical and horizontal conductivity
             * @param specificYield The default specific yield
             * @param specificStorage The default specific storage
             * @param confined If node is part of a confined layer?
             * @return
             */
            int
            readLandMask(NodeVector nodes, std::string path, int numberOfNodes, double defaultK, double aquiferDepth,
                         double anisotropy,
                         double specificYield,
                         double specificStorage,
                         double edgeLengthLeftRight,
                         double edgeLengthFrontBack,
                         bool confined) {
                io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "X", "Y", "cell_area");
                double x = 0;
                double y = 0;
                double area = 0;
                int spatID = 0;
                int i = 0;
                lookupGlobalIDtoID.reserve(numberOfNodes);
                Model::DensityProperties densityProperties = Model::DensityProperties::setDensityProperties();
                while (in.read_row(spatID, x, y, area)) {
                    //area is in km needs to be in m
                    //TODO implement a container for these parameters
                    nodes->emplace_back(new Model::StandardNode(nodes,
                                                                x,
                                                                y,
                                                                1e+6 * area * Model::si::square_meter,
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
                //Return number of total top nodes
                return i - 1;
            };

            /**
             * @brief Read in a custom defintion for the ocean boundary
             * @param path Where to read from
             */
            void readHeadBoundary(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    if (nodes->at(pos)->hasGHB()) {
                        auto flow = nodes->at(pos)->getExternalFlowByName(Model::GENERAL_HEAD_BOUNDARY);
                        nodes->at(pos)->removeExternalFlow(Model::GENERAL_HEAD_BOUNDARY);
                        nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                        0 * Model::si::meter,
                                                        data,
                                                        0 * Model::si::meter);
                    }
                });
            }

            /**
             * @brief Read in a custom river defintion file
             * Structured as: global_ID, Head, Bottom, Conduct
             * @param path Where to read the file from
             */
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
                        i = lookupGlobalIDtoID.at(spatID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(i)->addExternalFlow(Model::RIVER, head * Model::si::meter, conduct, bottom * Model::si::meter);
                }
            }


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
                    readTwoColumns(path, [this](double data, int pos) {
                        nodes->at(pos)->setElevation(data * Model::si::meter);
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
                readTwoColumns(path, [this](double data, int pos) {
                    nodes->at(pos)->setElevation(data * Model::si::meter);
                });
            };

            /**
             * @brief Read equilibrium water-table information used for the dynamic river computation
             * @param path Where to read the file from
             */
            void readEqWTD(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    nodes->at(pos)->setEqHead(data * Model::si::meter);
                });
            };

            /**
             * @brief Read difuse gw-recharge
             * @param path Where to read the file from
             */
            void readGWRecharge(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    nodes->at(pos)->addExternalFlow(Model::RECHARGE,
                                                    0 * Model::si::meter,
                                                    data *
                                                    nodes->at(pos)->getProperties().get<Model::quantity<Model::SquareMeter>,
                                                            Model::Area>().value(),
                                                    0 * Model::si::meter);
                });
            }

            /**
             * @brief Read cell conductance defintion
             * @note currently does check if val > 10 m/day
             * @param path Where to read the file from
             */
            void readConduct(std::string path) {
                readTwoColumns(path, [this](double data, int pos) {
                    if (data != 0) {
                        //0 is possible data error, known to occur with Gleeson based map
                        nodes->at(pos)->setK(data * (Model::si::meter / Model::day));
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
                in.read_header(io::ignore_no_column, "POINTID", "Lenght", "Bankfull", "spatID", "Width");

                int id = 0;
                double lenght = 0;
                double bankfull = 0;
                double width = 0;
                int spatID = 0;

                std::unordered_map<int, std::array<double, 3>> out;

                while (in.read_row(id, lenght, bankfull, spatID, width)) {
                    int i = 0;
                    try {
                        i = lookupGlobalIDtoID.at(spatID);
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
             * @brief Reads in river defintions based on a specific elevation data-set
             * @param file to read from
             * @param bankfull_depth A map with addition information @see calulateRiverStage
             */
            void readBlueCells(std::string file,
                               std::unordered_map<int, std::array<double, 3>> bankfull_depth) {
                LOG(debug) << "reading from" << file;
                io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(file);
                in.read_header(io::ignore_no_column, "spatID", "data");
                double spatID = 0;
                double elevation = 0;

                while (in.read_row(spatID, elevation)) {
                    int i = 0;
                    try {
                        i = lookupGlobalIDtoID.at((int) spatID);
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
             * @param pathLokalLakes
             * @param pathLokalWetlands
             */
            void readLakesandWetlands(std::string pathGlobalLakes,
                                      std::string pathGlobalWetlands,
                                      std::string pathLokalLakes,
                                      std::string pathLokalWetlands) {

                std::vector<std::string> paths = {pathGlobalLakes, pathGlobalWetlands, pathLokalLakes,
                                                  pathLokalWetlands};

                int itter = 0;
                for (std::string path : paths) {
                    io::CSVReader<2, io::trim_chars<'"', '\t'>, io::no_quote_escape<','>> in(path);
                    in.read_header(io::ignore_no_column, "spatID", "Value");

                    double percentage = 0;
                    float spatID = 0;
                    while (in.read_row(spatID, percentage)) {
                        int i = 0;
                        try {
                            i = lookupGlobalIDtoID.at((int) spatID);
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

                        double elevation = nodes->at(
                                i)->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>().value();
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
                            //Global lakes
                            //flowHead -= 10;
                            bottom -= 100;
                            nodes->at(i)->addExternalFlow(Model::LAKE,
                                                          flowHead * Model::si::meter,
                                                          conduct,
                                                          bottom * Model::si::meter);
                        } else if (itter == 1) {
                            //GLOBAL WETLANDS
                            bottom -= 2;
                            nodes->at(i)->addExternalFlow(Model::GLOBAL_WETLAND,
                                                          flowHead * Model::si::meter,
                                                          conduct,
                                                          bottom * Model::si::meter);
                        } else if (itter == 2) {
                            //Lokal lakes
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

            /**
             * @brief A drainage component similar to @cite de Graaf 2014
             */
            void addDrainageHack() {
                //iterate over all top nodes without a river above 300m elevation
                for (const std::unique_ptr<Model::NodeInterface> &node : *nodes.get()) {
                    if (node->getProperties().get<int, Model::Layer>() == 0) {
                        //Add negative special_flow
                        node->addExternalFlow(Model::FLOODPLAIN_DRAIN,
                                              node->getProperties().get<Model::quantity<Model::Meter>,
                                                      Model::Elevation>(),
                                              node->getK().value() *
                                              node->getProperties().get<Model::quantity<Model::Meter>,
                                                      Model::VerticalSize>().value(),
                                              node->getProperties().get<Model::quantity<Model::Meter>,
                                                      Model::EdgeLengthLeftRight>()); // todo adapt to rectangular nodes
                    }
                }
            }
        };

    }
}
