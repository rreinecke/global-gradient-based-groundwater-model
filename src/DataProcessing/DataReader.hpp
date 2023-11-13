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

#ifndef GLOBAL_FLOW_DATAREADER_HPP
#define GLOBAL_FLOW_DATAREADER_HPP

#include "../../lib/csv.h"
#include "../Model/Node.hpp"
#include "../Logging/Logging.hpp"
#include "../Simulation/Options.hpp"
#include <boost/filesystem.hpp>
#include <fstream>
#include <vector>


namespace GlobalFlow {
    using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<Model::NodeInterface>>>;
    namespace fs = boost::filesystem;
    using large_num = unsigned long int;

    /**
     * @class DataReader
     * @brief Interface that needs to be implemented for reading in required data for the model
     */
    class DataReader {
    protected:
        NodeVector nodes;
        /**
         * @var basePath the standard path for data to be read from
         * relative to executable or absolute path
         */
        std::string basePath{"data"};
        fs::path data_dir{basePath};
        /** @var lookupSpatIDtoNodeIDs <SpatID, Layer, RefID, NodeID>*/
        std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<large_num, large_num>>> lookupSpatIDtoNodeIDs;
        /** @var lookupArcIDtoSpatIDs <ArcID(0.5°), vector<SpatID(5')>>*/
        std::unordered_map<large_num, std::vector<large_num>> lookupArcIDtoSpatIDs;
    public:
        /** Virt destructor -> interface*/
        virtual ~DataReader() = default;

        /**
         * @brief Initialize internal ref to node vector
         * @param nodes The vector of nodes
         */
        void initNodes(NodeVector nodeVector) { this->nodes = nodeVector; }

        /**
         * @brief Entry point for reading simulation data
         * @attention This method needs to be implemented!
         * @note readData() is called by simulation at startup
         * @param op Options object
         */
        virtual void readData(Simulation::Options op) = 0;

        /**
         * @brief Generic method for looping through files inside a directory and applying a generic function
         * @param path the directory
         * @param files a vector of files
         * @param fun a function that is applied e.g. reading the data
         */
        template<class Fun>
        void loopFiles(std::string path, std::vector<std::string> files, Fun fun) {
            for (std::string file : files) {
                std::string real_path = path + '/' + file;
                fun(real_path);
            }
        }

        /**
         * @brief Generic method for looping through files inside a directory and layers, applying a generic function
         * @param path the directory
         * @param files a vector of files
         * @param numberOfLayers number of layers
         * @param fun a function that is applied e.g. reading the data
         */
        template<class Fun>
        void loopFilesAndLayers(std::string path, std::vector<std::string> files, int numberOfLayers, Fun fun) {
            for (std::string file : files) {
                std::string real_path = path + '/' + file;
                fun(real_path, numberOfLayers);
            }
        }

        /**
         * @brief Read data from a two-column csv file and apply function to data
         * @param path to the csv file
         * @param processData A processing function e.g. upscaling of data
         */
        template<class ProcessDataFunction>
        void readTwoColumns(std::string path, ProcessDataFunction processData) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "data");
            large_num spatID{0};
            int layer{0};
            large_num refID{0};
            double data = 0;
            std::unordered_map<large_num,large_num> nodeIDs;
            while (in.read_row(spatID, data)) {
                try {
                    nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                for (auto nodeID : nodeIDs){
                    if (nodes->at(nodeID.second)->getProperties().get<large_num, Model::SpatID>() != spatID) {
                        throw "Error in reading spatID";
                    }
                    processData(data, nodeID.second);
                }
            }
        }

        /**
         * @brief Creates a mapping of 0.5° SpatID to a list of contained 5' ArcIDs
         * @param path to file
         */
        void readSpatIDtoArcID(std::string path) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "arcID");
            large_num arcID = 0;
            large_num spatID = 0;

            while (in.read_row(spatID, arcID)) {
                lookupArcIDtoSpatIDs[arcID].push_back(std::move(spatID));
            }
        }


        /**
         * @brief provides acccess to mapping of different resolutions
         * @return <ArcID(0.5°), vector<SpatID(5')>>
         */
        const std::unordered_map<large_num, std::vector<large_num>> &getMappingArcIDtoSpatIDs() {
            return lookupArcIDtoSpatIDs;
        };

        /**
         * @brief provides access to mapping of data ids to position in node vector
         * @return <SpatID, vector<NodeID> (internal array id)>
         */
        void addMappingSpatIDtoNodeIDs(large_num spatID, int layer, large_num refID, large_num nodeID) {
            lookupSpatIDtoNodeIDs[spatID][layer][refID] = nodeID;
        };

        /**
         * @brief provides access to mapping of data ids to position in node vector
         * @return <SpatID, vector<NodeID> (internal array id)>
         */
        const std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<large_num, large_num>>>&
        getMappingSpatIDtoNodeIDs() {
            return lookupSpatIDtoNodeIDs;
        };

        /**
         * @brief Builds a correct path from the base dir
         * @param path The relative path from the config
         * @return A path based on the base dir
         */
        std::string buildDir(std::string path) {
            return (data_dir / fs::path(path)).string();
        };

        template<class T>
        using Matrix = std::vector<std::vector<T>>;

        static std::vector<Model::quantity<Model::Dimensionless>> calcNusInZones(std::vector<double> densityZones){
            double densityFresh = 1000.0;
            std::vector<Model::quantity<Model::Dimensionless>> nusInZones;
            for (double densityZone : densityZones) {
                // nus of zones is equal to nus of zeta surface below
                nusInZones.push_back(((densityZone - densityFresh) / densityFresh) * Model::si::si_dimensionless);
            }
            return nusInZones;
        };

        static std::vector<Model::quantity<Model::Dimensionless>> calcDelnus(std::vector<double> densityZones){
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
                     large_num numberOfNodesPerLayer,
                     double edgeLengthLeftRight,
                     double edgeLengthFrontBack,
                     int numberOfLayers,
                     double defaultK,
                     double initialHead,
                     double aquiferDepth,
                     double anisotropy,
                     double specificYield,
                     double specificStorage,
                     bool useEfolding,
                     bool confined,
                     large_num maxRefinement,
                     bool isDensityVariable,
                     double effPorosity,
                     double maxTipSlope,
                     double maxToeSlope,
                     double minDepthFactor,
                     double slopeAdjFactor,
                     double vdfLock,
                     std::vector<double> densityZones) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "lon", "lat", "area");
            double lon{0};
            double lat{0};
            double area{0};
            large_num spatID{0};
            large_num nodeID{0};
            large_num refID{0};

            lookupSpatIDtoNodeIDs.reserve(numberOfNodesPerLayer);
            std::vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
            std::vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

            while (in.read_row(spatID, lon, lat, area)) {

                if (edgeLengthLeftRight == edgeLengthFrontBack) {
                    edgeLengthLeftRight = edgeLengthFrontBack = std::sqrt(area);
                }

                nodes->emplace_back(new Model::StandardNode(nodes,
                                                            lat,
                                                            lon,
                                                            area * Model::si::square_meter,
                                                            edgeLengthLeftRight * Model::si::meter,
                                                            edgeLengthFrontBack * Model::si::meter,
                                                            spatID,
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
                                                            maxRefinement,
                                                            isDensityVariable,
                                                            delnus,
                                                            nusInZones,
                                                            effPorosity,
                                                            maxTipSlope,
                                                            maxToeSlope,
                                                            minDepthFactor,
                                                            slopeAdjFactor,
                                                            vdfLock * Model::si::meter));
                for (int layer = 0; layer < numberOfLayers; layer++) {
                    lookupSpatIDtoNodeIDs[spatID][layer][refID] = nodeID + (numberOfNodesPerLayer * layer);
                }
                nodeID++;
            }
            //Return number of total top nodes
            return nodeID - 1;
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
        readLandMaskRefined(NodeVector nodes,
                            std::string path,
                            large_num numberOfNodesPerLayer,
                            double edgeLengthLeftRight,
                            double edgeLengthFrontBack,
                            int numberOfLayers,
                            double defaultK,
                            double initialHead,
                            double aquiferDepth,
                            double anisotropy,
                            double specificYield,
                            double specificStorage,
                            bool useEfolding,
                            bool confined,
                            large_num maxRefinement,
                            bool isDensityVariable,
                            double effPorosity,
                            double maxTipSlope,
                            double maxToeSlope,
                            double minDepthFactor,
                            double slopeAdjFactor,
                            double vdfLock,
                            std::vector<double> densityZones) {
            io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "lon", "lat", "area", "refID"); // todo use new refIDs (1, 11, 111 for deeper levels) only for reading
            double lon{0};
            double lat{0};
            double area{0};
            large_num refID{0};
            large_num nodeID{0};
            large_num spatID{0};

            lookupSpatIDtoNodeIDs.reserve(numberOfNodesPerLayer);
            std::vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
            std::vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

            while (in.read_row(spatID, lon, lat, area, refID)) {

                if (edgeLengthLeftRight == edgeLengthFrontBack) {
                    edgeLengthLeftRight = edgeLengthFrontBack = std::sqrt(area);
                }

                nodes->emplace_back(new Model::StandardNode(nodes,
                                                            lat,
                                                            lon,
                                                            area * Model::si::square_meter,
                                                            edgeLengthLeftRight * Model::si::meter,
                                                            edgeLengthFrontBack * Model::si::meter,
                                                            spatID,
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
                                                            maxRefinement,
                                                            isDensityVariable,
                                                            delnus,
                                                            nusInZones,
                                                            effPorosity,
                                                            maxTipSlope,
                                                            maxToeSlope,
                                                            minDepthFactor,
                                                            slopeAdjFactor,
                                                            vdfLock * Model::si::meter));
                for (int layer = 0; layer < numberOfLayers; layer++) {
                    lookupSpatIDtoNodeIDs[spatID][layer][refID] = nodeID + (numberOfNodesPerLayer * layer);
                }
                nodeID++;
            }
            //Return number of total top nodes
            return nodeID - 1;
        };

        /**
         * @brief Read in a custom definition for the general head boundary
         * @param path Where to read from
         */
        virtual void readHeadBoundary(std::string path) {
            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "elevation", "conduct");
            large_num spatID{0};
            double elevation{0};
            double conduct{0};
            std::unordered_map<large_num, large_num> nodeIDs;
            int layer{0};
            large_num refID{0};

            while (in.read_row(spatID, elevation, conduct)) {
                try {
                    nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }

                for (auto nodeID: nodeIDs) {
                    nodes->at(nodeID.second)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                              elevation * Model::si::meter,
                                                              conduct,
                                                              elevation * Model::si::meter);
                }
            }
        };

        void setVariableDensityConditionsAtBoundary(large_num numZones, double aquiferDepth){
            std::vector<Model::quantity<Model::Meter>> zetas;
            Model::quantity<Model::Meter> zeta;
            for (int nodeID = 0; nodeID < nodes->size(); ++nodeID) {
                if (nodes->at(nodeID)->hasGHB()){
                    /*for (int zetaID = 0; zetaID <= numZones; ++zetaID){
                        double elevation = nodes->at(nodeID)->getElevation().value();

                        if (zetaID == numZones){
                            // set last zeta to bottom of node
                            zeta = (elevation - aquiferDepth) * Model::si::meter;
                        } else {
                            // set all but last zeta to top of node (most saline zone over full node height at boundary)
                            zeta = elevation * Model::si::meter;
                        }
                        nodes->at(nodeID)->setZeta(zetaID, zeta);
                    }*/
                    nodes->at(nodeID)->setZoneOfSinksAndSources(0, numZones-1, numZones);
                } else {
                    nodes->at(nodeID)->setZoneOfSinksAndSources(0, 0, numZones);

                }
            }
        };

        /**
         * @brief Read in a custom definition file for initial heads
         * @param path Where to read the file from
         */
        virtual void readInitialHeads(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                nodes->at(nodeID)->setHead_direct(data);
            });
        };

        /**
         * @brief Read in a custom river definition file
         * Structured as: spatID, Head, Bottom, Conduct
         * @param path Where to read the file from
         */
        virtual void readRiverConductance(std::string path) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "Head", "Bottom", "Conduct");
            large_num spatID{0};
            double head{0};
            double conduct{0};
            double bottom{0};
            std::unordered_map<large_num, large_num> nodeIDs;
            int layer{0};

            while (in.read_row(spatID, head, bottom, conduct)) {
                try {
                    nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                for (auto nodeID : nodeIDs) { // in case the grid is refined: loop over all nodes at refIDs
                    nodes->at(nodeID.second)->addExternalFlow(Model::RIVER,
                                                       head * Model::si::meter,
                                                       conduct,
                                                       bottom * Model::si::meter);
                }
            }
        };

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
        virtual void readElevation(std::string path) {
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
         * @note needs to be in m^3/day
         */
        virtual void readGWRecharge(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                double recharge_m3_per_day = ((data / 1000) * nodes->at(nodeID)->getArea().value()) / 365; // Todo remove / 365 (read data should be mm/day)
                nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                   0 * Model::si::meter,
                                                   recharge_m3_per_day,
                                                   0 * Model::si::meter);
            });
        };

        /**
        * @brief Read diffuse groundwater recharge from a file and map the value using a conversion function
        * @tparam ConversionFunction Allows the dynamic recalculation of recharge based on cell area
        * @param path Where to read the file from
        * @param convertToRate The conversion function
        * @note needs to be in m^3/day
        */
        template<typename ConversionFunction>
        void readGWRechargeMapping(std::string path, ConversionFunction convertToRate) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "arcID", "data");

            int arcID = -1;
            std::vector<large_num> spatIDs;
            std::unordered_map<large_num, large_num> nodeIDs;
            double recharge = 0;
            int layer{0};
            large_num refID{0};

            while (in.read_row(arcID, recharge)) {
                //lookup nodes which get the special_flow
                try {
                    spatIDs = this->lookupArcIDtoSpatIDs.at(arcID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }

                for (large_num spatID : spatIDs) {
                    try {
                        nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                    }
                    catch (const std::out_of_range &ex) { //if Node does not exist ignore entry
                        continue;
                    }

                    for (auto nodeID : nodeIDs) { // in case the grid is refined: loop over all nodes at refIDs

                        //LOG(debug) << nodeID << "," << spatID << "," << arcID;
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
                                                      nodes->at(nodeID.second)->getArea().value());
                        //NANChecker(QMax, "QMax Problem");
                        NANChecker(q_rech, "Recharge-init Problem");

                        nodes->at(nodeID.second)->addExternalFlow(Model::RECHARGE,
                                                                  0 * Model::si::meter,
                                                                  q_rech,
                                                                  0 * Model::si::meter);
                        //rechargeAdded++;
                    }
                }
                spatIDs.clear();
            }
            //LOG(debug) << "missing mapping of arcID to spatID(s): " << missingMapping;
            //LOG(debug) << "recharge added to " << rechargeAdded << " nodes";
        };

        /**
         * @brief Read cell conductance definition
         * @param path Where to read the file from
         */
        virtual void readConduct(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                if (data != 0) { // todo check the impact of the default value on the result! median of glhymps is 0.273
                    //0 is possible data error, known to occur with Gleeson based map
                    nodes->at(nodeID)->setK(data * (Model::si::meter / Model::day)); // Question: check if val > 10 m/day?
                }
            });
        };

        /**
         * @brief Helper function for @see readBlueCells
         * @param path Where to read the file from
         * @return a map of bankfull depth, stream width, and length
         */
        std::unordered_map<large_num, std::array<double, 3>> calculateRiverStage(std::string path) {
            io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "POINTID", "length", "bankfull", "spatID", "width");

            int id{0};
            double length{0};
            double Q_bankfull{0};
            double width{0};
            large_num spatID{0};
            std::unordered_map<large_num, large_num> nodeIDs;
            int layer{0};
            large_num refID{0};

            std::unordered_map<large_num, std::array<double, 3>> out;

            while (in.read_row(id, length, Q_bankfull, spatID, width)) {
                try {
                    nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) { //if Node does not exist ignore entry
                    continue;
                }
                double bankfull = 0.349 * std::pow(Q_bankfull, 0.341);

                for (auto nodeID : nodeIDs) { // in case the grid is refined: loop over all nodes at refIDs
                    out[nodeID.second] = {{bankfull, width, length * 1000}};
                }
            }
            return out;
        };

        /**
         * @brief Reads in river definitions based on a specific elevation data-set
         * @param file to read from
         * @param riverStage A map with additional information @see calculateRiverStage
         */
        void readBlueCells(std::string file,
                           std::unordered_map<large_num, std::array<double, 3>> riverStage) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(file);
            in.read_header(io::ignore_no_column, "spatID", "data");
            large_num spatID{0};
            double riverElevation{0};
            std::unordered_map<large_num, large_num> nodeIDs;
            int layer{0};
            large_num refID{0};

            while (in.read_row(spatID, riverElevation)) {
                try {
                    nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) { //if Node does not exist ignore entry
                    continue;
                }

                for (auto nodeID : nodeIDs) { // in case the grid is refined: loop over all nodes at refIDs

                    double bankfull_depth = riverStage[nodeID.second][0];
                    double width = riverStage[nodeID.second][1];
                    double length = riverStage[nodeID.second][2];

                    if (bankfull_depth <= 1) {
                        bankfull_depth = 1; // river depth is at least 1 meter
                    }
                    double riverBottom = riverElevation - bankfull_depth;
                    double K = nodes->at(nodeID.second)->getK().value();
                    // Conductance estimation following Harbaugh (2005)
                    double conduct = K * length * width / (riverElevation - riverBottom);
                    if (conduct <= 0) { conduct = 1; }


                    nodes->at(nodeID.second)->addExternalFlow(Model::RIVER_MM,
                                                       riverElevation * Model::si::meter,
                                                       conduct,
                                                       riverBottom * Model::si::meter);
                    //LOG(debug) << "nodeID = " << nodeID << ", conduct = " << conduct << ", bottom = " <<
                    //              riverBottom << ", riverElevation = " << riverElevation;
                    //LOG(debug) << "K = " << K << ", length = " << length << ", width = " << width << ", bankfull_depth = "
                    //           << bankfull_depth;
                }
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

            std::vector<std::string> paths = {std::move(pathGlobalLakes), std::move(pathGlobalWetlands),
                                              std::move(pathLocalLakes), std::move(pathLocalWetlands)};

            int itter = 0;
            for (const std::string path : paths) {
                io::CSVReader<2, io::trim_chars<'"', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "data");

                double percentage{0};
                large_num spatID{0};
                std::unordered_map<large_num, large_num> nodeIDs;
                int layer{0};
                large_num refID{0};

                while (in.read_row(spatID, percentage)) {
                    if (percentage == 0) {
                        continue;
                    }
                    percentage = percentage / 100;

                    try {
                        nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    for (auto nodeID: nodeIDs) {

                        if (nodes->at(nodeID.second)->getSpatID() != (int) spatID) {
                            throw "Error in reading spatID";
                        }

                        double elevation = nodes->at(nodeID.second)->getElevation().value();
                        try {
                            elevation = nodes->at(nodeID.second)->getExternalFlowElevation(Model::RIVER_MM).value();
                        } catch (const std::out_of_range &ex) {}

                        double flowElevation = elevation;
                        double bottom = elevation;
                        double K_s = nodes->at(nodeID.second)->getK().value();
                        if (itter == 1) {
                            //global wetlands
                            percentage = percentage * 0.8;
                        }
                        double A_s = nodes->at(nodeID.second)->getArea().value() * percentage;
                        double M = 5;
                        //Simple_ Criv=KLW/M, M is the thickness of the riverbed and K is the hydraulic conductivity of the riverbed
                        double conduct = (K_s * A_s) / M;
                        //if(conduct > 1e6){
                        //    std::cout << "To high K:" << K_s << "area: " << A_s << "at: " << spatID <<"\n";
                        //}

                        //LOG(debug) << "itter = " << itter << ", conduct = " << conduct;
                        if (itter == 0) {
                            //Global LAKE
                            //nodes->at(i)->removeExternalFlow(Model::RIVER_MM);
                            //flowElevation -= 10;
                            bottom -= 100;
                            nodes->at(nodeID.second)->addExternalFlow(Model::GLOBAL_LAKE,
                                                               flowElevation * Model::si::meter,
                                                               conduct,
                                                               bottom * Model::si::meter);
                        } else if (itter == 1) {
                            //Global WETLANDS
                            bottom -= 2;
                            nodes->at(nodeID.second)->addExternalFlow(Model::GLOBAL_WETLAND,
                                                               flowElevation * Model::si::meter,
                                                               conduct,
                                                               bottom * Model::si::meter);
                        } else if (itter == 2) {
                            //Local LAKE
                            //flowElevation -= 10;
                            bottom -= 10;
                            nodes->at(nodeID.second)->addExternalFlow(Model::LAKE,
                                                               flowElevation * Model::si::meter,
                                                               conduct,
                                                               bottom * Model::si::meter);
                        } else {
                            //Local WETLANDS
                            bottom -= 2;
                            nodes->at(nodeID.second)->addExternalFlow(Model::WETLAND,
                                                               flowElevation * Model::si::meter,
                                                               conduct,
                                                               bottom * Model::si::meter);
                        }
                    }
                }
                itter++;
            }
        };

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
        };

        void initializeZetas(large_num numberOfNodes){
            double topOfNode;
            double bottomOfNode;

            // add zeta surfaces to top and bottom of each node
            for (int nodeIter = 0; nodeIter < numberOfNodes; ++nodeIter) {
                topOfNode = nodes->at(nodeIter)->getElevation().value();
                bottomOfNode = nodes->at(nodeIter)->getBottom().value();

                nodes->at(nodeIter)->addZeta(0, topOfNode * Model::si::meter);
                nodes->at(nodeIter)->addZeta(1, bottomOfNode * Model::si::meter);
            }
        }

        /**
         * @brief Compute initial density surfaces using Ghyben-Herzberg relation for interface
         * @param numberOfLayers Number of layers in model
         * @param numberOfNodesPerLayer Number of nodes per layer
         * @param maxDistance Maximum distance from the Ghyben-Herzberg line
         * @param densityZones Density in density zones (from fresh to saline)
         * @note Ghyben-Herzberg: zeta surface = - (density_fresh / (density_saline) - density_fresh)) * GW_head
         */
        void setZetasGhybenHerzberg(int numberOfLayers, large_num numberOfNodesPerLayer, double maxDistance,
                             std::vector<double> densityZones) {
            double initial_head{0.0};
            double zetaGhybenHerzberg{0.0};
            double minDensity = densityZones.front();
            double maxDensity = densityZones.back();
            double meanDensity = (maxDensity + minDensity) * 0.5;
            double maxDifDensity = (maxDensity - minDensity) * 0.5; // calculate maximum difference from mean to min/max
            double zeta{0.0};
            std::vector<double> zetaDeltas;

            large_num numberOfNodes = numberOfLayers * numberOfNodesPerLayer;
            initializeZetas(numberOfNodes); // initialize zeta surface at top and bottom

            for (int localZetaID = 0; localZetaID < densityZones.size(); ++localZetaID) {
                zetaDeltas.push_back( ((meanDensity - densityZones[localZetaID]) / maxDifDensity) * maxDistance );
                //LOG(debug) << "zetaDeltas[localZetaID = " << localZetaID << "]: " << zetaDeltas[localZetaID];
            }
            // loop through nodes
            for (large_num nodeID = 0; nodeID < numberOfNodes; ++nodeID) {
                initial_head = nodes->at(nodeID)->getHead().value();
                // use Ghyben-Herzberg relation to derive zeta surface elevation
                zetaGhybenHerzberg = - (minDensity / (maxDensity - minDensity)) * initial_head;
                //LOG(debug) << "zetaGhybenHerzberg: " << zetaGhybenHerzberg << ", head: " << initial_head;

                for (int localZetaID = 1; localZetaID < densityZones.size(); ++localZetaID){
                    zeta = zetaGhybenHerzberg + zetaDeltas[localZetaID];
                    if (zeta > 0){ zeta = 0; }
                    //LOG(debug) << "zeta[localZetaID = " << localZetaID << "]: " << zeta;
                    nodes->at(nodeID)->addZeta(localZetaID, zeta * Model::si::meter);
                }
            }
        }

        /**
         * @brief Read initial data for density surface height ("zeta") from files
         * @param numberOfLayers Number of layers in model
         * @param path Path to files
         * @param files Vector of file names
         */
        void readInitialZetas(int numberOfLayers, large_num numberOfNodesPerLayer,
                              const std::string& path, std::vector<std::string> files) {

            large_num numberOfNodes = numberOfLayers * numberOfNodesPerLayer;
            initializeZetas(numberOfNodes); // initialize zeta surface at top and bottom

            // read initial data for density surfaces
            loopFilesAndLayers(path, files, numberOfLayers, [this] (std::string path, int numberOfLayers) {
                int localZetaID{0};
                large_num spatID{0};
                std::unordered_map<large_num, large_num> nodeIDs;
                double zeta{0};

                // loop through layers
                for (int layer = 0; layer < numberOfLayers; ++layer) {
                    io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(path);
                    inZetas.read_header(io::ignore_no_column, "spatID", "localZetaID", "zeta");

                    while (inZetas.read_row(spatID, localZetaID, zeta)) {
                        try {
                            nodeIDs = lookupSpatIDtoNodeIDs.at(spatID).at(layer);
                        }
                        catch (const std::out_of_range &ex) { // if node does not exist ignore entry
                            continue;
                        }
                        for (auto nodeID : nodeIDs) { // in case the grid is refined: loop over all nodes at refIDs
                            nodes->at(nodeID.second)->addZeta(localZetaID, zeta * Model::si::meter);
                        }
                    }
                }
            });

        };

        void readEffectivePorosity(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                nodes->at(nodeID)->setEffectivePorosity(data * Model::si::si_dimensionless);
            });
        };

        void readZonesSourcesSinks(std::string path, std::vector<double> densityZones) {
            /**
             * Here we use zoneOfSinks and zoneOfSources (containing values between 0 and number of density zones).
             * Thus, sources and sinks are associated to the respective zone. Rule: zoneOfSinks <= zoneOfSources
             * For simulation of submarine groundwater discharge:
             * - zoneOfSources: an integer between 1 and the number of zones (brackish/saline water)
             * - zoneOfSinks: 0 (fresh water)
             */

            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "zoneOfSinks", "zoneOfSources");
            large_num spatID{0};
            double zoneOfSinks{0};
            double zoneOfSources{0};
            std::unordered_map<int, std::unordered_map<large_num, large_num>> mapAtSpatID;

            int layer{0};
            large_num refID{0};
            while (in.read_row(spatID, zoneOfSinks, zoneOfSources)) {
                try {
                    mapAtSpatID = lookupSpatIDtoNodeIDs.at(spatID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                // loop through layers
                for (const auto& mapAtLayer : mapAtSpatID) {
                    // loop through refIDs
                    for (auto nodeIDs : mapAtLayer.second) { // apply to all layers and refIDs at this spatID
                        nodes->at(nodeIDs.second)->setZoneOfSinksAndSources(zoneOfSinks, zoneOfSources,
                                                                            densityZones.size());
                    }
                }
            }
        };
    };
}
#endif //DATAREADER_HPP
