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
        /** @var lookupSpatIDtoNodeIDs <SpatID, Layer, NodeID>*/
        std::unordered_map<large_num, std::unordered_map<int, large_num>> lookupSpatIDtoNodeID;
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
         * @brief Entry point for reading groundwater head simulation data
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
            double data = 0;
            large_num nodeID;

            int i{0};
            while (in.read_row(spatID, data)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                if (nodes->at(nodeID)->getProperties().get<large_num, Model::SpatID>() != spatID) {
                    throw "Error in reading spatID";
                }
                processData(data, nodeID);
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
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
        void addMappingSpatIDtoNodeID(large_num spatID, int layer, large_num nodeID) {
            lookupSpatIDtoNodeID[spatID][layer] = nodeID;
        };

        /**
         * @brief provides access to mapping of data ids to position in node vector
         * @return <SpatID, vector<NodeID> (internal array id)>
         */
        const std::unordered_map<large_num, std::unordered_map<int, large_num>>&
        getMappingSpatIDtoNodeID() {
            return lookupSpatIDtoNodeID;
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
                     double effPorosity,
                     double maxTipSlope,
                     double maxToeSlope,
                     double minDepthFactor,
                     double slopeAdjFactor,
                     double vdfLock,
                     std::vector<double> densityZones,
                     int sourceZoneGHB,
                     int sourceZoneRecharge) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "lon", "lat", "area");
            double lon{0};
            double lat{0};
            double area{0};
            large_num spatID{0};
            large_num nodeID{0};
            bool isSteadyState{true};
            bool isDensityVariable{false};

            lookupSpatIDtoNodeID.reserve(numberOfNodesPerLayer);
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
                                                            isSteadyState,
                                                            isDensityVariable,
                                                            delnus,
                                                            nusInZones,
                                                            effPorosity,
                                                            maxTipSlope,
                                                            maxToeSlope,
                                                            minDepthFactor,
                                                            slopeAdjFactor,
                                                            vdfLock * Model::si::meter,
                                                            sourceZoneGHB,
                                                            sourceZoneRecharge));
                for (int layer = 0; layer < numberOfLayers; layer++) {
                    lookupSpatIDtoNodeID[spatID][layer] = nodeID + (numberOfNodesPerLayer * layer);
                }
                nodeID++;
            }
            LOG(debug) << "    Added " << nodeID << " nodes";
            //Return number of total top nodes
            return nodeID - 1;
        };

        /**
         * @brief Read in a custom definition for the general head boundary using conductivity; multiplying it with
         * the (1) length of the GHB along the node edges without neighbour and (2) its height; and dividing by distance
         * @param path Where to read from
         * @note transforms hydraulic conductivity [m/day] to conductance [m^2/day]
         */
        virtual void readGHB_elevation_conductivity(std::string path) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "layer", "conductivity", "elevation");
            large_num spatID{0};
            int layer{0};
            double elevation{0};
            double conductivity{0};
            large_num nodeID;
            double conductance{0};
            double ghbLength{0};
            double ghbDistance{0};
            double ghbVerticalSize{0};

            std::vector<Model::NeighbourPosition> neigPos_LRFB;
            std::unordered_map<Model::NeighbourPosition, large_num> horizontal_neighbours;

            int i{0};
            while (in.read_row(spatID, layer, conductivity, elevation)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                // calculate the length and distance of the GHB
                neigPos_LRFB = nodes->at(nodeID)->getNeigPos_LRFB();
                horizontal_neighbours = nodes->at(nodeID)->getListOfHorizontalNeighbours();
                for (const auto &neigPos: neigPos_LRFB) {
                    auto got = horizontal_neighbours.find(neigPos);
                    if (got == horizontal_neighbours.end()) { // no neighbour at position
                        // use square root of area if edge lengths are equal (edge lengths could be 0 if read from file)
                        if (nodes->at(nodeID)->getEdgeLengthLeftRight().value() ==
                            nodes->at(nodeID)->getEdgeLengthFrontBack().value()) {
                            ghbLength += sqrt(nodes->at(nodeID)->getArea().value());
                            ghbDistance += sqrt(nodes->at(nodeID)->getArea().value()) * 0.5;
                        } else {
                            // use individual edge lengths if they are not equal
                            if (neigPos == Model::NeighbourPosition::LEFT or
                                neigPos == Model::NeighbourPosition::RIGHT) {
                                ghbLength += nodes->at(nodeID)->getEdgeLengthLeftRight().value();
                                ghbDistance += nodes->at(nodeID)->getEdgeLengthFrontBack().value() * 0.5;
                            } else {
                                ghbLength += nodes->at(nodeID)->getEdgeLengthFrontBack().value();
                                ghbDistance += nodes->at(nodeID)->getEdgeLengthLeftRight().value() * 0.5;
                            }
                        }
                    }
                }
                ghbVerticalSize = nodes->at(nodeID)->getVerticalSize().value();

                if (conductivity < 1e-11) {
                    conductivity = 1e-11; // set conductivity to a min of 1e-11 m/day
                }
                conductance = conductivity * ghbLength * ghbVerticalSize / ghbDistance; // m/day to m^2/day

                nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                   elevation * Model::si::meter,
                                                   conductance,
                                                   elevation * Model::si::meter);
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
        };

        /**
         * @brief Read in a custom definition for the general head boundary
         * @param path Where to read from
         */
        virtual void readGHB_elevation_conductance(std::string path) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "layer", "conductance", "elevation");
            large_num spatID{0};
            double elevation{0};
            double conductance{0};
            large_num nodeID;
            int layer{0};

            int i{0};
            while (in.read_row(spatID, layer, conductance, elevation)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }

                nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                   elevation * Model::si::meter,
                                                   conductance,
                                                   elevation * Model::si::meter);
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
        };

        /**
         * @brief Read in a custom definition file for initial heads
         * @param path Where to read the file from
         */
        virtual void readInitialHeads(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                if (nodes->at(nodeID)->hasGHB()){
                    double ghbElevation = nodes->at(nodeID)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY);
                    if (data < ghbElevation) {
                        data = ghbElevation;
                    }
                }
                nodes->at(nodeID)->setHead_allLayers(data * Model::si::meter);
                nodes->at(nodeID)->setHead_TZero_allLayers(data * Model::si::meter);
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
            large_num nodeID;
            int layer{0};

            int i{0};
            while (in.read_row(spatID, head, bottom, conduct)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->addExternalFlow(Model::RIVER,
                                                   head * Model::si::meter,
                                                   conduct,
                                                   bottom * Model::si::meter);
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
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
                io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "elevation");
                large_num spatID{0};
                int layer{0};
                double elevation{0};
                large_num nodeID;

                int i{0};
                while (in.read_row(spatID, elevation)) {
                    try {
                        nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }
                    nodes->at(nodeID)->setElevation_allLayers(elevation * Model::si::meter);
                    i++;
                }
                LOG(debug) << "    ... for " << i << " nodes";
            });
        };

        /**
         * @brief Read elevation data from a specified path
         * @note !Uses setElevation() function. Should only be called after all layers are build
         * as it affects layers below
         * @param path Where to read the file from
         */
        virtual void readElevation(std::string path) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "elevation");
            large_num spatID{0};
            int layer{0};
            double elevation{0};
            large_num nodeID{0};

            int i{0};
            while (in.read_row(spatID, elevation)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->setElevation_allLayers(elevation * Model::si::meter);
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
        };

        /**
         * @brief Read equilibrium water-table information used for the dynamic river computation
         * @param path Where to read the file from
         */
        void readEqWTD(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                if (nodes->at(nodeID)->hasGHB()){
                    double ghbElevation = nodes->at(nodeID)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY);
                    if (data < ghbElevation) {
                        data = ghbElevation;
                    }
                }
                nodes->at(nodeID)->setEqHead_allLayers(data * Model::si::meter);
                nodes->at(nodeID)->setHead_TZero_allLayers(data * Model::si::meter);
            });
        };

        /**
         * @brief Read diffuse gw-recharge (input is in mm/day, is transformed to m^3/day)
         * @param path Where to read the file from
         * @note needs to be in m^3/day
         */
        virtual void readGWRecharge(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                double recharge_m3_per_day = ((data / 1000) * nodes->at(nodeID)->getArea().value());
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
            large_num nodeID;
            double recharge{0};
            int layer{0};

            int i{0};
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
                        nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                    }
                    catch (const std::out_of_range &ex) { //if Node does not exist ignore entry
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
                                                  nodes->at(nodeID)->getArea().value());
                    //NANChecker(QMax, "QMax Problem");
                    NANChecker(q_rech, "Recharge-init Problem");

                    nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                       0 * Model::si::meter,
                                                       q_rech,
                                                       0 * Model::si::meter);
                    i++;
                }
                spatIDs.clear();
            }
            LOG(debug) << "    ... for " << i << " nodes";
        };

        /**
         * @brief Read cell conductivity definition
         * @param path Where to read the file from
         */
        virtual void readConductivity(std::string path) {
            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "layer", "conductivity");
            large_num spatID{0};
            int layer{0};
            double conductivity{0};
            large_num nodeID{0};
            //double threshold{1e-11}; // todo debug Equation: make preconditioner run with threshold (or delete such nodes)
            int i{0};
            int j{0};
            while (in.read_row(spatID, layer, conductivity)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                if (conductivity < 1e-11){
                    conductivity = 1e-11;
                }
                nodes->at(nodeID)->setK_allLayers(conductivity * (Model::si::meter / Model::day));

                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
        };

        /**
         * @brief Helper function for @see readBlueCells
         * @param path Where to read the file from
         * @return a map of bankfull depth, stream width, and length
         */
        std::array<double, 3> calculateRiverStage(std::string path) {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "length", "width", "bankfull");

            large_num spatID{0};
            double riverLength{0};
            double riverWidth{0};
            double Q_bankfull{0};
            int layer{0};
            large_num nodeID;
            std::array<double, 3> out;

            while (in.read_row(spatID, riverLength, riverWidth, Q_bankfull)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer); // layer is 0 since rivers are at surface
                }
                catch (const std::out_of_range &ex) { //if Node does not exist ignore entry
                    continue;
                }
                double bankfull_depth = 0.349 * std::pow(Q_bankfull, 0.341);

                out = {{bankfull_depth, riverWidth, riverLength * 1000}};
            }
            return out;
        };

        /**
         * @brief Reads in river definitions based on a specific elevation data-set
         * @param file to read from
         * @param riverStage A map with additional information @see calculateRiverStage
         */
        void readBlueCells(std::string file, std::array<double, 3> riverStage) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(file);
            in.read_header(io::ignore_no_column, "spatID", "data");
            large_num spatID{0};
            double riverElevation{0};
            double riverWidth{0};
            double riverLength{0};
            double bankfull_depth{0};
            double riverConductance{0};
            double riverDepth{0};
            double riverBottom{0};
            large_num nodeID;
            int layer{0};

            int i{0};
            while (in.read_row(spatID, riverElevation)) {
                try {
                    nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                }
                catch (const std::out_of_range &ex) { //if Node does not exist ignore entry
                    continue;
                }

                bankfull_depth = riverStage[0];
                riverWidth = riverStage[1];
                riverLength = riverStage[2];

                if (bankfull_depth <= 1) {
                    bankfull_depth = 1; // river depth is at least 1 meter
                }

                if (nodes->at(nodeID)->hasGHB()){
                    double ghbElevation = nodes->at(nodeID)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY);
                    if (riverElevation < ghbElevation){
                        riverElevation = ghbElevation; // river elevation cannot be below GHB elevation
                    }
                }
                riverBottom = riverElevation - bankfull_depth;
                double K = nodes->at(nodeID)->getK().value();
                // Conductance estimation following Harbaugh (2005)
                riverConductance = K * riverLength * riverWidth / bankfull_depth;
                if (riverConductance <= 1) { riverConductance = 1; } // river conductance is at least 1

                nodes->at(nodeID)->addExternalFlow(Model::RIVER_MM,
                                                          riverElevation * Model::si::meter,
                                                          riverConductance,
                                                          riverBottom * Model::si::meter);
                //LOG(debug) << "nodeID = " << nodeID << ", conduct = " << conduct << ", bottom = " <<
                //              riverBottom << ", riverElevation = " << riverElevation;
                //LOG(debug) << "K = " << K << ", length = " << length << ", width = " << width << ", bankfull_depth = "
                //           << bankfull_depth;
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
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
                large_num nodeID;
                int layer{0};

                while (in.read_row(spatID, percentage)) {
                    if (percentage == 0) {
                        continue;
                    }
                    percentage = percentage / 100;

                    try {
                        nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }

                    if (nodes->at(nodeID)->getSpatID() != (int) spatID) {
                        throw "Error in reading spatID";
                    }

                    double elevation = nodes->at(nodeID)->getExternalFlowElevation(Model::RIVER_MM);
                    if (std::isnan(elevation)){
                        elevation = nodes->at(nodeID)->getElevation().value();
                    }

                    double flowElevation = elevation;
                    double bottom = elevation;
                    double K_s = nodes->at(nodeID)->getK().value();
                    if (itter == 1) {
                        //global wetlands
                        percentage = percentage * 0.8;
                    }
                    double A_s = nodes->at(nodeID)->getArea().value() * percentage;
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
                        nodes->at(nodeID)->addExternalFlow(Model::GLOBAL_LAKE,
                                                           flowElevation * Model::si::meter,
                                                           conduct,
                                                           bottom * Model::si::meter);
                    } else if (itter == 1) {
                        //Global WETLANDS
                        bottom -= 2;
                        nodes->at(nodeID)->addExternalFlow(Model::GLOBAL_WETLAND,
                                                           flowElevation * Model::si::meter,
                                                           conduct,
                                                           bottom * Model::si::meter);
                    } else if (itter == 2) {
                        //Local LAKE
                        //flowElevation -= 10;
                        bottom -= 10;
                        nodes->at(nodeID)->addExternalFlow(Model::LAKE,
                                                           flowElevation * Model::si::meter,
                                                           conduct,
                                                           bottom * Model::si::meter);
                    } else {
                        //Local WETLANDS
                        bottom -= 2;
                        nodes->at(nodeID)->addExternalFlow(Model::WETLAND,
                                                           flowElevation * Model::si::meter,
                                                           conduct,
                                                           bottom * Model::si::meter);
                    }
                }
                itter++;
            }
        };


        void addRiverWhereNoSurfaceWaterBody(double swbElevationFactor, double riverConductivity){
            double riverElevation{0};
            double riverWidth{0};
            double riverLength{0};
            double riverConductance{0};
            double riverDepth{0};
            double riverBottom{0};
            int layer{0};
            large_num spatID{0};

            int i = 0;
            for (const auto &node : *nodes) {
                // if river/lake/wetland exists in node: move on to next node
                if (node->hasTypeOfExternalFlow(Model::RIVER_MM) or
                    node->hasTypeOfExternalFlow(Model::LAKE) or
                    node->hasTypeOfExternalFlow(Model::WETLAND) or
                    node->hasTypeOfExternalFlow(Model::GLOBAL_LAKE) or
                    node->hasTypeOfExternalFlow(Model::GLOBAL_WETLAND)){
                    continue;
                } else { // no river in node
                    spatID = node->getSpatID();
                    layer = 0;
                    // if we are on the top layer
                    if (node->getID() == lookupSpatIDtoNodeID.at(spatID).at(layer)){
                        riverWidth = 25.7; // based on mean of global data
                        riverLength = 10.9 * 1000; // based on mean of global data
                        riverElevation = swbElevationFactor * node->getElevation().value();
                        riverDepth = 2.1; // based on mean of Q_bankfull (=204); depth = 0.349 * std::pow(204, 0.341)
                        riverBottom = riverElevation - riverDepth;
                        // Conductance estimation following Harbaugh (2005)
                        riverConductance = riverConductivity * riverLength * riverWidth / (riverElevation - riverBottom);
                        node->addExternalFlow(Model::RIVER_MM,
                                                           riverElevation * Model::si::meter,
                                                           riverConductance,
                                                           riverBottom * Model::si::meter);
                        i++;
                    }
                }
            }
            LOG(debug) << "    Added river (since no surface water body was defined by data) to " << i << " nodes";
        }

        /**
         * @brief Adds an evapotranspiration module to cells
         * @deprecated
         */
        void addEvapo() {
            int i{0};
            for (const auto &node : *nodes.get()) {
                if (node->getProperties().get<int, Model::Layer>() == 0 and
                    node->getProperties().get<Model::quantity<Model::Meter>,
                            Model::Elevation>().value() <= 600) {
                    node->addExternalFlow(Model::EVAPOTRANSPIRATION,
                                          node->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>(),
                                          0.001 *
                                          node->getProperties().get<Model::quantity<Model::SquareMeter>, Model::Area>().value(),
                                          0.5 * Model::si::meter);
                }
                i++;
            }
            LOG(debug) << "    ... for " << i << " nodes";
        };


        /**
         * @brief Compute initial density surfaces using Ghyben-Herzberg relation for interface
         * @param numberOfLayers Number of layers in model
         * @param numberOfNodesPerLayer Number of nodes per layer
         * @param maxDistance Maximum distance from the Ghyben-Herzberg line
         * @param densityZones Density in density zones (from fresh to saline)
         * @note Ghyben-Herzberg: zeta surface = - (density_fresh / (density_saline) - density_fresh)) * GW_head
         */
        /*void setZetasGhybenHerzberg(int numberOfLayers, large_num numberOfNodesPerLayer,
                                    std::vector<double> densityZones) {
            double ghybenHerzberg{0.0};
            double minDensity = densityZones.front();
            double maxDensity = densityZones.back();
            double meanDensity = (maxDensity + minDensity) * 0.5;
            double maxDifDensity = (maxDensity - minDensity) * 0.5; // calculate maximum difference from mean to min/max
            double zeta{0.0};
            std::vector<double> zetaDeltas;
            double ghbElevation;
            large_num numberOfNodes = numberOfLayers * numberOfNodesPerLayer;

            double maxDistance = 10;
            for (double densityZone : densityZones) {
                zetaDeltas.push_back( ((meanDensity - densityZone) / maxDifDensity) * maxDistance );
                //LOG(debug) << "zetaDeltas[zetaID = " << zetaID << "]: " << zetaDeltas[zetaID];
            }

            // loop through nodes
            for (const auto &node : *nodes) {
                // 1. initialize zeta surface at top and bottom
                node->initializeZetas();

                // 2. use Ghyben-Herzberg relation to derive 50% interface elevation
                ghybenHerzberg = - (minDensity / (maxDensity - minDensity)) * node->getHead().value();
                //LOG(debug) << "zetaGhybenHerzberg: " << zetaGhybenHerzberg << ", head: " << initial_head;

                // 3. set zeta surfaces (potentially between top and bottom)
                for (int zetaID = 1; zetaID < densityZones.size(); ++zetaID){
                    // if node has a GHB
                    if (node->hasGHB()){
                        zeta = ghybenHerzberg + zetaDeltas[zetaID];
                        // if zeta above GHB: set to GHB elevation
                        ghbElevation = node->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY);
                        if (zeta > ghbElevation) { zeta = ghbElevation; }

                    // if node has no GHB: set to node bottom
                    } else {
                        zeta = node->getBottom().value();
                    }

                    node->addZeta(zetaID, zeta * Model::si::meter);
                }
            }
        }*/


        void setDefaultZetas(int numberOfZones) {
            for (const auto &node : *nodes) {
                // initialize zeta surface at top and bottom
                node->initializeZetas();
                for (int zetaID = 1; zetaID < numberOfZones; ++zetaID){
                    node->addZeta(zetaID, node->getBottom()); // add additional zetas at bottom
                }
            }
        }

        /**
         * @brief Read initial data for density surface height ("zeta") from files
         * @param numberOfLayers Number of layers in model
         * @param path Path to files
         * @param files Vector of file names
         */
        void readInitialZetas(large_num numberOfLayers, large_num numberOfNodesPerLayer,
                              const std::string& path, std::vector<std::string> files) {

            for (const auto &node : *nodes) {
                // initialize zeta surface at top and bottom
                node->initializeZetas();
            }

            // read initial data for density surfaces
            loopFilesAndLayers(path, files, numberOfLayers, [this] (std::string path, int numberOfLayers) {
                int zetaID{0};
                large_num spatID{0};
                large_num nodeID;
                double zeta{0};

                // loop through layers
                for (int layer = 0; layer < numberOfLayers; ++layer) {
                    io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> inZetas(path);
                    inZetas.read_header(io::ignore_no_column, "spatID", "zetaID", "zeta");

                    while (inZetas.read_row(spatID, zetaID, zeta)) {
                        try {
                            nodeID = lookupSpatIDtoNodeID.at(spatID).at(layer);
                        }
                        catch (const std::out_of_range &ex) { // if node does not exist ignore entry
                            continue;
                        }
                        nodes->at(nodeID)->addZeta(zetaID, zeta * Model::si::meter);
                    }
                }
            });

        };

        void readEffectivePorosity(std::string path) {
            readTwoColumns(path, [this](double data, int nodeID) {
                // snap very low porosity values to 0
                if (data < 1e-3){ data = 0; }
                nodes->at(nodeID)->setEffectivePorosity(data * Model::si::si_dimensionless);
            });
        };
    };
}
#endif //DATAREADER_HPP
