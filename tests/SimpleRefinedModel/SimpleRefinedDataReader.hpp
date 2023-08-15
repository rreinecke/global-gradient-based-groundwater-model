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
            readLandMaskRefined(nodes,
                            buildDir(op.getNodesDir()),
                            op.getNumberOfNodesPerLayer(),
                            op.getEdgeLengthLeftRight(),
                            op.getEdgeLengthFrontBack(),
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

            LOG(userinfo) << "Building grid by spatial ID (refined)";
            DataProcessing::buildBySpatID(nodes,
                                          this->getMappingSpatIDtoNodeIDs(),
                                          1, // resolution = 0.0833 <- input for global models
                                          10, // lonRange = 360
                                          10, // latRange = 180
                                          false, // isGlobal = true
                                          op.getNumberOfNodesPerLayer(),
                                          op.getGHBConduct(),
                                          op.getBoundaryCondition());

            if (op.getNumberOfLayers() > 1) {
                LOG(userinfo) << "Building the model layer(s) below";
                DataProcessing::buildBottomLayers(nodes,
                                                  op.getNumberOfLayers(),
                                                  op.getConfinements(),
                                                  op.getAquiferDepth(),
                                                  op.getInitialK(),
                                                  op.getAnisotropy());

                LOG(userinfo) << "Copying neighbours to bottom layer(s)";
                DataProcessing::copyNeighboursToBottomLayers(nodes, op.getNumberOfLayers());
            }

            LOG(userinfo) << "Reading lithology";
            readConduct(buildDir(op.getLithology()));

            LOG(userinfo) << "Reading elevation";
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Initializing head";
            readInitialHeads((buildDir(op.getInitialHeadsDir())));

            LOG(userinfo) << "Defining rivers";
            readRiverConductance(buildDir(op.getKRiver()));

            LOG(userinfo) << "Adding GHB conductance";
            //readHeadBoundary(buildDir(op.getKGHBDir()));

        }

    private:
        template<class T>
        using Matrix = std::vector<std::vector<T>>;

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
            io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "X", "Y", "area", "refID"); // todo use new refIDs (1, 11, 111 for deeper levels) only for reading
            double X{0};
            double Y{0};
            double area{0};
            int refID{0};
            large_num nodeID{0};
            large_num spatID{0};

            lookupSpatIDtoNodeIDs.reserve(numberOfNodesPerLayer);
            std::vector<Model::quantity<Model::Dimensionless>> delnus = calcDelnus(densityZones);
            std::vector<Model::quantity<Model::Dimensionless>> nusInZones = calcNusInZones(densityZones);

            while (in.read_row(spatID, X, Y, area, refID)) {
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
            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "data", "refID");
            large_num spatID{0};
            double data{0};
            int refID{0};
            int layer{0};
            large_num nodeID;

            while (in.read_row(spatID, data, refID)) {
                try {
                    nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(layer).at(refID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->setK(data * (Model::si::meter / Model::day));
            }
        };

        void readElevation(std::string path) {
            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "data", "refID");
            large_num spatID{0};
            double data{0};
            int refID{0};
            int layer{0};
            large_num nodeID;

            while (in.read_row(spatID, data, refID)) {
                try {
                    nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(layer).at(refID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->setElevation(data * Model::si::meter);
            }
        };

        void readRiverConductance(std::string path) {
            io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "Head", "Bottom", "Conduct", "refID");
            large_num spatID{0};
            double head{0};
            int refID{0};
            double conduct{0};
            double bottom{0};
            large_num nodeID;

            while (in.read_row(spatID, head, bottom, conduct, refID)) {
                try {
                    nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(0).at(refID); //layer = 0
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
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "layer", "recharge", "refID");
            large_num spatID{0};
            int layer{0};
            double recharge{0};
            int refID{0};
            large_num nodeID;

            while (in.read_row(spatID, layer,recharge, refID)) {
                try {
                    nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(layer).at(refID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                0 * Model::si::meter,
                                                   recharge * nodes->at(nodeID)->getProperties().get < Model::quantity <
                                                        Model::SquareMeter > ,
                                                        Model::Area > ().value(),
                                                0 * Model::si::meter);
            }
        }

    void readInitialHeads(std::string path) {
        io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
        in.read_header(io::ignore_no_column, "spatID", "data", "refID");
        large_num spatID{0};
        double data{0};
        int refID{0};
        large_num nodeID;

        while (in.read_row(spatID, data, refID)) {
            try {
                nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(0).at(refID); //layer = 0
            }
            catch (const std::out_of_range &ex) {
                //if Node does not exist ignore entry
                continue;
            }
            nodes->at(nodeID)->setHead_direct(data);
        }
    }

    void readHeadBoundary(std::string path) {
        io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
        in.read_header(io::ignore_no_column, "spatID", "head", "conduct", "refID");
        large_num spatID{0};
        double head{0};
        double conduct{0};
        large_num nodeID;
        int refID{0};

        while (in.read_row(spatID, head, conduct, refID)) {
            try {
                nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(0).at(refID); // layer = 0
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
