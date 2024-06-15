#ifndef TESTING_SIMPLEDATAREADER_HPP
#define TESTING_SIMPLEDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"

namespace GlobalFlow {
namespace DataProcessing {


class SimpleRefinedDataReader : public DataReader {
    public:
    SimpleRefinedDataReader() = default;

        void readData(Simulation::Options op) override {
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
                            op.getMaxRefinement(),
                            op.getEffectivePorosity(),
                            op.getMaxTipSlope(),
                            op.getMaxToeSlope(),
                            op.getMinDepthFactor(),
                            op.getSlopeAdjFactor(),
                            op.getVDFLock(),
                            op.getDensityZones(),
                            op.getSourceZoneGHB(),
                            op.getSourceZoneRecharge());

            if (op.getNumberOfLayers() > 1) {
                LOG(userinfo) << "Building the model layer(s) below";
                DataProcessing::buildBottomLayers(nodes,
                                                  op.getNumberOfLayers(),
                                                  op.getConfinements(),
                                                  op.getAquiferDepth(),
                                                  op.getInitialK(),
                                                  op.getAnisotropy());
            }

            LOG(userinfo) << "Building grid by spatial ID (refined)";
            DataProcessing::buildBySpatID(nodes,
                                          this->getMappingSpatIDtoNodeIDs(),
                                          op.getResolution(),
                                          op.getXRange(),
                                          op.getYRange(),
                                          op.isGlobal(),
                                          op.getNumberOfLayers(),
                                          op.getNumberOfNodesPerLayer(),
                                          op.getGHBConduct(),
                                          op.getBoundaryCondition());

            if (op.getNumberOfLayers() > 1) {
                LOG(userinfo) << "Copying neighbours to bottom layer(s)";
                DataProcessing::copyNeighboursToBottomLayers(nodes, op.getNumberOfLayers());
            }

            LOG(userinfo) << "Reading lithology";
            readConductivity(buildDir(op.getLithology()));

            LOG(userinfo) << "Reading elevation";
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Initializing head";
            readInitialHeads((buildDir(op.getInitialHeadsDir())));

            if(op.isKRiverFromFile()) {
                LOG(userinfo) << "Defining rivers";
                readRiverConductance(buildDir(op.getKRiver()));
            }

            if(op.isKGHBFromFile()) {
                LOG(userinfo) << "Reading the boundary condition";
                readGHB_elevation_conductance(buildDir(op.getKGHBDir()));
            }
        }

    private:
        void readElevation(std::string path) override {
            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "data", "refID");
            large_num spatID{0};
            double data{0};
            large_num refID{0};
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
                nodes->at(nodeID)->setElevation_allLayers(data * Model::si::meter);
            }
        };

        void readRiverConductance(std::string path) override {
            io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "Head", "Bottom", "Conduct", "refID");
            large_num spatID{0};
            double head{0};
            double conduct{0};
            double bottom{0};
            large_num nodeID;
            int layer{0};
            large_num refID{0};

            while (in.read_row(spatID, head, bottom, conduct, refID)) {
                try {
                    nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(layer).at(refID);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(nodeID)->addExternalFlow(Model::RIVER,
                                                   head * Model::si::meter,
                                                   conduct,
                                                   bottom * Model::si::meter);
            }
        }

        void readGWRecharge(std::string path) override {
            io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "spatID", "layer", "recharge", "refID");
            large_num spatID{0};
            int layer{0};
            double recharge{0};
            large_num refID{0};
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
                                                   recharge * nodes->at(nodeID)->getArea().value(),
                                                0 * Model::si::meter);
            }
        }

    void readInitialHeads(std::string path) override {
        io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
        in.read_header(io::ignore_no_column, "spatID", "data", "refID");
        large_num spatID{0};
        double data{0};
        int layer{0};
        large_num refID{0};
        large_num nodeID;

        while (in.read_row(spatID, data, refID)) {
            try {
                nodeID = lookupSpatIDtoNodeIDs.at(spatID).at(layer).at(refID);
            }
            catch (const std::out_of_range &ex) {
                //if Node does not exist ignore entry
                continue;
            }
            nodes->at(nodeID)->setHead_allLayers(data * Model::si::meter);
        }
    }
};
}
}
#endif //TESTING_SIMPLEDATAREADER_HPP
