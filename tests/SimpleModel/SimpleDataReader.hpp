#ifndef TESTING_SIMPLEDATAREADER_HPP
#define TESTING_SIMPLEDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"

namespace GlobalFlow {
namespace DataProcessing {


class SimpleDataReader : public DataReader {
    public:
        SimpleDataReader() = default;

        void readData(Simulation::Options op) override {
            LOG(userinfo) << "Building the initial model layer";
            readLandMask(nodes,
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

            LOG(userinfo) << "Building the bottom layers";
            DataProcessing::buildBottomLayers(nodes,
                                              op.getNumberOfLayers(),
                                              op.getConfinements(),
                                              op.getAquiferDepth(),
                                              op.getInitialK(),
                                              op.getAnisotropy());

            LOG(userinfo) << "Building grid by spatial ID";
            DataProcessing::buildBySpatID(nodes,
                                          this->getMappingSpatIDtoNodeIDs(),
                                          op.getResolution(),
                                          op.getLonRange(),
                                          op.getLatRange(),
                                          op.isGlobal(),
                                          op.getNumberOfLayers(),
                                          op.getNumberOfNodesPerLayer(),
                                          op.getGHBConduct(),
                                          op.getBoundaryCondition());

            LOG(userinfo) << "Reading hydraulic parameters";
            readConduct(buildDir(op.getLithology()));
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Initializing head";
            readInitialHeads((buildDir(op.getInitialHeadsDir())));

            LOG(userinfo) << "Defining rivers";
            readRiverConductance(buildDir(op.getKRiver()));

            //LOG(userinfo) << "Building grid by rows and columns (boundaries need to be specified in with a file)";
            //DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodesPerLayer(), op.getNumberOfLayers());
        }

    private:

};
}
}
#endif //TESTING_SIMPLEDATAREADER_HPP
