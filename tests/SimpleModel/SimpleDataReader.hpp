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
                LOG(userinfo) << "Building the bottom layers";
                DataProcessing::buildBottomLayers(nodes,
                                                  op.getNumberOfLayers(),
                                                  op.getConfinements(),
                                                  op.getAquiferDepth(),
                                                  op.getInitialK(),
                                                  op.getAnisotropy());
            }

            LOG(userinfo) << "Building grid by spatial ID";
            DataProcessing::buildBySpatID(nodes,
                                          this->getMappingSpatIDtoNodeID(),
                                          op.getResolution(),
                                          op.getXRange(),
                                          op.getYRange(),
                                          op.isGlobal(),
                                          op.getNumberOfLayers(),
                                          op.getNumberOfNodesPerLayer(),
                                          op.getGHBConduct(), // default conductance of general head boundary
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

            LOG(userinfo) << "Defining rivers";
            readRiverConductance(buildDir(op.getKRiver()));

            if(op.isKGHBFromFile()) {
                LOG(userinfo) << "Reading the boundary condition";
                readGHB_elevation_conductance(buildDir(op.getKGHBDir()));
            }
        }

    private:

};
}
}
#endif //TESTING_SIMPLEDATAREADER_HPP
