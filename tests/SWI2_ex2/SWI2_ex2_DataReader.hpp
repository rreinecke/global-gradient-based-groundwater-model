#ifndef TESTING_SWI2_EX2_DATAREADER_HPP
#define TESTING_SWI2_EX2_DATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"
#include "../../src/Model/Node.hpp"
#include "../../src/Model/Units.hpp"

namespace GlobalFlow {
    namespace DataProcessing {

        class SWI2_ex2_DataReader : public DataReader {
        public:
            SWI2_ex2_DataReader() = default;

            void readData(Simulation::Options op) override {
                LOG(userinfo) << "Reading land mask (with default values from config)";
                readLandMask(nodes, buildDir(op.getNodesDir()), op.getNumberOfNodesPerLayer(),
                             op.getEdgeLengthLeftRight(), op.getEdgeLengthFrontBack(),
                             op.getNumberOfLayers(), op.getInitialK()[0], op.getInitialHead(),op.getAquiferDepth()[0],
                             op.getAnisotropy()[0], op.getSpecificYield(), op.getSpecificStorage(), op.useEfolding(),
                             op.isConfined(0),
                             op.getEffectivePorosity(), op.getMaxTipSlope(), op.getMaxToeSlope(),
                             op.getMinDepthFactor(), op.getSlopeAdjFactor(), op.getVDFLock(), op.getDensityZones(),
                             op.getSourceZoneGHB(), op.getSourceZoneRecharge());

                LOG(userinfo) << "Building grid by spatial ID";
                DataProcessing::buildBySpatID(nodes,
                                              this->getMappingSpatIDtoNodeID(),
                                              op.getResolution(),
                                              op.getXRange(),
                                              op.getYRange(),
                                              op.isGlobal(),
                                              op.getNumberOfLayers(),
                                              op.getNumberOfNodesPerLayer(),
                                              op.getGHBConduct(),
                                              op.getBoundaryCondition());

                LOG(userinfo) << "Reading elevation";
                readElevation(buildDir(op.getElevation()));

                if (op.isInitialHeadFromFile()){
                    LOG(userinfo) << "Initializing head";
                    readInitialHeads((buildDir(op.getInitialHeadsDir())));
                }

                if (op.isKFromFile()) {
                    LOG(userinfo) << "Reading hydraulic conductivity";
                    readConductivity(buildDir(op.getLithology()));
                }

                LOG(userinfo) << "Reading the boundary condition";
                readGHB_elevation_conductance(buildDir(op.getKGHBDir()));

                if (op.isInitialZetasAsArray()) {
                    LOG(userinfo) << "Reading zetas";
                    readInitialZetas(op.getNumberOfLayers(), op.getDensityZones().size(),
                                     buildDir(op.getInitialZetas()), op.getInitialZetas_a());
                }

                if (op.isEffectivePorosityFromFile()){
                    LOG(userinfo) << "Reading effective porosity";
                    readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
                }
            }

        private:
        };
    }
}
#endif //TESTING_SWI2_EX2_2DATAREADER_HPP
