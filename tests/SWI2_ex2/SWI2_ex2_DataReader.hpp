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

                LOG(userinfo) << "Reading elevation";
                readElevation(buildDir(op.getElevation()));

                if (op.isKFromFile()) {
                    LOG(userinfo) << "Reading hydraulic conductivity";
                    readConduct(buildDir(op.getLithology()));
                }

                //LOG(userinfo) << "Reading the groundwater recharge";
                //readGWRecharge(buildDir(op.getRecharge()));

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

        };
    }
}
#endif //TESTING_SWI2_EX2_2DATAREADER_HPP
