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
         */
        GlobalDataReader() = default;

        /**
         * @brief overriding read data for the large model
         * @param op
         * @note Number of nodes per layer for global: 2161074, for North America: 396787, for New Zealand: 4603
         */
        void readData(Simulation::Options op) override {
            LOG(userinfo) << "Reading land mask (with default values from config)";
            readLandMask(nodes, buildDir(op.getNodesDir()), op.getNumberOfNodesPerLayer(),
                         op.getEdgeLengthLeftRight(), op.getEdgeLengthFrontBack(),
                         op.getNumberOfLayers(), op.getInitialK()[0], op.getInitialHead(),op.getAquiferDepth()[0],
                         op.getAnisotropy()[0], op.getSpecificYield(), op.getSpecificStorage(), op.useEfolding(),
                         op.isConfined(0), op.isDensityVariable(),
                         op.getEffectivePorosity(), op.getMaxTipSlope(), op.getMaxToeSlope(),
                         op.getMinDepthFactor(), op.getSlopeAdjFactor(), op.getVDFLock(), op.getDensityZones());

            if (op.getNumberOfLayers() > 1) {
                LOG(userinfo) << "Building the model layer(s) below";
                DataProcessing::buildBottomLayers(nodes,
                                                  op.getNumberOfLayers(),
                                                  op.getConfinements(),
                                                  op.getAquiferDepth(),
                                                  op.getInitialK(),
                                                  op.getAnisotropy());
            }

            LOG(userinfo) << "Setting neighbours by SpatID";
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

            if (op.getNumberOfLayers() > 1) {
                LOG(userinfo) << "Copying neighbours to bottom layer(s)";
                DataProcessing::copyNeighboursToBottomLayers(nodes, op.getNumberOfLayers());

                if (op.useEfolding()) {
                    LOG(userinfo) << "Reading e-folding";
                    readEfold(buildDir(op.getEfolding()), op.getEfolding_a());
                }
            }

            LOG(userinfo) << "Reading mapping of SpatID to ArcID";
            readSpatIDtoArcID(buildDir(op.getMapping()));

            if(op.isKGHBFromFile()) {
                LOG(userinfo) << "Reading the boundary condition";
                readHeadBoundary(buildDir(op.getKGHBDir()));
            }

            if (op.isKFromFile()) {
                LOG(userinfo) << "Reading hydraulic conductivity";
                readConduct(buildDir(op.getLithology()));
            }

            LOG(userinfo) << "Reading elevation";
            readElevation(buildDir(op.getElevation()));

            // read either initial head (default) or equilibrium water table depth from file, if available
            if (op.isInitialHeadFromFile()){
                LOG(userinfo) << "Reading initial head";
                readInitialHeads((buildDir(op.getInitialHeadsDir())));
            } else if (op.isEqWTDFromFile()){
                LOG(userinfo) << "Reading equal water table depth";
                readEqWTD(buildDir(op.getEqWTD())); // requires elevation to be set
            }

            LOG(userinfo) << "Reading groundwater recharge";
            //readGWRecharge(buildDir(op.getRecharge()));
            readGWRechargeMapping(buildDir(op.getRecharge()),
                                  [](const double &recharge, const double &area) {
                                      return (((recharge / 1000) * area) / 365);});

            LOG(userinfo) << "Reading rivers";
            if (op.isKRiverFromFile()) {
                readRiverConductance(buildDir(op.getKRiver()));
            } else {
                readBlueCells(buildDir(op.getRiverElevation()),
                              calculateRiverStage(buildDir(op.getRiverExtent())));
            }


            LOG(userinfo) << "Reading lakes and wetlands"; // should be placed after readBlueCells
            readLakesAndWetlands(buildDir(op.getGlobalLakes()),
                                 buildDir(op.getGlobalWetlands()),
                                 buildDir(op.getLocalLakes()),
                                 buildDir(op.getLocalWetlands()));

            if (op.isDensityVariable()) {
                LOG(userinfo) << "Reading initial zeta heights";
                readInitialZetas(op.getNumberOfNodesPerLayer(), op.getNumberOfLayers(),
                                 buildDir(op.getInitialZetas()), op.getInitialZetas_a()); // requires elevation to be set

                if (op.isEffectivePorosityFromFile()) {
                    LOG(userinfo) << "Reading effective porosity";
                    readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
                }

                if (op.isZonesSourcesSinksFromFile()) {
                    LOG(userinfo) << "Reading zones of sources and sinks";
                    readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinksDir()),
                                          op.getDensityZones());
                }

                LOG(userinfo) << "Setting the variable density conditions at general head boundaries";
                // Needs to be called after GHB was set (after buildByGrid/buildBySpatID and readHeadBoundary)
                setVariableDensityConditionsAtBoundary(op.getDensityZones().size(),
                                                       op.getAquiferDepth()[0]);
            }

        }
    };
}
}
