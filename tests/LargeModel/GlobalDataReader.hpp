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
            std::vector<int> stressPeriodSteps = op.getStressPeriodSteps();
            int n_totalSteps = std::accumulate(stressPeriodSteps.begin(), stressPeriodSteps.end(), 0);
            LOG(userinfo) << "Total number of time steps simulated: " << n_totalSteps;
            LOG(userinfo) << "Reading land mask (with default values from config)";

            /*
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             * %%% build grid (read grid, set neighbors, add lower layers) %%%
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             */
            readLandMask(nodes, buildDir(op.getNodesDir()), op.getNumberOfNodesPerLayer(),
                         op.getEdgeLengthLeftRight(), op.getEdgeLengthFrontBack(),
                         op.getNumberOfLayers(), op.getInitialK()[0], op.getInitialHead(), op.getAquiferDepth()[0],
                         op.getAnisotropy()[0], op.getSpecificYield(), op.getSpecificStorage(), op.useEfolding(),
                         op.isConfined(0),
                         op.getEffectivePorosity(), op.getMaxTipSlope(), op.getMaxToeSlope(),
                         op.getMinDepthFactor(), op.getSlopeAdjFactor(), op.getVDFLock(), op.getDensityZones(),
                         op.getSourceZoneGHB(), op.getSourceZoneRecharge());

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
                                          this->getMappingSpatIDtoNodeID(),
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
                DataProcessing::copyNeighboursToBottomLayers(nodes,
                                                             op.getNumberOfLayers()); // todo is it possible to include this in buildBySpatID

                if (op.useEfolding()) {
                    LOG(userinfo) << "Reading e-folding";
                    readEfold(buildDir(op.getEfolding()), op.getEfolding_a());
                }
            }

            /*
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             * %%% General Head Boundary %%%
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             */
            if (op.isKGHBFromFile()) {
                LOG(userinfo) << "Reading General Head Boundary";
                readGHB_elevation_conductivity(
                        buildDir(op.getKGHBDir())); // should be placed before reading recharge, heads, elevation
            }

            /*
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             * %%% elevation and initial heads %%%
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             */
            LOG(userinfo) << "Reading elevation";
            readElevation(buildDir(op.getElevation()));

            // read either initial head (default) or equilibrium water table depth from file, if available
            if (op.isInitialHeadFromFile()) {
                LOG(userinfo) << "Reading initial head";
                readInitialHeads((buildDir(op.getInitialHeadsDir())));
            } else if (op.isEqWTDFromFile()) {
                LOG(userinfo) << "Reading equal water table depth";
                readEqWTD(buildDir(op.getEqWTD())); // requires elevation to be set
            }

            /*
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             * %%% read hydraulic conductivity of nodes %%%
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             */
            if (op.isKFromFile()) {
                LOG(userinfo) << "Reading hydraulic conductivity (applying to all layers)";
                readConductivity(buildDir(op.getLithology()));
            }

            /*
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             * %%% read recharge, rivers, lakes, wetlands %%%
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             */
            LOG(userinfo) << "Reading groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            if (op.isKRiverFromFile()) {
                LOG(userinfo) << "Reading river conductance";
                readRiverConductance(buildDir(op.getKRiver()));
            } else {
                LOG(userinfo) << "Reading river properties (elevation, length, width, depth), calculating conductance";
                readBlueCells(buildDir(op.getRiverElevation()),
                              calculateRiverStage(buildDir(op.getRiverExtent())));
            }

            LOG(userinfo) << "Reading lakes and wetlands"; // should be placed after reading rivers
            readLakesAndWetlands(buildDir(op.getGlobalLakes()),
                                 buildDir(op.getGlobalWetlands()),
                                 buildDir(op.getLocalLakes()),
                                 buildDir(op.getLocalWetlands()));

            LOG(userinfo) << "Adding river at nodes without surface water body";
            // should be placed after readBlueCells and readLakesAndWetlands
            addRiverWhereNoSurfaceWaterBody(op.getSWBElevationFactor(), op.getRiverConductivity());

            /*
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             * %%% read data for variable density %%%
             * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             */
            if (op.isEffectivePorosityFromFile()) { // needs to be placed before reading zetas
                LOG(userinfo) << "Reading effective porosity";
                readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
            }

            if (op.getStressPeriodVariableDensity()[0]) {
                if (op.isInitialZetasAsArray()) { // needs to be placed after reading effective porosity
                    LOG(userinfo) << "Reading initial heights of " << op.getDensityZones().size() - 1 <<
                                  " active zeta surfaces from file"; // requires elevation to be set
                    readInitialZetas(op.getNumberOfLayers(), op.getDensityZones().size(),
                                     buildDir(op.getInitialZetas()), op.getInitialZetas_a());
                } else {
                    LOG(userinfo) << "Set initial zetas to default (bottom of nodes)";
                    setDefaultZetas(op.getDensityZones().size());
                }
            }
        }
    };
}
}
