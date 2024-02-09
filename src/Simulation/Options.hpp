/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef GLOBAL_FLOW_OPTIONS_HPP
#define GLOBAL_FLOW_OPTIONS_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional/optional.hpp>

namespace GlobalFlow {
    namespace Simulation {

/**
 * @class Options
 * Reads simulation options from a JSON file
 * Defines getters and setters for options
 */
        class Options {

            //++Input data++//
            std::vector<std::string> ELEVATION_a;
            std::vector<std::string> EFOLDING_a;
            std::vector<std::string> EQUAL_WATER_TABLE_DEPTH_a;
            std::vector<std::string> RIVER_ELEVATION_a;
            std::vector<std::string> INITIAL_ZETAS_a;

            std::string ELEVATION{""};
            std::string EFOLDING{""};
            std::string EQUAL_WATER_TABLE_DEPTH{""};
            std::string RIVER_ELEVATION{""};
            std::string INITIAL_ZETAS{""};

            std::string RECHARGE{""};
            std::string ZONES_SOURCES_SINKS_FILE{""};
            std::string PSEUDO_SOURCE_FLOW{""};
            std::string LITHOLOGY{""};
            std::string RIVER{""};
            std::string GLOBAL_LAKES{""};
            std::string GLOBAL_WETLANDS{""};
            std::string LOCAL_LAKES{""};
            std::string LOCAL_WETLANDS{""};
            std::string K_DIR{""};
            std::string RIVER_K{""};
            std::string GHB_K_DIR{""};
            std::string SS_FILE{""};
            std::string SY_FILE{""};
            std::string AQ_DEPTH{""};
            std::string INITIAL_HEAD_FILE{""};
            std::string INITIAL_ZETAS_FILE{""};
            std::string INITIAL_ZONES{""};
            std::string EFFECTIVE_POROSITY_FILE{""};

            //++Special mappings++//
            std::string SPATID_ARCID{""};

            //++General configuration++//
            unsigned long int NUMBER_OF_NODES_PER_LAYER{0};
            long Y_RANGE{0};
            long X_RANGE{0};
            double RESOLUTION_IN_DEGREE{0.0};
            double EDGE_LENGTH_LEFT_RIGHT{0.0};
            double EDGE_LENGTH_FRONT_BACK{0.0};
            int LAYERS{0};
            bool USE_EFOLDING{false};
            int IITER{0};
            int I_ITTER{0};
            double RCLOSE_HEAD{0.1};
            double RCLOSE_ZETA{0.1};

            std::string SOLVER{"PCG"};
            std::vector<int> STEADY_STATE_STRESS_PERIOD_STEPS{0};
            std::vector<int> TRANSIENT_STRESS_PERIOD_STEPS{0};
            std::vector<std::string> STEADY_STATE_STRESS_PERIOD_STEPSIZES{""};
            std::vector<std::string> TRANSIENT_STRESS_PERIOD_STEPSIZES{""};
            std::string NODES{""};
            int THREADS{0};
            bool CACHE{false};
            bool ADAPTIVE_STEP_SIZE{false};
            double INITIAL_HEAD{0.0};
            std::vector<double> K{0.001};
            double GHB_K{0.1};
            double RIVER_CONDUCTIVITY{10.0};
            double SWB_ELEVATION_FACTOR{0.8};
            std::vector<int> AQUIFER_DEPTH{100};
            std::vector<double> ANISOTROPY{10};
            double SPECIFIC_YIELD{0.15};
            double SPECIFIC_STORAGE{0.000015};
            std::string DEFAULT_BOUNDARY_CONDITION{"GeneralHeadBoundary"};
            bool SENSITIVITY{false};
            std::vector<bool> CONFINED{};
            // refinement information
            bool GRID_REFINED{false};
            int MAX_REFINEMENT{1};
            // density information
            bool DENSITY_VARIABLE{false};
            std::vector<double> DENSITY_ZONES{1000.0};
            double EFFECTIVE_POROSITY{0.0};
            double MAX_TIP_SLOPE{0.2};
            double MAX_TOE_SLOPE{0.2};
            double MIN_DEPTH_FACTOR{0.1};
            double SLOPE_ADJ_FACTOR{0.1};
            double VDF_LOCK{0.001};
            std::vector<int> ZONES_SOURCES_SINKS{0};

            bool k_from_file{false};
            bool k_ghb_from_file{false};
            bool specificstorage_from_file{false};
            bool specificyield_from_file{false};
            bool INITIAL_ZETAS_AS_ARRAY{false};
            bool k_river_from_file{false};
            bool aquifer_depth_from_file{false};
            bool eq_wtd_from_file{false};
            bool initial_head_from_file{false};
            bool effective_porosity_from_file{false};
            bool zones_sources_sinks_from_file{false};

            bool IS_GLOBAL{false};
            double MAX_HEAD_CHANGE{0.01};
            double MAX_ZETA_CHANGE{0.01};
            bool DAMPING{false};
            double MIN_DAMP{0.01};
            double MAX_DAMP{0.5};

        public:

            enum BoundaryCondition {
                GENERAL_HEAD_BOUNDARY,
                GENERAL_HEAD_NEIGHBOUR,
                STATIC_HEAD_SEA_LEVEL,
                NONE
            };

            std::vector<int> getSteadyStateStressPeriodSteps() { return STEADY_STATE_STRESS_PERIOD_STEPS; }

            std::vector<int> getTransientStressPeriodSteps() { return TRANSIENT_STRESS_PERIOD_STEPS; }

            std::vector<std::string> getSteadyStateStressPeriodStepsizes() {
                return STEADY_STATE_STRESS_PERIOD_STEPSIZES;
            }

            std::vector<std::string> getTransientStressPeriodStepsizes() {
                return TRANSIENT_STRESS_PERIOD_STEPSIZES;
            }

            void setClosingCritHead(double crit_head) { RCLOSE_HEAD = crit_head; }

            void setClosingCritZeta(double crit_zeta) { RCLOSE_ZETA = crit_zeta; }


            void setDamping(bool set) { DAMPING = set; }

            bool isDampingEnabled() { return DAMPING; }

            bool useEfolding() { return USE_EFOLDING; }

            double getMinDamp() { return MIN_DAMP; }

            double getMaxDamp() { return MAX_DAMP; }

            double getMaxHeadChange() { return MAX_HEAD_CHANGE; }

            double getMaxZetaChange() { return MAX_ZETA_CHANGE; }

            bool isConfined(int layer) { return CONFINED[layer]; }


            std::vector<bool> getConfinements() { return CONFINED; }

            BoundaryCondition getBoundaryCondition() {
                if (DEFAULT_BOUNDARY_CONDITION == "GeneralHeadBoundary") {
                    return BoundaryCondition::GENERAL_HEAD_BOUNDARY;
                }
                if (DEFAULT_BOUNDARY_CONDITION == "GeneralHeadNeighbour") {
                    return BoundaryCondition::GENERAL_HEAD_NEIGHBOUR;
                }
                if (DEFAULT_BOUNDARY_CONDITION == "StaticSeaLevel"){
                    return BoundaryCondition::STATIC_HEAD_SEA_LEVEL;
                }
                return BoundaryCondition::NONE;
            }

            bool isSensitivity() { return SENSITIVITY; }

            bool isKFromFile() { return k_from_file; }

            bool isKGHBFromFile() { return k_ghb_from_file; }

            bool isSpecificStorageFile() { return specificstorage_from_file; }

            bool isSpecificYieldFile() { return specificyield_from_file; }

            bool isKRiverFromFile() { return k_river_from_file; }

            bool isAquiferDepthFile() { return aquifer_depth_from_file; }

            bool isEqWTDFromFile() { return eq_wtd_from_file;}

            bool isInitialHeadFromFile() { return initial_head_from_file;}

            bool isEffectivePorosityFromFile() { return effective_porosity_from_file;}

            bool isZonesSourcesSinksFromFile() { return zones_sources_sinks_from_file;}

            bool isInitialZetasAsArray(){ return INITIAL_ZETAS_AS_ARRAY; }

            std::string getKDir() { return K_DIR; }

            std::string getKRiver() { return RIVER_K; }

            std::string getKGHBDir() { return GHB_K_DIR; }

            std::string getSSDir() { return SS_FILE; }

            std::string getSYDir() { return SY_FILE; }

            std::string getAQDepthDir() { return AQ_DEPTH; }

            std::string getInitialHeadsDir() {return INITIAL_HEAD_FILE;}

            std::string getEffectivePorosityDir() {return EFFECTIVE_POROSITY_FILE;}

            bool isGlobal() { return IS_GLOBAL; }

            int getInnerItter() { return I_ITTER; }

            unsigned long int
            getNumberOfNodesPerLayer() {
                return NUMBER_OF_NODES_PER_LAYER;
            };

            long
            getYRange() {
                return Y_RANGE;
            };

            long
            getXRange() {
                return X_RANGE;
            };

            double
            getResolution() {
                return RESOLUTION_IN_DEGREE;
            }

            double
            getEdgeLengthLeftRight() {
                return EDGE_LENGTH_LEFT_RIGHT;
            };

            double
            getEdgeLengthFrontBack() {
                return EDGE_LENGTH_FRONT_BACK;
            };

            int
            getNumberOfLayers() {
                return LAYERS;
            }

            bool isGridRefined() { return GRID_REFINED; }  // todo could be removed

            int getMaxRefinement() { return MAX_REFINEMENT; }

            bool isDensityVariable() { return DENSITY_VARIABLE; }

            std::vector<double>
            getDensityZones() {
                return DENSITY_ZONES;
            }

            double
            getEffectivePorosity() {
                return EFFECTIVE_POROSITY;
            }

            double
            getMaxTipSlope() {
                return MAX_TIP_SLOPE;
            }

            double
            getMaxToeSlope() {
                return MAX_TOE_SLOPE;
            }

            double
            getMinDepthFactor() {
                return MIN_DEPTH_FACTOR;
            }

            double
            getSlopeAdjFactor() {
                return SLOPE_ADJ_FACTOR;
            }

            double
            getVDFLock() {
                return VDF_LOCK;
            }

            int
            getMaxIterations() {
                return IITER;
            }

            double
            getConverganceCriteriaHead() {
                return RCLOSE_HEAD;
            }

            double
            getConverganceCriteriaZeta() {
                return RCLOSE_ZETA;
            }

            std::string
            getSolverName() {
                return SOLVER;
            }

            std::string
            getNodesDir() {
                return NODES;
            }

            std::string
            getElevation() {
                return ELEVATION;
            }

            std::string
            getEfolding() {
                return EFOLDING;
            }

            std::string
            getEqWTD() {
                return EQUAL_WATER_TABLE_DEPTH;
            }

            std::string
            getInitialZetas() {
                return INITIAL_ZETAS;
            }

            std::string getRiverElevation() {
                return RIVER_ELEVATION;
            }

            std::vector<std::string>
            getElevation_A() {
                return ELEVATION_a;
            }

            std::vector<std::string>
            getEfolding_a() {
                return EFOLDING_a;
            }

            std::vector<std::string>
            getEqWTD_a() {
                return EQUAL_WATER_TABLE_DEPTH_a;
            }

            std::vector<std::string>
            getRiverElevation_a() {
                return RIVER_ELEVATION_a;
            }

            std::vector<std::string>
            getInitialZetas_a() {
                return INITIAL_ZETAS_a;
            }

            std::string
            getRecharge() {
                return RECHARGE;
            }

            std::string
            getZonesOfSourcesAndSinksDir() {
                return ZONES_SOURCES_SINKS_FILE;
            }

            std::string
            getPseudoSourceFlow() {
                return PSEUDO_SOURCE_FLOW;
            }

            std::string
            getLithology() {
                return LITHOLOGY;
            }

            std::string
            getRiverExtent() {
                return RIVER;
            }

            std::string
            getGlobalLakes() {
                return GLOBAL_LAKES;
            }

            std::string
            getGlobalWetlands() {
                return GLOBAL_WETLANDS;
            }

            std::string
            getLocalLakes() {
                return LOCAL_LAKES;
            }

            std::string
            getLocalWetlands() {
                return LOCAL_WETLANDS;
            }

            std::string
            getMapping() {
                return SPATID_ARCID;
            }

            int
            getThreads() {
                return THREADS;
            }

            const bool
            adaptiveStepSizeEnabled() {
                return ADAPTIVE_STEP_SIZE;
            }

            bool
            cacheEnabled() {
                return CACHE;
            }

            double
            getInitialHead() {
                return INITIAL_HEAD;
            }

            std::vector<double>
            getInitialK() {
                return K;
            }

            double
            getGHBConduct() {
                return GHB_K;
            }

            double
            getRiverConductivity() {
                return RIVER_CONDUCTIVITY;
            }

            double getSWBElevationFactor() { return SWB_ELEVATION_FACTOR; }

            std::vector<int>
            getAquiferDepth() {
                return AQUIFER_DEPTH;
            }

            std::vector<double>
            getAnisotropy() {
                return ANISOTROPY;
            }

            double
            getSpecificYield() {
                return SPECIFIC_YIELD;
            }

            double
            getSpecificStorage() {
                return SPECIFIC_STORAGE;
            }

            void
            load(const std::string &filename);

            //FIXME implement node serialization
            void
            save(const std::string &filename);
        };

    }
}//ns
#endif //COVERAGE_OPTIONS_H