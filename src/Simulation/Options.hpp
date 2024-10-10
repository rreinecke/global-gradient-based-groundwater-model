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
            std::vector<std::string> EFOLDING_a;
            std::vector<std::string> EQUAL_WATER_TABLE_DEPTH_a;
            std::vector<std::string> RIVER_ELEVATION_a;
            std::vector<std::string> INITIAL_ZETAS_a;
            std::vector<std::string> VERTICAL_SIZE_a;

            std::string EFOLDING_DIR{""};
            std::string INITIAL_ZETAS_DIR{""};
            std::string VERTICAL_SIZE_DIR{""};

            std::string RIVER_ELEVATION_FILE{""};
            std::string ELEVATION_FILE{""};
            std::string EQUAL_WATER_TABLE_DEPTH_FILE{""};
            std::string RECHARGE_FILE{""};
            std::string ZONES_SOURCES_FILE{""};
            std::string LITHOLOGY_FILE{""};
            std::string RIVER_FILE{""};
            std::string GLOBAL_LAKES_FILE{""};
            std::string GLOBAL_WETLANDS_FILE{""};
            std::string LOCAL_LAKES_FILE{""};
            std::string LOCAL_WETLANDS_FILE{""};
            std::string K_FILE{""};
            std::string RIVER_K_FILE{""};
            std::string GHB_FILE{""};
            std::string SS_FILE{""};
            std::string SY_FILE{""};
            std::string INITIAL_HEAD_FILE{""};
            std::string EFFECTIVE_POROSITY_FILE{""};

            //++Special mappings++//
            std::string SPATID_ARCID{""};

            //++General configuration++//
            // model config
            std::vector<bool> STRESS_PERIOD_STEADY_STATE{true};
            std::vector<int> STRESS_PERIOD_STEPS{0};
            std::vector<std::string> STRESS_PERIOD_STEP_SIZES{""};
            std::vector<bool> STRESS_PERIOD_VARIABLE_DENSITY{false};
            std::string NODES;
            unsigned long int NUMBER_OF_NODES_PER_LAYER{0};
            long Y_RANGE{0};
            long X_RANGE{0};
            bool IS_GLOBAL{false};
            double RESOLUTION_IN_DEGREE{0.0};
            double EDGE_LENGTH_LEFT_RIGHT{0.0};
            double EDGE_LENGTH_FRONT_BACK{0.0};
            int LAYERS{0};
            std::vector<bool> CONFINED{};
            bool USE_EFOLDING{false};
            std::string DEFAULT_BOUNDARY_CONDITION{"GeneralHeadBoundary"};
            bool SENSITIVITY{false};

            // vdf config
            std::vector<double> DENSITY_ZONES{1000.0};
            double MAX_TIP_SLOPE{0.2};
            double MAX_TOE_SLOPE{0.2};
            double MIN_DEPTH_FACTOR{0.1};
            double SLOPE_ADJ_FACTOR{0.1};
            double VDF_LOCK{0.001};
            int VDF_STEPS_PER_HEAD_STEP{1};

            // numerics
            int THREADS{0};
            std::string SOLVER{"PCG"};
            int MAX_OUTER_ITERATIONS_HEAD{0};
            int MAX_OUTER_ITERATIONS_ZETA{0};
            int MAX_INNER_ITERATIONS{0};
            double RCLOSE_HEAD{0.1};
            double RCLOSE_ZETA{0.1};
            double MAX_HEAD_CHANGE{0.01};
            double MAX_ZETA_CHANGE{0.01};
            bool DAMPING{false};
            double MIN_DAMP{0.01};
            double MAX_DAMP{0.5};

            // input: data config
            bool effective_porosity_from_file{false};
            bool efold_as_array{false};
            bool eq_wtd_from_file{false};
            bool ghb_from_file{false};
            bool initial_head_from_file{false};
            bool initial_zetas_as_array{false};
            bool k_from_file{false};
            bool k_river_from_file{false};
            bool specificstorage_from_file{false};
            bool specificyield_from_file{false};
            bool vertical_size_as_array{false};
            bool zones_sources_from_file{false};

            // input: default data
            double EFFECTIVE_POROSITY{0.0};
            double INITIAL_HEAD{0.0};
            bool ADAPTIVE_STEP_SIZE{false};
            std::vector<double> K{0.001};
            double GHB_K{0.1};
            std::vector<double> VERTICAL_SIZES{100};
            std::vector<double> ANISOTROPY{10};
            double RIVER_CONDUCTIVITY{10.0};
            double SWB_ELEVATION_FACTOR{0.8};
            double SPECIFIC_YIELD{0.15};
            double SPECIFIC_STORAGE{0.000015};
            int SOURCE_ZONE_GHB{0};
            int SOURCE_ZONE_RECHARGE{0};
            double MIN_GHB_K{0};
            double MIN_K{0};
            double MIN_EFFECTIVE_POROSITY{0};
            double MIN_VERTICAL_SIZE{0};

        public:

            enum BoundaryCondition {
                GENERAL_HEAD_BOUNDARY,
                GENERAL_HEAD_NEIGHBOUR,
                STATIC_HEAD_SEA_LEVEL,
                NONE
            };

            std::vector<bool> getStressPeriodSteadyState() { return STRESS_PERIOD_STEADY_STATE; }

            std::vector<int> getStressPeriodSteps() { return STRESS_PERIOD_STEPS; }

            std::vector<std::string> getStressPeriodStepSizes() { return STRESS_PERIOD_STEP_SIZES; }

            std::vector<bool> getStressPeriodVariableDensity() {return STRESS_PERIOD_VARIABLE_DENSITY; }

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

            bool isGHBFromFile() { return ghb_from_file; }

            bool isSpecificStorageFile() { return specificstorage_from_file; }

            bool isSpecificYieldFile() { return specificyield_from_file; }

            bool isKRiverFromFile() { return k_river_from_file; }

            bool isVerticalSizeAsArray() { return vertical_size_as_array; }

            bool isEqWTDFromFile() { return eq_wtd_from_file;}

            bool isInitialHeadFromFile() { return initial_head_from_file;}

            bool isEffectivePorosityFromFile() { return effective_porosity_from_file;}

            bool isZonesSourcesFromFile() { return zones_sources_from_file;}

            bool isInitialZetasAsArray(){ return initial_zetas_as_array; }

            std::string getKDir() { return K_FILE; }

            std::string getKRiver() { return RIVER_K_FILE; }

            std::string getGHBDir() { return GHB_FILE; }

            std::string getSSDir() { return SS_FILE; }

            std::string getSYDir() { return SY_FILE; }

            std::string getInitialHeadsDir() {return INITIAL_HEAD_FILE;}

            std::string getEffectivePorosityDir() {return EFFECTIVE_POROSITY_FILE;}

            bool isGlobal() { return IS_GLOBAL; }

            int getMaxInnerIterations() { return MAX_INNER_ITERATIONS; }

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

            bool isDensityVariable() { return DENSITY_VARIABLE; }

            std::vector<double>
            getDensityZones() {
                return DENSITY_ZONES;
            }

            double getEffectivePorosity() { return EFFECTIVE_POROSITY; }

            double getMaxTipSlope() { return MAX_TIP_SLOPE; }

            double getMaxToeSlope() { return MAX_TOE_SLOPE; }

            double getMinDepthFactor() { return MIN_DEPTH_FACTOR; }

            double
            getSlopeAdjFactor() {
                return SLOPE_ADJ_FACTOR;
            }

            double
            getVDFLock() {
                return VDF_LOCK;
            }

            int getVDFStepsPerHeadStep() { return VDF_STEPS_PER_HEAD_STEP; }

            int getSourceZoneGHB() { return SOURCE_ZONE_GHB; }

            int getSourceZoneRecharge() { return SOURCE_ZONE_RECHARGE; }

            int
            getMaxOuterIterationsHead() {
                return MAX_OUTER_ITERATIONS_HEAD;
            }

            int
            getMaxOuterIterationsZeta() {
                return MAX_OUTER_ITERATIONS_ZETA;
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
                return ELEVATION_FILE;
            }

            std::string
            getEfolding() {
                return EFOLDING_DIR;
            }

            std::string
            getEqWTD() {
                return EQUAL_WATER_TABLE_DEPTH_FILE;
            }

            std::string
            getInitialZetasDir() {
                return INITIAL_ZETAS_DIR;
            }

            std::string getRiverElevation() {
                return RIVER_ELEVATION_FILE;
            }

            std::vector<std::string>
            getVerticalSize_a() {
                return VERTICAL_SIZE_a;
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
                return RECHARGE_FILE;
            }

            std::string
            getVerticalSizeDir() {
                return VERTICAL_SIZE_DIR;
            }

            std::string
            getLithology() {
                return LITHOLOGY_FILE;
            }

            std::string
            getRiverExtent() {
                return RIVER_FILE;
            }

            std::string
            getGlobalLakes() {
                return GLOBAL_LAKES_FILE;
            }

            std::string
            getGlobalWetlands() {
                return GLOBAL_WETLANDS_FILE;
            }

            std::string
            getLocalLakes() {
                return LOCAL_LAKES_FILE;
            }

            std::string
            getLocalWetlands() {
                return LOCAL_WETLANDS_FILE;
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

            std::vector<double>
            getVerticalSizes() {
                return VERTICAL_SIZES;
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

            double
            getMinGHBK() {
                return MIN_GHB_K;
            }

            double
            getMinK() {
                return MIN_K;
            }

            double
            getMinEffectivePorosity() {
                return MIN_EFFECTIVE_POROSITY;
            }

            double
            getMinVerticalSize() {
                return MIN_VERTICAL_SIZE;
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