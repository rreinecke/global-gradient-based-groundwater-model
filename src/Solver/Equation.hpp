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

#ifndef GLOBAL_FLOW_EQUATION_HPP
#define GLOBAL_FLOW_EQUATION_HPP

//This disables boundary checks!
//#define EIGEN_NO_DEBUG

#include <stdio.h>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <cfenv>
#include <future>
#include <iostream>
#include <unordered_set>

#include "../../lib/Eigen/Sparse"

#include "../../lib/Eigen/Core"
//#include "../../lib/Eigen/PardisoSupport"
#include "../Model/Node.hpp"
#include "../Simulation/Options.hpp"
#include "Numerics.hpp"

namespace GlobalFlow {
    namespace Solver {
        using namespace boost::units;
        using namespace Eigen;

	using pr_t = double; //change here if other precision should be used e.g. long double
	using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;
	using large_num = unsigned long int;
	using long_vector = Matrix<pr_t, Dynamic, 1>;

/**
 * @class Equation The internal finite difference equation
 * Should only be accessed through the stepper
 */
        class Equation {
        public:
            Equation(NodeVector nodes, Simulation::Options options);

            ~Equation();

            /**
             * Solve the current iteration step
             */
            void solve();

            /**
             * Solve Zeta Surface Equation
             */
            void solve_zetas();

            /**
             * @return The number of iterations groundwater flow solution
             */
            int getItter();

            /**
             * @return The number of iterations for variable density solution
             */
            int getItter_zetas();

            /**
             * @return The current residual error
             */
            double getError();

            //No copy and copy assign for Equations
            Equation(const Equation &) = delete;

            Equation &
            operator=(const Equation &) = delete;

            /**
             * Helper to write out current residuals
             * @param os
             * @param eq
             * @return
             */
            friend std::ostream &
            operator<<(std::ostream &os, Equation &eq) {
                IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
                os << eq.getResults().format(CleanFmt);
                return os;
            };

            long_vector getResults() {
                return this->x;
            }

            long_vector getRHS(){
                return this->b;
            }

            long_vector getResults_zetas() {
                return this->x_zetas;
            }

            long_vector getRHS_zetas(){
                return this->b_zetas;
            }

            /**
             * Set simulation step settings
             */
            void updateIsSteadyState(bool is_steadyState) {
                isSteadyState = is_steadyState;
                std::for_each(nodes->begin(),
                              nodes->end(),
                              [is_steadyState](std::unique_ptr<Model::NodeInterface> const &node) {
                                  node->updateIsSteadyState(is_steadyState);
                              });
            }

            void updateIsDensityVariable(bool is_densityVariable) {
                isDensityVariable = is_densityVariable;
                std::for_each(nodes->begin(),
                              nodes->end(),
                              [is_densityVariable](std::unique_ptr<Model::NodeInterface> const &node) {
                                  node->updateIsDensityVariable(is_densityVariable);
                              });
            }

            /**
             * Set the correct stepsize (default is DAY) //
             * @param mod
             */
            void updateStepSize(double stepSize) {
                std::for_each(nodes->begin(),
                              nodes->end(),
                              [stepSize](std::unique_ptr<Model::NodeInterface> const &node) {
                    node->updateStepSize(stepSize);
                });
            }

            typedef typename Eigen::Matrix<pr_t, -1, 1, 0, -1, 1>::Scalar Scalar;
            typedef Matrix<Scalar, Dynamic, 1> VectorType;

            VectorType& getResiduals() const {
                return cg.getResiduals();
            }

            VectorType& getResiduals_zetas() const {
                return cg_zetas.getResiduals();
            }

            void updateClosingCrit(double crit) { cg.setTolerance(crit); }

            //void updateAllowedMaxHeadChange(double head){
            //    maxAllowedHeadChange = head;
            //}

            void updateMaxInnerItter(int iter){
                max_inner_iterations = iter;
                cg.setMaxIterations(iter);
            }

            /**
             * @note resests dampening object and counters
             */
            void enableDamping() {
                isAdaptiveDamping = true;
            }

    private:
        bool initalized = false;

        large_num numberOfNodesPerLayer;

        large_num numberOfLayers;

        long numberOfNodesTotal;

        long numberOfActiveZetas;

        /**
         * _var_ only used if disabling of cells is required
         */
        NodeVector nodes;

        long_vector x;
        long_vector b;
        SparseMatrix<pr_t> A;

        long_vector x_zetas;
        long_vector b_zetas;
        SparseMatrix<pr_t> A_zetas;

        long_vector x_zetas_t0;
        long_vector zetaChanges;
        long_vector oldZetaChanges;

        const Simulation::Options options;

        bool isAdaptiveDamping{true};
        AdaptiveDamping adaptiveDamping;

        long MAX_OUTER_ITERATIONS{0};
        pr_t RCLOSE_HEAD{0};
        pr_t RCLOSE_ZETA{0};
        Index threads;
        long max_inner_iterations{0};
        long max_inner_iterations_zetas{0};

        //From current run
        long __itter{0};
        long __itter_zetas{0};
        double __error{0};

        double currentMaxHeadChange{0};
        double maxAllowedHeadChange{0.01};

        double currentMaxZetaChange{0};
        double maxAllowedZetaChange{0.01};

        double dampMin{0.01};
        double dampMax{0.01};

        std::unordered_map<large_num, std::unordered_map<large_num,long long>> nodeID_to_zetaID_to_rowID;

        ConjugateGradient<SparseMatrix<pr_t>, Lower | Upper, IncompleteLUT<SparseMatrix<pr_t>::Scalar>> cg;

        ConjugateGradient<SparseMatrix<pr_t>, Lower | Upper, IncompleteLUT<SparseMatrix<pr_t>::Scalar>> cg_zetas;

        bool isSteadyState = false;

        bool isGNC{false};

        bool isDensityVariable{false};

        int numberOfZones{0};

        int maxRefinement{1};

        /**
         * Update the matrix for the current iteration
         */

        void inline updateEquation();

        void inline preconditionA();

        void inline updateEquation_zetas(const int layer);

        void inline addToA(std::unique_ptr<Model::NodeInterface> const &node);

        void inline addToA_zetas(std::unique_ptr<Model::NodeInterface> const &node, int zetaID);

        void inline prepareEquation_zetas(const int layer);

        bool inline isHeadChangeGreater();

        /**
         * Update heads in inner iteration
         */
        void inline updateHeadAndHeadChange();

        /**
         * Update zetas in inner iteration
         */
        void inline updateZetaIter(int layer);

        /**
         * Update zone change
         */
        void inline updateZoneChange();

        /**
        * Calculate the ghost node correction flow budget
        */
        void inline updateGNCBudget();

        /**
         * Calculate the final budget
         */
        void inline updateBudget();

        /**
         * Update head change of previous time step
         */
        void inline updateHeadChangeTZero();

        /**
         * Update head of previous time step
         */
        void inline updateHeadTZero();

        /**
         * Write the final zeta surface heights to the nodes
         */
        void inline updateZetas(const int layer);

        /**
         * Write the final zeta surface heights to the nodes
         */
        void inline clipZetas();

        void inline setZetasIter();

        void inline setZetasTZero();

        void inline checkAllZetaSlopes();

        void inline adjustZetaHeights();

        void inline setUnconvergedZetasToZetas_TZero(int layer);

        };
}
}

#endif
