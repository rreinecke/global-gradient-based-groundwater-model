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
            Equation(large_num numberOfNodesPerLayer, NodeVector nodes, Simulation::Options options);

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
             * @return The number of iterations
             */
            int getItter();

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
             * Toogle the steady-state in all nodes
             * @return
             */
            bool toggleSteadyState() {
                SteadyState = !SteadyState;
                bool state = SteadyState;
                std::for_each(nodes->begin(),
                              nodes->end(),
                              [state](std::unique_ptr<Model::NodeInterface> const &node) {
                                  node->toggleSteadyState(state);
                              });
                return SteadyState;
            }

            /**
             * Set the correct stepsize (default is DAY) // todo remove if not needed
             * @param mod
             */
            void updateStepModifier(double mod) {
                std::for_each(nodes->begin(),
                              nodes->end(),
                              [mod](std::unique_ptr<Model::NodeInterface> const &node) { node->updateStepModifier(mod); });
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

        large_num numberOfNodesTotal;

        int initialHead;

        /**
         * _var_ only used if disabling of cells is required
         */
        NodeVector nodes;

        long_vector x;
        long_vector _x_;
        long_vector b;
        long_vector _b_;
        SparseMatrix<pr_t> A;
        SparseMatrix<pr_t> _A_;

        long_vector x_zetas;
        long_vector b_zetas;
        SparseMatrix<pr_t> A_zetas;

        const Simulation::Options options;

        bool isAdaptiveDamping{true};
        AdaptiveDamping adaptiveDamping;
        AdaptiveDamping adaptiveDamping_zetas;

        int IITER{0};//FIXME this is used as outer iterations
        pr_t RCLOSE{0};
	    int inner_iterations{0};

        //From current run
        int __itter{0};
        double __error{0};

        bool isCached{false};
        bool isCached_zetas{false};

        double maxHeadChange{0.01};
        double maxZetaChange{0.01};
        double dampMin{0.01};
        double dampMax{0.01};

        bool disable_dry_cells{false};
        //Maybe rename me :D
        std::unordered_set<large_num> disabled_nodes;
        std::unordered_set<large_num> inactive_nodes;
        //Real -> Current
        std::unordered_map<large_num, long long> index_mapping; // Question: do we need an index_mapping_zetas or is this enough?

        bool dry_have_changed{true};
        bool inactive_have_changed{true};

        template<typename Set>
        bool set_compare(Set const &lhs, Set const &rhs) {
            return lhs.size() == rhs.size()
                   && std::equal(lhs.begin(), lhs.end(), rhs.begin());
        }

        ConjugateGradient<SparseMatrix<pr_t>, Lower | Upper, IncompleteLUT<SparseMatrix<pr_t>::Scalar>> cg;

        ConjugateGradient<SparseMatrix<pr_t>, Lower | Upper, IncompleteLUT<SparseMatrix<pr_t>::Scalar>> cg_zetas;

        BiCGSTAB<SparseMatrix<pr_t>, IncompleteLUT<SparseMatrix<pr_t>::Scalar>> bicgstab;

        BiCGSTAB<SparseMatrix<pr_t>, IncompleteLUT<SparseMatrix<pr_t>::Scalar>> bicgstab_zetas;
        //Used for NWT
        bool nwt{false};

        bool vdf{false};

        int numberOfZones{0};
        /**
         * Helper for updating the matrix
         * @param node
         * @param cached
         */
        void addToA(std::unique_ptr<Model::NodeInterface> const &node, bool cached);

        void addToA_zeta(large_num nodeIter, large_num iterOffset, int localZetaID, bool cached);

        /**
         * Update the matrix for the current iteration
         */
        void inline updateMatrix();

        void inline updateMatrix_zetas(large_num iterOffset, int localZetaID);

        /**
         * Reallocate matrix and vectors based on dried nodes
         * @bug This is currently missing re-enabling of deactivated nodes!
         * Re-enable if:
         * 1) head in cell below needs to be higher than threshold
         * 2) head in one of 4 neighbours higher than threshold
         */
        void inline reallocateMatrix();

        /**
         * Run the preconditioner for heads
         */
        void inline preconditioner();

        /**
         * Run the preconditioner for zetas
         */
        void inline preconditioner_zetas();

        /**
         * Update heads in inner iteration
         */
        void inline updateIntermediateHeads();

        /**
         * Update zetas in inner iteration
         */
        void inline updateIntermediateZetas(large_num iterOffset, int localZetaID);

        /**
         * Calculate the final budget
         */
        void inline updateBudget();

        /**
         * Write the final head to the nodes
         */
        void inline updateFinalHeads();

        /**
         * Write the final zeta surface heights to the nodes
         */
        void inline updateZetaChanges(int localZetaID);

        /**
         * Write the final zeta surface heights to the nodes
         */
        void inline updateTopZetasToHeads();

        void inline setZetasPosInNodes();

        void inline checkAllZetaSlopes(int localZetaID);

        void inline updateZetasAfterEquation();

        void inline adjustZetaHeights();

        bool SteadyState = false;
        //Only for testing purposes
        bool simpleHead = true;

        };
}
}

#endif
