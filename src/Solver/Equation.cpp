#include "Equation.hpp"

namespace GlobalFlow {
    namespace Solver {

Equation::Equation(NodeVector nodes, Simulation::Options options) : options(options) {
    this->numberOfNodesPerLayer = options.getNumberOfNodesPerLayer();
    this->numberOfLayers = options.getNumberOfLayers();
    this->numberOfNodesTotal = numberOfNodesPerLayer * numberOfLayers;
    LOG(userinfo) << "Setting up Equation for " << numberOfNodesPerLayer << " nodes"
                  << " on " << numberOfLayers << " layer(s) (in total " << numberOfNodesTotal << " nodes)";

    this->MAX_OUTER_ITERATIONS_HEAD = options.getMaxOuterIterationsHead();
    this->MAX_OUTER_ITERATIONS_ZETA = options.getMaxOuterIterationsZeta();

    this->RCLOSE_HEAD = options.getConverganceCriteriaHead();
    this->RCLOSE_ZETA = options.getConverganceCriteriaZeta();
    this->maxAllowedHeadChange = options.getMaxHeadChange();
    this->isAdaptiveDamping = options.isDampingEnabled();
    this->dampMin = options.getMinDamp();
    this->dampMax = options.getMaxDamp();
    this->threads = options.getThreads();
    this->max_inner_iterations = options.getMaxInnerIterations();
    this->nodes = std::move(nodes);

    this->numberOfZones = options.getDensityZones().size();
    this->maxAllowedZetaChange = options.getMaxZetaChange();

    //set inner iterations
    cg.setMaxIterations(max_inner_iterations);
    cg.setTolerance(RCLOSE_HEAD);
    cg_zetas.setMaxIterations(max_inner_iterations);
    cg_zetas.setTolerance(RCLOSE_ZETA);

    Eigen::SparseMatrix<pr_t> sparseMatrix(numberOfNodesTotal, numberOfNodesTotal);
    A = std::move(sparseMatrix);
    int numberOfEntries = (int) 4 + 2 + 1; // +2 for top/down, + 1 for this node
    A.reserve(long_vector::Constant(numberOfNodesTotal, numberOfEntries));
    long_vector __b(numberOfNodesTotal);
    b = std::move(__b);
    long_vector __x(numberOfNodesTotal);
    x = std::move(__x);
}

Equation::~Equation() {
    LOG(debug) << "Destroying equation\n" << std::endl;
}

void inline
Equation::addToA(std::unique_ptr<Model::NodeInterface> const &node) {
    large_num nodeID = node->getID();
    for (const auto &[nodeID_neig, conductance] : node->getMatrixEntries()) {
        A.coeffRef(long(nodeID), long(nodeID_neig)) = conductance.value();
    }
}

void inline
Equation::addToA_zetas(std::unique_ptr<Model::NodeInterface> const &node, int zetaID) {
    large_num nodeID = node->getID();
    for (const auto &[nodeID_neig, zoneConductance]: nodes->at(nodeID)->getMatrixEntries(zetaID)) {
        A_zetas.coeffRef(nodeID_zetaID_rowID[nodeID][zetaID],
                         nodeID_zetaID_rowID[nodeID_neig][zetaID]) = zoneConductance.value();
    }
}

void inline
Equation::updateEquation() {
    //LOG(debug) << "Updating equation";
#ifdef EIGEN_HAS_OPENMP
    Eigen::initParallel();
#endif
#pragma omp parallel for schedule(dynamic, (numberOfNodesTotal/(threads * 4))) num_threads(threads) default(none)
    for (int nodeID = 0; nodeID < numberOfNodesTotal; ++nodeID) {
        addToA(nodes->at(nodeID));
        x(nodeID) = nodes->at(nodeID)->getHead().value(); // cannot just use x, since damping might alter head values
        b(nodeID) = nodes->at(nodeID)->getRHS().value();
    }

    if (!A.isCompressed()) { A.makeCompressed(); }
    //LOG(debug) << "Compressed A";
    cg.compute(A);
    //LOG(debug) << "Computed conjugate gradients for A";
    if (cg.info() != Success) {
        LOG(numerics) << "Fail in preconditioning matrix";
        throw "Fail in preconditioning matrix";
    }
}


void inline
Equation::updateEquation_zetas(const int layer) {

//#pragma omp parallel for if(numberOfActiveZetas > threads) schedule(dynamic, (numberOfNodesPerLayer/threads)) num_threads(threads) default(none)
    for (large_num rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        auto zetaID = rowID_to_zetaID[rowID];
        addToA_zetas(nodes->at(nodeID), zetaID);
        x_zetas(rowID) = nodes->at(nodeID)->getZeta(zetaID).value();
        b_zetas(rowID) = nodes->at(nodeID)->getRHS(zetaID).value();
    }

    //LOG(debug) << "A_zetas.block:\n" << A_zetas.block(0,0,numberOfActiveZetas,numberOfActiveZetas); // startRow, startCol, numRows, numCol
    //LOG(debug) << "b_zetas.block:\n" << b_zetas.block(0,0,numberOfActiveZetas,1); // startRow, startCol, numRows, numCol
    //LOG(debug) << "x_zetas.block:\n" << x_zetas.block(0,0,numberOfActiveZetas,1); // startRow, startCol, numRows, numCol

    //LOG(numerics) << "Preconditioning matrix before iteration (zetas)";
    if (A_zetas.size() != 0) {
        //LOG(numerics) << "Compressing Matrix (zetas)";
        A_zetas.makeCompressed();
        cg_zetas.compute(A_zetas);
        if (cg_zetas.info() != Success) {
            LOG(userinfo) << "Fail in preconditioning matrix (zetas)";
            throw "Fail in preconditioning matrix (zetas)";
        }
    }
}

void inline
Equation::updateHeadAndHeadChange() {
// no parallel here
    for (long rowID = 0; rowID < numberOfNodesTotal; rowID++) {
        nodes->at(rowID)->setHeadAndHeadChange(headChanges[rowID] * si::meter);
    }
}


void inline
Equation::updateZetas() {
// todo parallel here
    for (large_num rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        auto zetaID = rowID_to_zetaID[rowID];
        nodes->at(nodeID)->setZeta(zetaID, x_zetas[long(rowID)] * si::meter);
    }
}


void inline
Equation::updateHeadChangeTZero() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->updateHeadChange_TZero();
    }
}

void inline
Equation::updateHeadTZero() {
#pragma omp parallel for num_threads(threads) default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->updateHead_TZero();
        }
    }

void inline
Equation::clipFrontZeta() {
//#pragma omp parallel for num_threads(threads) default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->clipFrontZeta();
        }
}


void inline
Equation::adjustZetaHeights() {
    LOG(debug) << "Vertical zeta movement";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->zetaMovementBetweenLayers();
    }

    LOG(debug) << "Horizontal zeta movement";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->horizontalZetaMovement();
    }

    LOG(debug) << "Clipping inner zetas";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->clipInnerZetas();
    }

    LOG(debug) << "Preventing zeta locking";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->preventZetaLocking();
    }

    LOG(debug) << "Correct crossing zetas";
# pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->correctCrossingZetas();
    }

    LOG(debug) << "Check zeta order and whether front and back are in correct position";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->checkZetas();
    }

    LOG(debug) << "Set Zetas_TZero for next step";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num nodeID = 0; nodeID < numberOfNodesTotal; ++nodeID) {
        auto zetas = nodes->at(nodeID)->getZetas();
        nodes->at(nodeID)->setZetas_TZero(zetas);
    }
}

void inline
Equation::updateZoneChange() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveZoneChange();
    }
}

void inline
Equation::updateBudget() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveMassBalance();
    }
}

/**
 * Solve Equation
 *
 */
void
Equation::solve() {
    updateEquation();
    LOG(debug) << "Initialized A, x and b";
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxAllowedHeadChange, x);

    double oldMaxHeadChange{0};
    int innerIterAddon{0};
    long outerIteration{0};
    long innerIterations{0};
    bool headFail{false};
    char smallHeadChangeCounter{0};
    bool headConverged{false};
    double currentMaxHead{0};
    while (outerIteration < MAX_OUTER_ITERATIONS_HEAD) {
        x = cg.solveWithGuess(b, x); // solving inner iterations
        headChanges = adaptiveDamping.getChanges(getResiduals(), x, isAdaptiveDamping);
        innerIterations = cg.iterations();
        //LOG(numerics) << "Inner iterations: " << innerIterations;
        if (innerIterations == 0 and outerIteration == 0) {
            LOG(numerics) << "Convergence criterion met without iterations.";
            break;
        }
        updateHeadAndHeadChange(); // needs to be before "head change convergence"
        currentMaxHeadChange = std::max(std::abs(headChanges.maxCoeff()), std::abs(headChanges.minCoeff()));
        LOG(numerics) << "Max Absolute Head Change: " << currentMaxHeadChange;
        currentMaxHead = std::max(std::abs(x.maxCoeff()), std::abs(x.minCoeff()));
        LOG(numerics) << "Max Absolute Head: " << currentMaxHead;

        /**
         * @brief head change convergence
         */
        if (currentMaxHeadChange > maxAllowedHeadChange) { //convergence is not reached
            smallHeadChangeCounter = 0; //reset counter
        } else { //converged with head criterion
            if (headConverged and smallHeadChangeCounter == 0) {
                LOG(numerics) << "Conditional convergence - check mass balance";
                break;
            } else {
                headConverged = true;
                smallHeadChangeCounter++;
                if (smallHeadChangeCounter >= 2) {
                    LOG(numerics) << "Reached head change convergence";
                    //LOG(debug) << "x.block:\n" << x.block(0,0,30,1); // startRow, startCol, numRows, numCol
                    break;
                }
            }
        }

        /**
         * @brief residual norm convergence
         */
        if (cg.info() == Success and outerIteration != 0) {
            LOG(numerics) << "cg solver success";
            //LOG(debug) << "x:\n" << x;
            break;
        }

        if(currentMaxHeadChange == oldMaxHeadChange){
            //The head change is really the same -> increase inner iterations
            innerIterAddon += 10;
            cg.setMaxIterations(max_inner_iterations + innerIterAddon);
        }
        oldMaxHeadChange = currentMaxHeadChange;

        updateEquation();
        //LOG(debug) << "Updated A, x and b";
        outerIteration++;
    }

    if (outerIteration == MAX_OUTER_ITERATIONS_HEAD) {
        std::cerr << "Fail in solving matrix with max iterations\n";
        LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
        LOG(numerics) << "|Residual|_l2: " << cg.error();
    }

    __itter = outerIteration;
    __error = cg.error_inf();

    /**
     * ###############################
     * # Solve Zeta Surface Equation #
     * ###############################
     */
     if(isDensityVariable) {
         __itter_zetas = 0;
         // Clipping top zeta to current groundwater level
         clipFrontZeta();
         LOG(debug) << "Check zeta order and whether front and back are in correct position";
#pragma omp parallel for num_threads(threads) default(none)
         for (large_num k = 0; k < numberOfNodesTotal; ++k) {
             nodes->at(k)->checkZetas();
         }

         for (int layer = 0; layer < numberOfLayers; layer++) {
             LOG(numerics) << "Finding zeta surface heights in layer " << layer;
             prepareEquation_zetas(layer);
             if (A_zetas.size() == 0) { continue; } // if matrix empty: continue with next iteration
             solve_zetas(layer, false); // isAdditionalStep=false
         }
         updateZoneChange(); // needs to be before adjustZetaHeights to get zone change without horizontal tip/toe movement
         LOG(numerics) << "Adjusting zeta heights (after zeta height convergence)";
         adjustZetaHeights();
     }

    /**
    * ###############################
    * # Update budgets #
    * ###############################
    */

    LOG(numerics) << "Updating head change and head of previous time step";
    updateHeadChangeTZero();
    updateHeadTZero();
    LOG(numerics) << "Updating budget";
    updateBudget();
}

void inline
Equation::resetZetas(int layer) {
    large_num offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for num_threads(threads) default(none) shared(offset)
    for (large_num rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        auto zetaID = rowID_to_zetaID[rowID];
        auto zetaTZero = nodes->at(nodeID)->getZeta_TZero(zetaID);
        nodes->at(nodeID)->setZeta(zetaID, zetaTZero);
    }
}

void inline
Equation::updateZetaTimeStep(int layer, double additionalSteps) {
#pragma omp parallel for num_threads(threads) default(none) shared(additionalSteps)
    for (large_num rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        auto zetaID = rowID_to_zetaID[rowID];
        nodes->at(nodeID)->updateZetaStepSize(nodes->at(nodeID)->getStepSize().value() / additionalSteps);
    }
}

void inline
Equation::alignZetaTimeStep(int layer) {
    large_num offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for num_threads(threads) default(none) shared(offset)
    for (large_num rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        nodes->at(nodeID)->alignZetaStepSize();
    }
}




/**
 * Solve Zeta Surface Equation
 */
void
Equation::solve_zetas(int layer, bool isAdditionalStep){
    int outerIteration{0};
    double maxZetaChange{0};
    int minorChangeCount{0};
    double curMaxAllowedZetaChange{0};
    if (isAdditionalStep) {
        curMaxAllowedZetaChange = maxAllowedZetaChange / stepSize;
    } else {
        curMaxAllowedZetaChange = maxAllowedZetaChange;
    }
    while (outerIteration < MAX_OUTER_ITERATIONS_ZETA) {
        outerIteration++;
        updateEquation_zetas(layer);
        x_zetas_t0 = x_zetas;
        x_zetas = cg_zetas.solveWithGuess(b_zetas, x_zetas); // solving inner iterations
        zetaChanges = x_zetas - x_zetas_t0;
        updateZetas();

        /**
         * @brief zeta change convergence
         */
        maxZetaChange = std::max(std::abs(zetaChanges.maxCoeff()), std::abs(zetaChanges.minCoeff()));
        LOG(numerics) << "Outer iteration: " << outerIteration << ", max zeta change: " << maxZetaChange;
        if (maxZetaChange < curMaxAllowedZetaChange) {
            minorChangeCount++;
            if (minorChangeCount >= 2) {
                LOG(numerics) << "Reached zeta change convergence";
                break;
            }
            if (cg_zetas.info() == Success) {
                LOG(numerics) << "cg_zetas solver success";
                //LOG(debug) << "x_zetas.block:\n" << x_zetas.block(0,0,numberOfActiveZetas,1); // startRow, startCol, numRows, numCol
                break;
            }
        } else {
            minorChangeCount = 0;
        }

        if (outerIteration == MAX_OUTER_ITERATIONS_ZETA) { //  or std::abs(currentMaxZetaChange) > 100*maxZetaChange
            resetZetas(layer); // reset zetas to values before non-converging outer iteration
            if (!isAdditionalStep) {
                LOG(numerics) << "Restart solving zetas with reduced step size";
                outerIteration = 0; // set outer iteration to 0
                updateZetaTimeStep(layer, stepSize); // change to smaller zeta time step
                solve_zetas(layer, true); // calling solve_zetas with smaller time step, with isAdditionalStep=true
                alignZetaTimeStep(layer); // re-align zeta time step with head time step
                break;
            } else {
                // Deactivate zetas at non-converging node: set effective porosity to 0, and zetas to bottom
                for (long rowID = 0; rowID < zetaChanges.size(); ++rowID) {
                    if (std::abs(zetaChanges(rowID)) > curMaxAllowedZetaChange) {
                        LOG(debug) << "Deactivating zetas for nodeID = " << rowID_to_nodeID[rowID];
                        nodes->at(rowID_to_nodeID[rowID])->deactivateZetas();
                    }
                }
                // rerun solve_zeta without deactivated nodes
                prepareEquation_zetas(layer);
                solve_zetas(layer, false); // isAdditionalStep=false
            }
        }
    } // end of outer iteration loop
    __itter_zetas += outerIteration;
}


void inline
Equation::prepareEquation_zetas(const int layer) {
    auto numberOfActiveZetas_TZero = rowID_to_nodeID.size();
    rowID_to_nodeID.clear();
    large_num offset = layer * numberOfNodesPerLayer;
    numberOfActiveZetas = 0;
    int count_no_new_nodes{0};
    long rowID{0};

    // finding nodes with active/inactive interfaces
    for (int zetaID = 1; zetaID < numberOfZones; zetaID++) {
        for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer + offset; nodeID++) {
            if (nodes->at(nodeID)->isZetaTZeroActive(zetaID)) {
                rowID_to_nodeID[rowID] = nodeID;
                rowID_to_zetaID[rowID] = zetaID;
                nodeID_zetaID_rowID[nodeID][zetaID] = rowID;  // for addToA_zetas()
                ++rowID;
            }
        }
    }
    numberOfActiveZetas = long(rowID_to_nodeID.size());
    LOG(debug) << "Number of active zetas on layer " << layer << ":  " << numberOfActiveZetas;

    int changeOfActiveZetas = int(numberOfActiveZetas - numberOfActiveZetas_TZero);
    LOG(debug) << "Number of active zetas changed by: " << changeOfActiveZetas;

    Eigen::SparseMatrix<pr_t> __A_zetas(numberOfActiveZetas, numberOfActiveZetas);
    A_zetas = std::move(__A_zetas);
    int numberOfEntries = (int) 4 + 1; // + 1 for this node
    A_zetas.reserve(long_vector::Constant(numberOfActiveZetas, numberOfEntries));
    long_vector __b_zetas(numberOfActiveZetas);
    b_zetas = std::move(__b_zetas);
    long_vector __x_zetas(numberOfActiveZetas);
    x_zetas = std::move(__x_zetas);
}


int
Equation::getItter() {
    return __itter;
}

int Equation::getItter_zetas(){
    return __itter_zetas;
}

double
Equation::getError() {
    return __error;
}

}
}//ns