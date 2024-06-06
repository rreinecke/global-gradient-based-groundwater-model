#include "Equation.hpp"

namespace GlobalFlow {
namespace Solver {

Equation::Equation(NodeVector nodes, Simulation::Options options) : options(options) {
    this->numberOfNodesPerLayer = options.getNumberOfNodesPerLayer();
    this->numberOfLayers = options.getNumberOfLayers();
    this->numberOfNodesTotal = numberOfNodesPerLayer * numberOfLayers;
    LOG(userinfo) << "Setting up Equation for " << numberOfNodesPerLayer << " nodes"
                  << " on " << numberOfLayers << " layer(s) (in total " << numberOfNodesTotal << " nodes)";

    this->MAX_OUTER_ITERATIONS = options.getMaxIterations();
    this->RCLOSE_HEAD = options.getConverganceCriteriaHead();
    this->RCLOSE_ZETA = options.getConverganceCriteriaZeta();
    this->maxAllowedHeadChange = options.getMaxHeadChange();
    this->isAdaptiveDamping = options.isDampingEnabled();
    this->dampMin = options.getMinDamp();
    this->dampMax = options.getMaxDamp();
    this->isGNC = options.isGridRefined(); // gnc = Ghost Node Correction
    this->threads = options.getThreads();
    this->max_inner_iterations = options.getInnerItter();
    this->maxRefinement = options.getMaxRefinement();
    this->nodes = std::move(nodes);

    this->numberOfZones = options.getDensityZones().size();
    this->maxAllowedZetaChange = options.getMaxZetaChange();

    //set inner iterations
    cg.setMaxIterations(max_inner_iterations);
    cg.setTolerance(RCLOSE_HEAD);
    max_inner_iterations_zetas = max_inner_iterations * 2;
    cg_zetas.setMaxIterations(max_inner_iterations_zetas);
    cg_zetas.setTolerance(RCLOSE_ZETA);

    Eigen::SparseMatrix<pr_t> sparseMatrix(numberOfNodesTotal, numberOfNodesTotal);
    A = std::move(sparseMatrix);
    int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 2 + 1; // +2 for top/down, + 1 for this node
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
        A_zetas.coeffRef(nodeID_to_zetaID_to_rowID[nodeID][zetaID],
                         nodeID_to_zetaID_to_rowID[nodeID_neig][zetaID]) = zoneConductance.value();
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
}

void inline
Equation::preconditionA() {
    if (!A.isCompressed()) { A.makeCompressed(); }
    //LOG(debug) << "Compressed A";
    cg.compute(A);
    //LOG(debug) << "Computed conjugate gradients for A";
    if (cg.info() != Success) {
        LOG(numerics) << "Fail in preconditioning matrix";
        throw "Fail in preconditioning matrix";
    }
    //LOG(debug) << "Equation is ready for solver";
}

void inline
Equation::updateEquation_zetas(const int layer) {
    large_num offset = layer * numberOfNodesPerLayer;

#pragma omp parallel for if(numberOfActiveZetas > threads) schedule(dynamic, (numberOfNodesPerLayer/threads)) default(none) shared(offset)
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer + offset; ++nodeID) {
        for (int zetaID = 1; zetaID < numberOfZones; zetaID++) {
            if (nodeID_to_zetaID_to_rowID[nodeID][zetaID] != -1) {
                addToA_zetas(nodes->at(nodeID), zetaID);
                x_zetas(nodeID_to_zetaID_to_rowID[nodeID][zetaID]) = nodes->at(nodeID)->getZetaIter(zetaID).value();
                b_zetas(nodeID_to_zetaID_to_rowID[nodeID][zetaID]) = nodes->at(nodeID)->getRHS(zetaID).value();
            }
        }
    }

    LOG(debug) << "A_zetas.block:\n" << A_zetas.block(0,0,numberOfActiveZetas,numberOfActiveZetas); // startRow, startCol, numRows, numCol
    LOG(debug) << "b_zetas.block:\n" << b_zetas.block(0,0,numberOfActiveZetas,1); // startRow, startCol, numRows, numCol
    LOG(debug) << "x_zetas.block:\n" << x_zetas.block(0,0,numberOfActiveZetas,1); // startRow, startCol, numRows, numCol

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
    //LOG(debug) << "Equation is ready for solver (zetas)";
}

void inline
Equation::updateHeadAndHeadChange() {
    long_vector changes = adaptiveDamping.getChanges(getResiduals(), x, isAdaptiveDamping);
// no parallel here
    for (long rowID = 0; rowID < numberOfNodesTotal; rowID++) {
        nodes->at(rowID)->setHeadAndHeadChange(changes[rowID] * si::meter);
    }
}

void inline
Equation::updateZetaIter(int layer) {
    int numAboveMaxZetaChange{0};
    auto offset = layer * numberOfNodesPerLayer;
    large_num spatID;
// no parallel here
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer + offset; nodeID++) {
        for (int zetaID = 1; zetaID < numberOfZones; zetaID++) {
            if (nodeID_to_zetaID_to_rowID[nodeID][zetaID] != -1) {
                nodes->at(nodeID)->setZetaIter(zetaID,
                                               zetaChanges[nodeID_to_zetaID_to_rowID[nodeID][zetaID]] * si::meter);
                if (std::abs(zetaChanges[nodeID_to_zetaID_to_rowID[nodeID][zetaID]]) > maxAllowedZetaChange) {
                    ++numAboveMaxZetaChange;
                    spatID = nodes->at(nodeID)->getSpatID();
                    //LOG(numerics) << "zetaChange (" << nodeID << ", zetaID: " << zetaID << "): " << zetaChanges[nodeID_to_zetaID_to_rowID[nodeID][zetaID]];
                    //LOG(numerics) << "x_zetas (" << nodeID << ", " << zetaID << "): " << x_zetas(nodeID_to_zetaID_to_rowID[nodeID][zetaID]);
                    //LOG(numerics) << "zetaIter (" << nodeID << ", " << zetaID+1 << "): " << nodes->at(nodeID)->getZetaIter(zetaID+1).value();
                }
            }
        }
    }

    if (std::abs(zetaChanges.maxCoeff()) > std::abs(zetaChanges.minCoeff())) {
        currentMaxZetaChange = zetaChanges.maxCoeff();
    } else {
        currentMaxZetaChange = zetaChanges.minCoeff();
    }
    LOG(numerics) << "Zeta Change larger than allowed value: " << numAboveMaxZetaChange
                  << " times. MAX Zeta Change: " << currentMaxZetaChange;
    if (numAboveMaxZetaChange == 1) {
        LOG(numerics) << "SpatID of that one node with too high zeta change:" << spatID;
    }
}


void inline
Equation::updateZetas(const int layer) {
    auto offset = layer * numberOfNodesPerLayer;
# pragma omp parallel for default(none) shared(offset)
    for (large_num k = offset; k < numberOfNodesPerLayer + offset; ++k) {
        for (int zetaID = 1; zetaID < numberOfZones; zetaID++) {
            auto zeta = nodes->at(k)->getZetaIter(zetaID);
            nodes->at(k)->setZeta(zetaID, zeta);
        }
    }
}


void inline
Equation::updateHeadChangeTZero() {
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->updateHeadChange_TZero();
    }
}

void inline
Equation::updateHeadTZero() {
#pragma omp parallel for default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->updateHead_TZero();
        }
    }

void inline
Equation::clipZetas() {
#pragma omp parallel for default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->clipZetas();
        }
}

void inline
Equation::setZetasTZero() {
#pragma omp parallel for default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->setZetas_TZero();
        }
    }

void inline
Equation::setZetasIter() {
#pragma omp parallel for default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            auto zetas = nodes->at(k)->getZetas_TZero();
            nodes->at(k)->setZetas_Iter(zetas);
        }
    }

void inline
Equation::adjustZetaHeights() {
    LOG(debug) << "Calculating vertical zeta movement";
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->zetaMovementBetweenLayers();
    }

    LOG(debug) << "Calculating horizontal zeta movement";
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->horizontalZetaMovement();
    }

    LOG(debug) << "Clipping inner zetas";
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->clipInnerZetas();
    }

    LOG(debug) << "Preventing zeta locking";
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->preventZetaLocking();
    }

    LOG(debug) << "Correct corssing zetas";
# pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->correctCrossingZetas();
    }

    LOG(debug) << "Check zeta order and whether front and back are in correct position";
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->checkZetas();
    }
}

void inline
Equation::updateZoneChange() {
#pragma omp parallel for default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveZoneChange();
    }
}

void inline
Equation::updateGNCBudget() {
#pragma omp parallel for default(none)
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->saveGNCMassBalance();
        }
    }

void inline
Equation::updateBudget() {
#pragma omp parallel for default(none)
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
    preconditionA();
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxAllowedHeadChange, x); // among other: set x_t0 = x (initial)

    double oldMaxHeadChange{0};
    int innerIterAddon{0};
    long outerIteration{0};
    long innerIterations{0};
    bool headFail{false};
    char smallHeadChangeCounter{0};
    bool headConverged{false};
    while (outerIteration < MAX_OUTER_ITERATIONS) {
        //LOG(debug) << "A.block:\n" << A.block(0,0,50,50); // startRow, startCol, numRows, numCol
        //LOG(debug) << "b.block:\n" << b.block(0,0,50,1); // startRow, startCol, numRows, numCol
        //LOG(debug) << "x.block:\n" << x.block(0,0,50,1); // startRow, startCol, numRows, numCol
        x = cg.solveWithGuess(b, x); // solving inner iterations
        innerIterations = cg.iterations();
        //LOG(numerics) << "Inner iterations: " << innerIterations;
        if (innerIterations == 0 and outerIteration == 0) {
            LOG(numerics) << "Convergence criterion met without iterations.";
            break;
        }
        updateHeadAndHeadChange(); // needs to be before "head change convergence"

        /**
         * @brief head change convergence
         */
        if (isHeadChangeGreater()) { //convergence is not reached
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
                    break;
                }
            }
        }

        /**
         * @brief residual norm convergence
         */
        if (cg.info() == Success and outerIteration != 0) {
            LOG(numerics) << "cg solver success";
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
        preconditionA();
        outerIteration++;
    }

    if (outerIteration == MAX_OUTER_ITERATIONS) {
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
         solve_zetas();
     }

    /**
    * ###############################
    * # Update budgets #
    * ###############################
    */
    if(isGNC) {
        updateGNCBudget();
    }

    LOG(numerics) << "Updating head change and head of previous time step";
    updateHeadChangeTZero();
    updateHeadTZero();
    LOG(numerics) << "Updating budget";
    updateBudget();
}

bool inline
Equation::isHeadChangeGreater(){
    currentMaxHeadChange = 0;
    double absHeadMax{0.0};
    double headSum{0.0};
    double headMean{0.0};
    double countAbsHeadAbove10k{0.0};
    double changeAtNode{0.0};
    double absHeadAtNode{0.0};
    large_num nodeID_headMax{0};
    large_num nodeID_changeMax{0};
//#pragma omp parallel for
    for (large_num nodeID = 0; nodeID < numberOfNodesTotal; nodeID++) {
        changeAtNode = std::abs(
                nodes->at(nodeID)->getProperties().get<quantity<Model::Meter>, Model::HeadChange>().value());
        nodeID_changeMax = (changeAtNode > currentMaxHeadChange) ? nodeID : nodeID_changeMax;
        currentMaxHeadChange = (changeAtNode > currentMaxHeadChange) ? changeAtNode : currentMaxHeadChange;
        absHeadAtNode = std::abs(nodes->at(nodeID)->getHead().value());
        nodeID_headMax = (absHeadAtNode > absHeadMax) ? nodeID : nodeID_headMax;
        absHeadMax = (absHeadAtNode > absHeadMax) ? absHeadAtNode : absHeadMax;
        headSum += nodes->at(nodeID)->getHead().value();
        if (absHeadAtNode > 10000) {
            //LOG(debug) << "nodeID: " << nodeID << "-> head = " << nodes->at(nodeID)->getHead().value() << ", x = " << x[nodeID];
            countAbsHeadAbove10k++;}
    }
    headMean = headSum / numberOfNodesTotal;
    //LOG(debug) << "MAX Absolute Head: " << absHeadMax << " (nodeID: " << nodeID_headMax << "), "
    //              << "Absolute head > 10,000 at " << countAbsHeadAbove10k << " nodes";
    //LOG(debug) << "MEAN Head: " << headMean;
    LOG(numerics) << "MAX Absolute Head Change: " << currentMaxHeadChange << " (nodeID: " << nodeID_changeMax << ")";
    return currentMaxHeadChange > maxAllowedHeadChange;
}

void inline
Equation::setUnconvergedZetasToZetas_TZero(int layer) {
    large_num offset = layer * numberOfNodesPerLayer;
    double zetaTZero{0.0};
    for (int zetaID = 1; zetaID < numberOfZones; ++zetaID) {
        for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer + offset; ++nodeID) {
            if (nodeID_to_zetaID_to_rowID[nodeID][zetaID] != -1) {
                if (std::abs(zetaChanges[nodeID_to_zetaID_to_rowID[nodeID][zetaID]]) > maxAllowedZetaChange){
                    zetaTZero = nodes->at(nodeID)->getZetaTZero(zetaID).value();
                    nodes->at(nodeID)->setZetaIter(zetaID, zetaTZero * si::meter);
                }
            }
        }
    }
}

/**
 * Solve Zeta Surface Equation
 */
void
Equation::solve_zetas(){
    __itter_zetas = 0;
    // If unconfined: clipping top zeta to current groundwater level
    clipZetas();
    setZetasTZero();
    setZetasIter();
    for (int layer = 0; layer < numberOfLayers; layer++) {
        LOG(numerics) << "Finding zeta surface heights in layer " << layer;

        prepareEquation_zetas(layer);

        if (A_zetas.size() == 0) { continue; } // if matrix empty continue with next iteration

        int outerIteration{0};
        long innerIteration{0};
        int innerIterAddon{0};
        char smallZetaChanges{0};
        bool zetaConverged{false};
        int outOfBoundsCount{0};
        while (outerIteration < MAX_OUTER_ITERATIONS) {
            LOG(numerics) << "Outer iteration: " << outerIteration;
            x_zetas_t0 = x_zetas;
            //Solve inner iterations
            x_zetas = cg_zetas.solveWithGuess(b_zetas, x_zetas);
            zetaChanges = x_zetas - x_zetas_t0;

            updateZetaIter(layer);
            innerIteration = cg_zetas.iterations();
            /*if (innerIteration == 0 and outerIteration == 0) {
                LOG(numerics) << "Convergence criterion for zeta (=" << RCLOSE_ZETA << ") too small - no iterations";
                break;
            }*/

            /**
         * @brief residual norm convergence
         */
            if (cg_zetas.info() == Success and outerIteration != 0) {
                LOG(numerics) << "cg_zetas solver success";
                break;
            }


            /**
             * @brief zeta change convergence // todo make function of this
             */
            if (std::abs(currentMaxZetaChange) > maxAllowedZetaChange) {
                //convergence is not reached -> reset counter
                smallZetaChanges = 0;
            } else {
                //itter converged with zeta criterion
                if (zetaConverged and smallZetaChanges == 0) {
                    LOG(numerics) << "Conditional convergence - check mass balance (zetas)";
                    break;
                } else {
                    zetaConverged = true;
                    smallZetaChanges++;

                    if (smallZetaChanges >= 2) {
                        LOG(numerics) << "Reached zeta change convergence";
                        break;
                    }
                }
            }

            updateEquation_zetas(layer);
            outerIteration++;
        } // end of outer iteration loop

        if (outerIteration == MAX_OUTER_ITERATIONS) {
            LOG(userinfo) << "Fail in solving zeta matrix with max iteration (layer: " << layer << ")";
            LOG(userinfo) << "Setting unconverged zetas to their value before iteration (layer: " << layer << ")";
            setUnconvergedZetasToZetas_TZero(layer);
        }

        __itter_zetas += outerIteration;
        // %%%%%%%%%%%%%%%%%%%%%%
        // % Update final zetas %
        // %%%%%%%%%%%%%%%%%%%%%%
        updateZetas(layer);
    } // end of layer loop



    updateZoneChange();
    LOG(numerics) << "Adjusting zeta heights (after zeta height convergence)";
    adjustZetaHeights();
}

void inline
Equation::prepareEquation_zetas(const int layer) {
    nodeID_to_zetaID_to_rowID.clear();
    large_num offset = layer * numberOfNodesPerLayer;
    numberOfActiveZetas = 0;
    long rowID{0};
    // finding nodes with active/inactive interfaces
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer + offset; nodeID++) {
        for (int zetaID = 1; zetaID < numberOfZones; zetaID++) {
            if (nodes->at(nodeID)->isZetaTZeroActive(zetaID)) {
                nodeID_to_zetaID_to_rowID[nodeID][zetaID] = rowID;
                ++rowID;
                ++numberOfActiveZetas; // tracking how many have been set inactive at this zetaID
            } else {
                nodeID_to_zetaID_to_rowID[nodeID][zetaID] = -1; // these entries will be ignored ( e.g. in loop filling A_zeta, x_zeta and b_zeta)
            }
        }
    }

    Eigen::SparseMatrix<pr_t> __A_zetas(numberOfActiveZetas, numberOfActiveZetas);
    A_zetas = std::move(__A_zetas);
    int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 1; // + 1 for this node
    A_zetas.reserve(long_vector::Constant(numberOfActiveZetas, numberOfEntries));
    long_vector __b_zetas(numberOfActiveZetas);
    b_zetas = std::move(__b_zetas);
    long_vector __x_zetas(numberOfActiveZetas);
    x_zetas = std::move(__x_zetas);

    updateEquation_zetas(layer);
    //oldZetaChanges = x_zetas * 0.0; // set to 0 at all nodes
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
