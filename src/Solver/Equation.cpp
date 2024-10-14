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
    this->vdfStepsPerHeadStep = options.getVDFStepsPerHeadStep();

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
    LOG(numerics) << "Destroying equation\n" << std::endl;
}

/**
 * @brief Add entries (conductances to neighbours and storage change within node) to matrix A
 * @param node interface for the current model node within the grid
 */
void inline
Equation::addToA(const std::unique_ptr<Model::NodeInterface>& node) {
    large_num nodeID = node->getID();
    for (const auto &[nodeID_neig, conductance] : node->getMatrixEntries()) {
        A.coeffRef(long(nodeID), long(nodeID_neig)) = conductance.value();
    }
}

/**
 * @brief Add entries (conductances to neighbours and storage change within node) to matrix A
 * @param node interface for the current model node within the grid
 * @param zetaID identifier of the density surface
 */
void inline
Equation::addToA_zetas(std::unique_ptr<Model::NodeInterface> const &node, large_num zetaID) {
    large_num nodeID = node->getID();
    for (const auto &[nodeID_neig, zoneConductance]: nodes->at(nodeID)->getMatrixEntries(zetaID)) {
        A_zetas.coeffRef(nodeID_zetaID_locID[nodeID][zetaID],
                         nodeID_zetaID_locID[nodeID_neig][zetaID]) = zoneConductance.value();
    }
}

/**
 * @brief Update the groundwater flow equation:
 * A - the matrix, b - external flows & storage changes, x - the groundwater heads
 */
void inline
Equation::updateEquation() {
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
    cg.compute(A);
    if (cg.info() != Eigen::Success) {
        LOG(numerics) << "Fail in preconditioning. Perhaps A is asymmetric or a row is all zero!";
        throw "Fail in preconditioning. Perhaps A is asymmetric or a row is all zero!";
    }
}

/**
 * @brief Update the variable density equation:
 * A_zetas - the matrix, b_zetas - external flows & storage changes, x_zetas - the density surface heights
 * @param layer aquifer layer number (increases with depth)
 */
void inline
Equation::updateEquation_zetas(const int& layer) {
#pragma omp parallel for if(numberOfActiveZetas > threads) schedule(dynamic, (numberOfNodesPerLayer/threads)) num_threads(threads) default(none)
    for (long rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        auto zetaID = rowID_to_zetaID[rowID];
        for (const auto &[nodeID_neig, matrixEntry]: nodes->at(nodeID)->getMatrixEntries(zetaID)) {
            A_zetas.coeffRef(nodeID_zetaID_locID[nodeID][zetaID],
                             nodeID_zetaID_locID[nodeID_neig][zetaID]) = matrixEntry.value();
        }
        x_zetas(rowID) = nodes->at(nodeID)->getZeta(zetaID).value();
        b_zetas(rowID) = nodes->at(nodeID)->getRHS(zetaID).value();
    }
    //LOG(debug) << "A_zetas.block:\n" << A_zetas.block(0,0,10,10); // startRow, startCol, numRows, numCol
    //LOG(debug) << "b_zetas.block:\n" << b_zetas.block(0,0,10,1); // startRow, startCol, numRows, numCol
    //LOG(debug) << "x_zetas.block:\n" << x_zetas.block(0,0,10,1); // startRow, startCol, numRows, numCol

    if (A_zetas.size() != 0) {
        A_zetas.makeCompressed();
        cg_zetas.compute(A_zetas);
        if (cg_zetas.info() != Eigen::Success) {
            // https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html
            LOG(userinfo) << "Fail in preconditioning. Perhaps A_zetas is asymmetric or a row is all zero!";
            throw "Fail in preconditioning. Perhaps A_zetas is asymmetric or a row is all zero!";
        }
    }
}

/**
 * @brief Update head change and head after each outer iteration
 */
void inline
Equation::updateHeadAndHeadChange() {
#pragma omp parallel for num_threads(threads) default(none)
    for (long rowID = 0; rowID < numberOfNodesTotal; ++rowID) {
        nodes->at(rowID)->setHeadAndHeadChange(headChanges[rowID] * si::meter);
    }
}

/**
 * @brief Update density surface heights after each outer iteration
 */
void inline
Equation::updateZetas() {
#pragma omp parallel for num_threads(threads) default(none)
    for (long rowID = 0; rowID < rowID_to_nodeID.size(); ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        auto zetaID = rowID_to_zetaID[rowID];
        nodes->at(nodeID)->setZeta_direct(zetaID, x_zetas[rowID] * si::meter);
    }
}

/**
 * @brief After new heads were found: save the new heads as previous heads (are used in next time step)
 */
void inline
Equation::updateHeadChangeTZero() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->updateHeadChange_TZero();
    }
}

/**
 * @brief After new density surfaces were found: save the new surfaces as previous surfaces (are used in next time step)
 */
void inline
Equation::updateHeadTZero() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->updateHead_TZero();
    }
}

/**
 * @brief Apply limits to the front density surface after groundwater heads were found and before new density
 * surface heights are computed
 * @
 */
void inline
Equation::clipFrontZeta() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->clipFrontZeta();
    }
}

/**
 * @brief Perform density surface adjustments (e.g. surface movements between layers, activating of surface in
 * neighboring nodes)
 */
void inline
Equation::adjustZetaHeights() {
    //LOG(numerics) << "    Apply zeta limits (after iteration)";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->limitZetas();
    }

    //LOG(numerics) << "    Save zone change without horizontal tip/toe movement";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveCurrentZoneChange();
    }

    //LOG(numerics) << "    Vertical zeta movement";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->zetaMovementBetweenLayers();
    }

    //LOG(numerics) << "    Horizontal zeta movement";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->horizontalZetaMovement();
    }

    //LOG(numerics) << "    Clipping inner zetas";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->clipInnerZetas();
    }

    //LOG(numerics) << "    Preventing zeta locking";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->preventZetaLocking();
    }

    //LOG(numerics) << "    Correct crossing zetas";
# pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->correctCrossingZetas();
    }

    //LOG(numerics) << "    Check zeta order and whether front and back are in correct position";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->checkZetas();
    }

    //LOG(numerics) << "    Apply zeta limits (after adjustment) and set Zetas_TZero for next step";
    bool setZetasTZero{true};
#pragma omp parallel for num_threads(threads) default(none) shared(setZetasTZero)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->limitZetas(setZetasTZero);
    }

    //LOG(numerics) << "    Save zone change without horizontal tip/toe movement";
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveVDFMassBalance();
    }
}

/**
 * @brief Update the mass balance of groundwater flow after the groundwater flow equation was solved
 */
void inline
Equation::updateBudget() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveMassBalance();
    }
}


void inline
Equation::resetVDFBudget() {
#pragma omp parallel for num_threads(threads) default(none)
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->resetVDFBudget();
    }
}

/**
 * @brief Solve flow equation and (if variable density is activated) solve variable density equation
 */
void
Equation::solve() {
    updateEquation();
    LOG(numerics) << "Initialized A, x and b";
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxAllowedHeadChange, x);

    double oldMaxHeadChange{0};
    int innerIterAddon{0};
    long outerIteration{0};
    long innerIterations{0};
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
        LOG(numerics) << "  Outer iteration: " << outerIteration << ", max absolute head change: " << currentMaxHeadChange;
        /*if (std::abs(x.maxCoeff()) > std::abs(x.minCoeff())) {
            LOG(numerics) << "  Max absolute head: " << x.maxCoeff();
        } else {
            LOG(numerics) << "  Max absolute head: " << x.minCoeff();
        }*/

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
            break;
        }

        if(currentMaxHeadChange == oldMaxHeadChange){
            //The head change is really the same -> increase inner iterations
            innerIterAddon += 10;
            cg.setMaxIterations(max_inner_iterations + innerIterAddon);
        }
        oldMaxHeadChange = currentMaxHeadChange;

        updateEquation();
        outerIteration++;
    }

    if (outerIteration == MAX_OUTER_ITERATIONS_HEAD) {
        LOG(userinfo) << "Fail in solving groundwater flow equation with max iterations";
        LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
        LOG(numerics) << "|Residual|_l2: " << cg.error();
        std::cerr << "Fail in solving matrix with max iterations\n";
    }

    __itter = outerIteration;
    __error = cg.error_inf();

    /**
     * ###################################
     * # Solve Variable Density Equation #
     * ###################################
     */
     if(isDensityVariable) {
         __itter_zetas = 0;
         LOG(numerics) << "Finding zeta surface heights";
         resetVDFBudget();
         for (int vdfStep = 1; vdfStep <= vdfStepsPerHeadStep; ++vdfStep) {
             LOG(numerics) << "  Sub-step " << vdfStep << " of " << vdfStepsPerHeadStep;

             // Clipping top zeta to current groundwater level
             clipFrontZeta();
             std::unordered_map<int, large_num> numberOfActiveZetas_TZero;
             for (int layer = 0; layer < numberOfLayers; layer++) {
                 numberOfActiveZetas_TZero[layer] = rowID_to_nodeID.size();
                 prepareSolveZetas(layer);
                 //int changeOfActiveZetas = int(numberOfActiveZetas - numberOfActiveZetas_TZero[layer]);
                 //LOG(numerics) << "    Number of active zetas on layer " << layer << ":  " << numberOfActiveZetas
                 //               << " (changed by: " << changeOfActiveZetas << ")";
                 if (A_zetas.size() == 0) { // if matrix empty: continue with next layer
                     LOG(debug) << "A_zetas is empty on layer " << layer;
                     continue;
                 }
                 solveZetas(layer);
             }
             LOG(numerics) << "    Adjusting zeta heights (after zeta height convergence)";
             adjustZetaHeights();
         }
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

/**
 * Solve density surface equation
 * @param layer aquifer layer number (increases with depth)
 * @param isAdditionalStep info whether surface heights are solved using smaller additional time steps
 */
void
Equation::solveZetas(const int& layer){
    int outerIteration{0};
    double currentMaxZetaChange{0.0};
    int minorChangeCount{0};

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
        currentMaxZetaChange = std::max(std::abs(zetaChanges.maxCoeff()), std::abs(zetaChanges.minCoeff()));
        LOG(numerics) << "    Outer iteration: " << outerIteration << ", max absolute zeta change: " << currentMaxZetaChange;
        /*if (std::abs(x_zetas.maxCoeff()) > std::abs(x_zetas.minCoeff())) {
            LOG(numerics) << "    Max absolute zeta: " << x_zetas.maxCoeff();
        } else {
            LOG(numerics) << "    Max absolute zeta: " << x_zetas.minCoeff();
        }*/
        if (currentMaxZetaChange < maxAllowedZetaChange) {
            minorChangeCount++;
            if (minorChangeCount >= 2) {
                LOG(numerics) << "    Reached zeta change convergence";
                break;
            }
        } else {
            minorChangeCount = 0;
        }
    } // end of outer iteration loop

    if (outerIteration >= MAX_OUTER_ITERATIONS_ZETA) {
        LOG(userinfo) << "Fail in solving variable density equation with max iterations";
        LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg_zetas.error_inf();
        LOG(numerics) << "|Residual|_l2: " << cg_zetas.error();
    }

    __itter_zetas += outerIteration;
}

/**
 * @brief Prepare density surface for variable density solver
 * @param layer aquifer layer number (increases with depth)
 */
void inline
Equation::prepareSolveZetas(const int& layer) {
    rowID_to_nodeID.clear();
    rowID_to_zetaID.clear();
    nodeID_zetaID_locID.clear();
    large_num offset = layer * numberOfNodesPerLayer;
    long id{0};
    std::unordered_map<large_num, long> zetaID_to_locID;
    // finding nodes with active/inactive surface
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer + offset; ++nodeID) {
        zetaID_to_locID.clear();
        for (large_num zetaID = 1; zetaID < numberOfZones; zetaID++) {
            if (nodes->at(nodeID)->isZetaTZeroActive(zetaID)) {
                rowID_to_nodeID.insert(std::pair<long, large_num>(id, nodeID));
                rowID_to_zetaID.insert(std::pair<long, large_num>(id, zetaID));
                zetaID_to_locID.insert(std::pair<large_num, long>(zetaID, id)); // building map for nodeID_zetaID_locID
                ++id;
            }
        }
        if (!zetaID_to_locID.empty()) {
            nodeID_zetaID_locID.insert(std::pair<large_num, std::unordered_map<large_num,long>>(nodeID, zetaID_to_locID));
        }
    }
    numberOfActiveZetas = long(rowID_to_nodeID.size());
    // LOG(debug) << "numberOfActiveZetas " << numberOfActiveZetas;
    Eigen::SparseMatrix<pr_t> empty_A_zetas(numberOfActiveZetas, numberOfActiveZetas);
    A_zetas = std::move(empty_A_zetas);
    int numberOfEntries = (int) 4 + 1; // + 1 for this node
    A_zetas.reserve(long_vector::Constant(numberOfActiveZetas, numberOfEntries));
    long_vector empty_b_zetas(numberOfActiveZetas);
    b_zetas = std::move(empty_b_zetas);
    long_vector empty_x_zetas(numberOfActiveZetas);
    x_zetas = std::move(empty_x_zetas);
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