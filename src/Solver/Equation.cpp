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
    this->initialHead = options.getInitialHead();
    this->maxAllowedHeadChange = options.getMaxHeadChange();
    this->isAdaptiveDamping = options.isDampingEnabled();
    this->dampMin = options.getMinDamp();
    this->dampMax = options.getMaxDamp();
    this->vdf = options.isDensityVariable(); // vdf = Variable Density Flow
    this->gnc = options.isGridRefined(); // gnc = Ghost Node Correction

    this->max_inner_iterations = options.getInnerItter();
    this->maxRefinement = options.getMaxRefinement();
    this->nodes = nodes;

    //set inner iterations
    cg.setMaxIterations(max_inner_iterations);
    cg.setTolerance(RCLOSE_HEAD);

    long numInactive{0};
    // finding inactive nodes
//#pragma omp parallel for
    large_num rowID{0};
    large_num colID{0};
    for (large_num nodeID = 0; nodeID < numberOfNodesTotal; ++nodeID) {
        if ( nodes->at(nodeID)->getHeadActive()) {
            // active node
            rowID = nodeID - numInactive;
            rowID_to_nodeID[rowID] = nodeID;
            colID = nodeID - numInactive; // todo test this
            neigNodeID_to_colID[nodeID] = colID; // used in addToA(): neigNodeIDs coming bank from getMatrixEntries()
        } else {
            // inactive node
            numInactive++; // tracking how many have been set inactive
        }
    }
    numberOfActiveNodes = numberOfNodesTotal - numInactive;
    LOG(userinfo) << "Number of inactive nodes (K very close to 0): " << numInactive;
    LOG(userinfo) << "Number of active nodes: " << numberOfActiveNodes;

    Eigen::SparseMatrix<pr_t> __A(numberOfActiveNodes, numberOfActiveNodes);
    A = std::move(__A);
    int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 2 + 1; // +2 for top/down, + 1 for this node
    A.reserve(long_vector::Constant(numberOfActiveNodes, numberOfEntries));
    long_vector __b(numberOfActiveNodes);
    b = std::move(__b);
    long_vector __x(numberOfActiveNodes);
    x = std::move(__x);

    if(vdf) {
        LOG(userinfo) << "Simulating variable density flow" << std::endl;
        this->numberOfZones = options.getDensityZones().size();
        this->maxAllowedZetaChange = options.getMaxZetaChange();

#pragma omp parallel for
        for (int i = 0; i < numberOfNodesTotal; ++i) {
            nodes->at(i)->initZetasTZero();
        }

        //set inner iterations
        cg_zetas.setMaxIterations(max_inner_iterations);
        cg_zetas.setTolerance(RCLOSE_ZETA);
    }
}

Equation::~Equation() {
    LOG(debug) << "Destroying equation\n" << std::endl;
}

void inline
Equation::addToA(large_num &rowID) {
    for (const auto &[neigNodeID, conductance] : nodes->at(rowID_to_nodeID[rowID])->getMatrixEntries()) {
        if (neigNodeID_to_colID.find(neigNodeID) != neigNodeID_to_colID.end()) {
            NANChecker(conductance.value(), "Matrix entry");
            //LOG(debug) << "colID = " << neigNodeID_to_colID[neigNodeID] << ", rowID = " << rowID;
            A.coeffRef(long(rowID), long(neigNodeID_to_colID[neigNodeID])) = conductance.value();
        }
    }
}


void inline
Equation::addToA_zeta(large_num iter, large_num offset, int localZetaID) {
    large_num nodeID = iter + offset;
    long long rowID = iter_to_rowID[iter];
    long long colID;

    for (const auto &[nodeID_neig, zoneConductance] : nodes->at(nodeID)->getMatrixEntries(localZetaID)) {
        colID = iter_to_rowID[nodeID_neig - offset]; // todo rename iter_to_rowID (in this line it's something to colID)
        if (colID != -1) {
            NANChecker(zoneConductance.value(), "Matrix entry (zetas)");
            //LOG(debug) << "colID = " << colID << ", rowID = " << rowID << ", iterOffset = " << iterOffset;
            A_zetas.coeffRef(rowID, colID) = zoneConductance.value();
        }
    }
}

void inline
Equation::updateEquation() {
#ifdef EIGEN_HAS_OPENMP
    Eigen::initParallel();
    Index threads = Eigen::nbThreads();
#endif
    int chunkSize = int(A.outerSize()/threads);
    large_num colID{0};
    large_num nodeID{0};
#pragma omp parallel for schedule(dynamic,chunkSize) num_threads(threads)
    for (large_num rowID = 0; rowID < numberOfActiveNodes; rowID++) {
        nodeID = rowID_to_nodeID[rowID];
        //---------------------Left
        //addToA(rowID);
        auto matrixEntries = nodes->at(rowID_to_nodeID[rowID])->getMatrixEntries();
        for (const auto &[neigNodeID, conductance] : matrixEntries) {
            if (neigNodeID_to_colID.find(neigNodeID) != neigNodeID_to_colID.end()) {
                NANChecker(conductance.value(), "Matrix entry");
                colID = neigNodeID_to_colID[neigNodeID];
                A.coeffRef(long(rowID), long(colID)) = conductance.value();
            }
        }
        x(long(rowID)) = nodes->at(nodeID)->getHead().value();
        //---------------------Right
        b(long(rowID)) = nodes->at(nodeID)->getRHS().value();
        NANChecker(b[long(rowID)], "Right hand side");
    }

    //LOG(debug) << "A:\n" << A;
    //LOG(debug) << "b (= rhs):\n" << b << std::endl;

    LOG(numerics) << "Preconditioning matrix before iteration";
    preconditionMatrix();
    LOG(numerics) << "Equation is ready for solver";
}

void inline
Equation::updateEquation_zetas(large_num layer, int localZetaID) {
    large_num numInactive{0};
    iter_to_rowID.clear();

    large_num offset = layer * numberOfNodesPerLayer;
    // finding inactive nodes
//#pragma omp parallel for
    for (large_num iter = 0; iter < numberOfNodesPerLayer; ++iter) {
        if ( nodes->at(iter + offset)->isZetaActive(localZetaID) ) {
            iter_to_rowID[iter] = iter - numInactive;
        } else {
            ++numInactive; // tracking how many have been set inactive
            iter_to_rowID[iter] = -1; // these entries will be ignored ( e.g. in loop filling A_zeta, x_zeta and b_zeta)
        }
    }

    const long numberOfActiveZetas = numberOfNodesPerLayer - numInactive;
    Eigen::SparseMatrix<pr_t> __A_zetas(numberOfActiveZetas, numberOfActiveZetas);
    A_zetas = std::move(__A_zetas);
    //int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 1; // + 1 for this node
    //A_zetas.reserve(long_vector::Constant(numberOfActiveZetas, numberOfEntries));
    long_vector __b_zetas(numberOfActiveZetas);
    b_zetas = std::move(__b_zetas);
    long_vector __x_zetas(numberOfActiveZetas);
    x_zetas = std::move(__x_zetas);

//#pragma omp parallel for
    for (large_num iter2 = 0; iter2 < numberOfNodesPerLayer; ++iter2) {
        auto nodeID = iter2 + offset;
        auto row_id = iter_to_rowID[iter2];
        if (row_id != -1) {
            //---------------------Left: fill matrix A_zeta and initiate x_zetas
            addToA_zeta(iter2, offset, localZetaID);
            x_zetas(row_id) = nodes->at(nodeID)->getZeta(localZetaID).value();
            //---------------------Right
            b_zetas(row_id) = nodes->at(nodeID)->getRHS(localZetaID).value();
        }
    }
    LOG(numerics) << "Preconditioning matrix before iteration (zetas)";
    preconditionMatrix_zetas();
    LOG(debug) << "Equation is ready for solver (zetas)";
}

void inline
Equation::preconditionMatrix() {
    if (!A.isApprox(A.adjoint())) {
        LOG(userinfo) << "Matrix NOT self-adjoint. Preconditioner will fail.";
    }

    if (!A.isCompressed()) {
        LOG(numerics) << "Compressing Matrix";
        A.makeCompressed();
        LOG(numerics) << "... done";
    }
    cg.compute(A);
    if (cg.info() != Success) {
        LOG(numerics) << "Fail in preconditioning matrix";
        throw "Fail in preconditioning matrix";
    }
}

void inline
Equation::preconditionMatrix_zetas() {
    /*if (!A_zetas.isApprox(A_zetas.adjoint())) {
        LOG(userinfo) << "Matrix NOT self-adjoint (zetas). Terminating.";
        throw "Matrix NOT self-adjoint. Terminating (zetas).";
    }*/
    if (A_zetas.size() != 0 and !A_zetas.isCompressed()) {
        LOG(numerics) << "Compressing Matrix (zetas)";
        A_zetas.makeCompressed();
        LOG(numerics) << "... done";
    }
    cg_zetas.compute(A_zetas);
    if (cg_zetas.info() != Success) {
        LOG(userinfo) << "Fail in preconditioning matrix (zetas)";
        throw "Fail in preconditioning matrix (zetas)";
    }
}

void inline
Equation::updateIntermediateHeads() {
    long_vector changes = adaptiveDamping.getChanges(getResiduals(), x, isAdaptiveDamping);
    //LOG(debug) << "changes:\n" << changes;
#pragma omp parallel for
    for (large_num rowID = 0; rowID < numberOfActiveNodes; ++rowID) {
        auto nodeID = rowID_to_nodeID[rowID];
        // set new head (= old head + change) and headChange
        nodes->at(nodeID)->setHeadChange(changes[rowID] * si::meter);
    }
}

void inline
Equation::updateIntermediateZetas(large_num layer, int localZetaID) {
    auto offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for
    for (large_num iter = 0; iter < numberOfNodesPerLayer; ++iter) {
        auto nodeID = iter + offset;
        auto rowID = iter_to_rowID[iter];
        if (rowID != -1) {
            nodes->at(nodeID)->setZetaChange(localZetaID, (double) x_zetas[rowID] * si::meter);  // before setZeta
            nodes->at(nodeID)->setZeta(localZetaID, (double) x_zetas[rowID] * si::meter);
        }
    }
}

void inline
Equation::updateHeadChangeTZero() {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->updateHeadChange_TZero();
    }
}

void inline
Equation::updateHeadTZero() {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->updateHead_TZero();
        }
    }

void inline
Equation::updateTopZetasToHeads() {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->setTopZetaToHead();
        }
}

void inline
Equation::updateZetas_TZero() {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->updateZetasTZero();
        }
    }


void inline
Equation::adjustZetaHeights() {
    LOG(debug) << "Calculating vertical zeta movement";
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        //LOG(debug) << "nodeID: " << k << ", zeta 0: " << nodes->at(k)->getZeta(0).value() << ", zeta 1: " << nodes->at(k)->getZeta(1).value() << ", zeta 2: " << nodes->at(k)->getZeta(2).value();
        nodes->at(k)->zetaMovementBetweenLayers();
    }

    LOG(debug) << "Calculating horizontal zeta movement";
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        //LOG(debug) << "nodeID: " << k << ", zeta 0: " << nodes->at(k)->getZeta(0).value() << ", zeta 1: " << nodes->at(k)->getZeta(1).value() << ", zeta 2: " << nodes->at(k)->getZeta(2).value();
        nodes->at(k)->horizontalZetaMovement();
    }

    LOG(debug) << "Clipping inner zetas";
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        //LOG(debug) << "nodeID: " << k << ", zeta 0: " << nodes->at(k)->getZeta(0).value() << ", zeta 1: " << nodes->at(k)->getZeta(1).value() << ", zeta 2: " << nodes->at(k)->getZeta(2).value();
        nodes->at(k)->clipInnerZetas();
    }
    LOG(debug) << "Correcting crossing zetas";
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        //LOG(debug) << "nodeID: " << k << ", zeta 0: " << nodes->at(k)->getZeta(0).value() << ", zeta 1: " << nodes->at(k)->getZeta(1).value() << ", zeta 2: " << nodes->at(k)->getZeta(2).value();
        nodes->at(k)->correctCrossingZetas();
    }

    LOG(debug) << "Preventing zeta locking";
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        //LOG(debug) << "nodeID: " << k << ", zeta 0: " << nodes->at(k)->getZeta(0).value() << ", zeta 1: " << nodes->at(k)->getZeta(1).value() << ", zeta 2: " << nodes->at(k)->getZeta(2).value();
        nodes->at(k)->preventZetaLocking();
    }
}

void inline
Equation::updateZoneChange() {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->saveZoneChange();
    }
}

void inline
Equation::updateGNCBudget() {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->saveGNCMassBalance();
        }
    }

void inline
Equation::updateBudget() {
#pragma omp parallel for
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
    LOG(numerics) << "Initializing Matrix (before iteration)";
    updateEquation(); // updating matrix before iteration
    //LOG(debug) << "A.block<0,0>(10,10):\n" << A.block(0,0,100,100);
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxAllowedHeadChange, x); // among other: set x_t0 = x (initial)

    LOG(numerics) << "Running Time Step";

    double maxHeadChange{0};
    double oldMaxHeadChange{0};
    int itterScale{0};

    // Returns true if max headchange is greater than defined val
    auto isHeadChangeGreater = [this,&maxHeadChange]() -> bool {
        double changeMax{0.0};
        double headMax{0.0};
        double headSum{0.0};
        double headMean{0.0};
        double countAbsHeadAbove10k{0.0};
        double changeAtNode{0.0};
        double absHeadAtNode{0.0};
        large_num nodeID_headMax{0};
        large_num nodeID_changeMax{0};
#pragma omp parallel for
        for (large_num rowID = 0; rowID < numberOfActiveNodes; ++rowID) {
            auto nodeID = rowID_to_nodeID[rowID];
            if (rowID != -1) {
                changeAtNode = std::abs(
                        nodes->at(nodeID)->getProperties().get<quantity<Model::Meter>, Model::HeadChange>().value());
                nodeID_changeMax = (changeAtNode > changeMax) ? nodeID : nodeID_changeMax;
                changeMax = (changeAtNode > changeMax) ? changeAtNode : changeMax;
                absHeadAtNode = std::abs(nodes->at(nodeID)->getHead().value());
                nodeID_headMax = (absHeadAtNode > headMax) ? nodeID : nodeID_headMax;
                headMax = (absHeadAtNode > headMax) ? absHeadAtNode : headMax;
                headSum += nodes->at(nodeID)->getHead().value();
                if (absHeadAtNode > 10000) { countAbsHeadAbove10k++; }
            }
        }
        headMean = headSum / numberOfActiveNodes;
        LOG(numerics) << "MAX Absolute Head: " << headMax << " (nodeID: " << nodeID_headMax << "), "
                      << "Absolute head > 10,000 at " << countAbsHeadAbove10k << " nodes";
        LOG(numerics) << "MEAN Head: " << headMean;
        maxHeadChange = changeMax;
        LOG(numerics) << "MAX Absolute Head Change: " << changeMax << " (nodeID: " << nodeID_changeMax << ")";
        return changeMax > maxAllowedHeadChange;
    };

    long outerIterations{0};
    long innerIterations{0};
    bool headFail{false};
    char smallHeadChangeCounter{0};
    bool headConverged{false};
    while (outerIterations < MAX_OUTER_ITERATIONS) {
        LOG(numerics) << "Outer iteration: " << outerIterations;
        x = cg.solveWithGuess(b, x); // solving inner iterations
        //LOG(debug) << "x:\n" << x << std::endl;

        innerIterations = cg.iterations();
        LOG(numerics) << "Inner iterations: " << innerIterations;
        if (innerIterations == 0 and outerIterations == 0) {
            LOG(numerics) << "Convergence criterion met without iterations.";
            break;
        }
        updateIntermediateHeads();

        /**
         * @brief head change convergence
         */
        headFail = isHeadChangeGreater();
        if (headFail) { //convergence is not reached
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
        if(maxHeadChange == oldMaxHeadChange){
            //The head change is really the same -> increase inner iterations
            itterScale = itterScale + 10;
            cg.setMaxIterations(max_inner_iterations + itterScale);
        }
        oldMaxHeadChange = maxHeadChange;

        /**
         * @brief residual norm convergence
         */
        /*if (cg.info() == Success and iterations != 0) {
            LOG(numerics) << "cg solver success"; // this is always reached in outer iteration 1 (starts at 0)
            break;
        }*/
        LOG(numerics) << "Head change bigger: " << headFail;
        LOG(numerics) << "Updating Equation";
        updateEquation();
        outerIterations++;
    }
    //LOG(debug) << "matrix (A):\n" << A << std::endl;
    //LOG(debug) << "head (x):\n" << x << std::endl;
    //LOG(debug) << "rhs (b):\n" << b << std::endl;

    if (outerIterations == MAX_OUTER_ITERATIONS) {
        std::cerr << "Fail in solving matrix with max iterations\n";
        LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
        LOG(numerics) << "|Residual|_l2: " << cg.error();
    }

    __itter = outerIterations;
    __error = cg.error_inf();

    /**
     * ###############################
     * # Solve Zeta Surface Equation #
     * ###############################
     */
     if(vdf) {
         solve_zetas();
     }

    /**
    * ###############################
    * # Update budgets #
    * ###############################
    */
    if(gnc) {
        updateGNCBudget();
    }

    LOG(numerics) << "Updating head change and head of previous time step";
    updateHeadChangeTZero();
    updateHeadTZero();
    LOG(numerics) << "Updating budget";
    updateBudget();
}

/**
 * Solve Zeta Surface Equation
 */
void
Equation::solve_zetas(){
    LOG(numerics) << "Solving for zeta surfaces";
    LOG(numerics) << "If unconfined: clipping top zeta to new surface heights";
    updateTopZetasToHeads();
    for (large_num layer = 0; layer < numberOfLayers; layer++) {
        LOG(numerics) << "Finding zeta surface heights in layer " << layer;
        for (int localZetaID = 1; localZetaID < numberOfZones; localZetaID++) {
            LOG(numerics) << "Solving for zeta surface " << localZetaID << "";
            LOG(numerics) << "Updating Matrix (zetas)";
            updateEquation_zetas(layer, localZetaID);
            if (A_zetas.size() == 0){ // if matrix empty continue with next iteration
                continue;
            }

            double maxZetaChange{0};
            double oldMaxZetaChange{0};
            int itterScale{0};

            // Returns true if max zetachange is greater than defined val
            auto isZetaChangeGreater = [this, &maxZetaChange, &layer, &localZetaID]() -> bool {
                double absZetaChangeMax{0.0};
                double absZetaChange{0.0};
                double absZeta{0.0};
                double absZetaMax{0.0};
                large_num offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for
                for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
                    large_num nodeID = k + offset;
                    auto rowID = iter_to_rowID[k];
                    if (rowID != -1) {
                        absZetaChange = std::abs(nodes->at(nodeID)->getZetaChange(localZetaID).value());
                        absZetaChangeMax = (absZetaChange > absZetaChangeMax) ? absZetaChange : absZetaChangeMax;
                        absZeta = std::abs(nodes->at(nodeID)->getZeta(localZetaID).value());
                        absZetaMax = (absZeta > absZetaMax) ? absZeta : absZetaMax;
                    }
                }
                LOG(numerics) << "MAX Absolute Zeta: " << absZetaMax;
                maxZetaChange = absZetaChangeMax;
                LOG(numerics) << "MAX Absolute Zeta Change: " << maxZetaChange;
                return maxZetaChange > maxAllowedZetaChange;
            };

            long outerIteration{0};
            long innerIteration{0};
            bool zetaFail{false};
            char smallZetaChanges{0};
            bool zetaConverged{false};

            //LOG(debug) << "A_zetas (before iteration):\n" << A_zetas << std::endl;
            //LOG(debug) << "b_zetas (before iteration):\n" << b_zetas << std::endl;
            //LOG(debug) << "x_zetas (before iteration):\n" << x_zetas << std::endl;

            while (outerIteration < MAX_OUTER_ITERATIONS) {
                LOG(numerics) << "Outer iteration (zetas): " << outerIteration;
                //Solve inner iterations
                x_zetas = cg_zetas.solveWithGuess(b_zetas, x_zetas);
                //LOG(debug) << "x_zetas:\n" << x_zetas << std::endl;
                updateIntermediateZetas(layer, localZetaID);

                innerIteration = cg_zetas.iterations();
                if (innerIteration == 0 and outerIteration == 0) {
                    LOG(numerics) << "Convergence criterion for zeta (=" << RCLOSE_ZETA << ") too small - no iterations";
                    break;
                }

                /**
                 * @brief zeta change convergence // todo make function of this
                 */
                zetaFail = isZetaChangeGreater();
                if (zetaFail) {
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

                if (maxZetaChange == oldMaxZetaChange) {
                    //The zeta change is really the same -> increase inner iterations
                    itterScale = itterScale + 10;
                    cg_zetas.setMaxIterations(max_inner_iterations + itterScale);
                }
                oldMaxZetaChange = maxZetaChange;

                /**
                 * @brief residual norm convergence // Question: make function of this?
                 */
                /*LOG(numerics) << "Inner iterations (zetas): " << innerItter;
                if (cg_zetas.info() == Success and iterations != 0) {
                    LOG(numerics) << "cg_zetas success";
                    break;
                }*/
                LOG(numerics) << "Zeta change bigger: " << zetaFail;

                LOG(numerics) << "Updating Equation";
                updateEquation_zetas(layer, localZetaID);
                outerIteration++;
            }

            if (outerIteration == MAX_OUTER_ITERATIONS) {
                LOG(userinfo) << "Fail in solving zeta matrix with max iterations"
                              << "(localZetaID = " << localZetaID << ", layer = " << layer << ")";
                LOG(numerics) << "|Residual|_inf / |RHS|_inf (zetas): " << cg_zetas.error_inf();
                LOG(numerics) << "|Residual|_l2 (zetas): " << cg_zetas.error();
            }

            //LOG(debug) << "x_zetas[" << localZetaID << "] on layer " << layer << ":\n" << x_zetas << std::endl;

        }
    }
    updateZoneChange();

    LOG(numerics) << "Adjusting zeta heights (after zeta height convergence)";
    adjustZetaHeights();
    LOG(numerics) << "Updating zetasTZero";
    updateZetas_TZero();
}


int
Equation::getItter() {
    return __itter;
}

double
Equation::getError() {
    return __error;
}

}
}//ns
