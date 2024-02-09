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
    this->vdf = options.isDensityVariable(); // vdf = Variable Density Flow
    this->gnc = options.isGridRefined(); // gnc = Ghost Node Correction

    this->max_inner_iterations = options.getInnerItter();
    this->maxRefinement = options.getMaxRefinement();
    this->nodes = nodes;

    //set inner iterations
    cg.setMaxIterations(max_inner_iterations);
    cg.setTolerance(RCLOSE_HEAD);

    long_vector __b(numberOfNodesTotal);
    b = std::move(__b);
    long_vector __x(numberOfNodesTotal);
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
Equation::addToA(large_num &nodeID){
    for (const auto &entry: nodes->at(nodeID)->getMatrixEntries()) {
        //tripletList.push_back(T(nodeID,entry.first,entry.second.value()));
    }
}

void inline
Equation::addToA_zeta(large_num nodeID, large_num offset, int localZetaID) {
    for (const auto &[nodeID_neig, zoneConductance] : nodes->at(nodeID)->getMatrixEntries(localZetaID)) {
        if (nodeID_to_rowID[nodeID_neig] != -1) {
            NANChecker(zoneConductance.value(), "Matrix entry (zetas)");
            //LOG(debug) << "colID = " << colID << ", rowID = " << rowID << ", iterOffset = " << iterOffset;
            A_zetas.coeffRef(nodeID_to_rowID[nodeID], nodeID_to_rowID[nodeID_neig]) = zoneConductance.value();
        }
    }
}

void inline
Equation::updateEquation() {
    Eigen::SparseMatrix<pr_t> sparseMatrix(numberOfNodesTotal, numberOfNodesTotal);
    A = std::move(sparseMatrix);
    int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 2 + 1; // +2 for top/down, + 1 for this node
    A.reserve(long_vector::Constant(numberOfNodesTotal, numberOfEntries));
//#ifdef EIGEN_HAS_OPENMP
    //Eigen::initParallel();
    //Index threads = Eigen::nbThreads();
//#endif
//#pragma omp parallel for //schedule(dynamic,(numberOfNodesTotal+threads*4-1)/(threads*4)) num_threads(threads)
    for (large_num nodeID = 0; nodeID < numberOfNodesTotal; ++nodeID) {
        //---------------------Left
        const auto matrixEntries = nodes->at(nodeID)->getMatrixEntries();
        for (const auto &[nodeID_neig, conductance] : matrixEntries) {
            //NANChecker(conductance.value(), "NAN in conductance");
            A.coeffRef(nodeID, nodeID_neig) = conductance.value();
        }
    }

    //if (!A.isApprox(A.adjoint())) {
    //    LOG(userinfo) << "Matrix NOT self-adjoint. Preconditioner will fail.";
    //}

    if (!A.isCompressed()) {
        A.makeCompressed();
    }

    LOG(numerics) << "Preconditioning matrix before solver";
    cg.compute(A);
    if (cg.info() != Success) {
        LOG(numerics) << "Fail in preconditioning matrix";
        throw "Fail in preconditioning matrix";
    }
    //preconditionMatrix();

#pragma omp parallel for //schedule(dynamic,(numberOfNodesTotal+threads*4-1)/(threads*4)) num_threads(threads)
    for (large_num nodeID = 0; nodeID < numberOfNodesTotal; ++nodeID) {
        x(nodeID) = nodes->at(nodeID)->getHead().value();
        //---------------------Right
        b(nodeID) = nodes->at(nodeID)->getRHS().value();
        NANChecker(b[nodeID], "Right hand side");
    }
    LOG(numerics) << "Equation is ready for solver";
}

void inline
Equation::updateEquation_zetas(large_num layer, int localZetaID) {
    large_num numInactive{0};
    nodeID_to_rowID.clear();

    large_num offset = layer * numberOfNodesPerLayer;
    // finding inactive nodes
//#pragma omp parallel for
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer+offset; nodeID++) {
        if ( nodes->at(nodeID)->isZetaActive(localZetaID) ) {
            nodeID_to_rowID[nodeID] = nodeID - offset - numInactive;
        } else {
            ++numInactive; // tracking how many have been set inactive
            nodeID_to_rowID[nodeID] = -1; // these entries will be ignored ( e.g. in loop filling A_zeta, x_zeta and b_zeta)
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
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer+offset; nodeID++) {
        if (nodeID_to_rowID[nodeID] != -1) {
            //---------------------Left: fill matrix A_zeta and initiate x_zetas
            //addToA_zeta(nodeID, offset, localZetaID);
            for (const auto &[nodeID_neig, zoneConductance] : nodes->at(nodeID)->getMatrixEntries(localZetaID)) {
                if (nodeID_to_rowID[nodeID_neig] != -1) {
                    NANChecker(zoneConductance.value(), "Matrix entry (zetas)");
                    //LOG(debug) << "rowID = " << nodeID_to_rowID[nodeID] << ", colID = " << nodeID_to_rowID[nodeID_neig];
                    A_zetas.coeffRef(nodeID_to_rowID[nodeID], nodeID_to_rowID[nodeID_neig]) = zoneConductance.value();
                }
            }
            x_zetas(nodeID_to_rowID[nodeID]) = nodes->at(nodeID)->getZeta(localZetaID).value();
            //---------------------Right
            b_zetas(nodeID_to_rowID[nodeID]) = nodes->at(nodeID)->getRHS(localZetaID).value();
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
        LOG(userinfo) << "Matrix NOT self-adjoint (zetas). Preconditioner will fail.";
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
Equation::updateHeads() {
    //LOG(debug) << "Updating heads";
#pragma omp parallel for
    for (long rowID = 0; rowID < numberOfNodesTotal; rowID++) {
        nodes->at(rowID)->setHead(x[rowID] * si::meter);
    }
}

void inline
Equation::updateHeadChange() {
    //LOG(debug) << "Updating head change";
    long_vector changes = adaptiveDamping.getChanges(getResiduals(), x, isAdaptiveDamping);
        //LOG(debug) << "changes:\n" << changes;
#pragma omp parallel for
        for (long rowID = 0; rowID < numberOfNodesTotal; rowID++) {
            nodes->at(rowID)->setHeadChange(changes[rowID] * si::meter);
        }
    }

void inline
Equation::updateZetas(large_num layer, int localZetaID) {
    auto offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer+offset; nodeID++) {
        if (nodeID_to_rowID[nodeID] != -1) {
            nodes->at(nodeID)->setZeta(localZetaID, x_zetas[nodeID_to_rowID[nodeID]] * si::meter);
        }
    }
}

void inline
Equation::updateZetaChange(large_num layer, int localZetaID) {
    auto offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer+offset; nodeID++) {
        if (nodeID_to_rowID[nodeID] != -1) {
            nodes->at(nodeID)->setZetaChange(localZetaID, x_zetas[nodeID_to_rowID[nodeID]] * si::meter);
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
    LOG(numerics) << "Initializing Equation (before iteration)";
    updateEquation();
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxAllowedHeadChange, x); // among other: set x_t0 = x (initial)

    double oldMaxHeadChange{0};
    int itterScale{0};
    long outerIterations{0};
    long innerIterations{0};
    bool headFail{false};
    char smallHeadChangeCounter{0};
    bool headConverged{false};
    while (outerIterations < MAX_OUTER_ITERATIONS) {
        LOG(userinfo) << "Outer iteration: " << outerIterations;
        //LOG(debug) << "x.block before:\n" << x.block(0,0,5,1); // numberOfNodesPerLayer
        x = cg.solveWithGuess(b, x); // solving inner iterations
        //LOG(debug) << "A.block:\n" << A.block(0,0,5,5); // startRow, startCol, numRows, numCol
        //LOG(debug) << "b.block:\n" << b.block(0,0,5,1);
        //LOG(debug) << "x.block (layer 0):\n" << x.block(0,0,5,1);
        //LOG(debug) << "x.block (layer 1):\n" << x.block(numberOfNodesPerLayer,0,5,1);

        innerIterations = cg.iterations();
        //LOG(numerics) << "Inner iterations: " << innerIterations;
        if (innerIterations == 0 and outerIterations == 0) {
            LOG(numerics) << "Convergence criterion met without iterations.";
            break;
        }
        updateHeads(); // needs to be before "headFail = isHeadChangeGreater();"
        updateHeadChange(); // needs to be before "headFail = isHeadChangeGreater();"

        /**
         * @brief head change convergence
         */
        headFail = isHeadChangeGreater();
        LOG(numerics) << "Head change bigger: " << headFail;

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
        if(maxCurrentHeadChange == oldMaxHeadChange){
            //The head change is really the same -> increase inner iterations
            itterScale = itterScale + 10;
            cg.setMaxIterations(max_inner_iterations + itterScale);
        }
        oldMaxHeadChange = maxCurrentHeadChange;

        /**
         * @brief residual norm convergence
         */
        /*if (cg.info() == Success and iterations != 0) {
            LOG(numerics) << "cg solver success"; // this is always reached in outer iteration 1 (starts at 0)
            break;
        }*/
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

bool inline
Equation::isHeadChangeGreater(){
    maxCurrentHeadChange = 0;
    double headMax{0.0};
    double headSum{0.0};
    double headMean{0.0};
    double countAbsHeadAbove10k{0.0};
    double changeAtNode{0.0};
    double absHeadAtNode{0.0};
    large_num nodeID_headMax{0};
    large_num nodeID_changeMax{0};
#pragma omp parallel for
    for (large_num nodeID = 0; nodeID < numberOfNodesTotal; nodeID++) {
        changeAtNode = std::abs(
                nodes->at(nodeID)->getProperties().get<quantity<Model::Meter>, Model::HeadChange>().value());
        nodeID_changeMax = (changeAtNode > maxCurrentHeadChange) ? nodeID : nodeID_changeMax;
        maxCurrentHeadChange = (changeAtNode > maxCurrentHeadChange) ? changeAtNode : maxCurrentHeadChange;
        absHeadAtNode = std::abs(nodes->at(nodeID)->getHead().value());
        nodeID_headMax = (absHeadAtNode > headMax) ? nodeID : nodeID_headMax;
        headMax = (absHeadAtNode > headMax) ? absHeadAtNode : headMax;
        headSum += nodes->at(nodeID)->getHead().value();
        if (absHeadAtNode > 10000) {
            //LOG(debug) << "nodeID: " << nodeID << "-> head = " << nodes->at(nodeID)->getHead().value() << ", x = " << x[nodeID];
            countAbsHeadAbove10k++;}
    }
    headMean = headSum / numberOfNodesTotal;
    LOG(userinfo) << "MAX Absolute Head: " << headMax << " (nodeID: " << nodeID_headMax << "), "
                  << "Absolute head > 10,000 at " << countAbsHeadAbove10k << " nodes";
    LOG(userinfo) << "MEAN Head: " << headMean;
    LOG(userinfo) << "MAX Absolute Head Change: " << maxCurrentHeadChange << " (nodeID: " << nodeID_changeMax << ")";
    return maxCurrentHeadChange > maxAllowedHeadChange;
}


bool inline
Equation::isZetaChangeGreater(large_num layer, int localZetaID){
    maxCurrentZetaChange = 0;
    double absZetaChange{0.0};
    double absZeta{0.0};
    double absZetaMax{0.0};
    large_num offset = layer * numberOfNodesPerLayer;
#pragma omp parallel for
    for (large_num nodeID = offset; nodeID < numberOfNodesPerLayer+offset; nodeID++) {
        if (nodeID_to_rowID[nodeID] != -1) {
            absZetaChange = std::abs(nodes->at(nodeID)->getZetaChange(localZetaID).value());
            maxCurrentZetaChange = (absZetaChange > maxCurrentZetaChange) ? absZetaChange : maxCurrentZetaChange;
            absZeta = std::abs(nodes->at(nodeID)->getZeta(localZetaID).value());
            absZetaMax = (absZeta > absZetaMax) ? absZeta : absZetaMax;
        }
    }
    LOG(numerics) << "MAX Absolute Zeta: " << absZetaMax;
    LOG(numerics) << "MAX Absolute Zeta Change: " << maxCurrentZetaChange;
    return maxCurrentZetaChange > maxAllowedZetaChange;
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

            long outerIteration{0};
            long innerIteration{0};
            bool zetaFail{false};
            char smallZetaChanges{0};
            bool zetaConverged{false};

            //LOG(debug) << "A_zetas (before iteration):\n" << A_zetas << std::endl;
            //LOG(debug) << "b_zetas (before iteration):\n" << b_zetas << std::endl;
            //LOG(debug) << "x_zetas (before iteration):\n" << x_zetas << std::endl;

            while (outerIteration < MAX_OUTER_ITERATIONS) {
                LOG(userinfo) << "Outer iteration (zetas): " << outerIteration;
                //Solve inner iterations
                x_zetas = cg_zetas.solveWithGuess(b_zetas, x_zetas);
                //LOG(debug) << "x_zetas:\n" << x_zetas << std::endl;
                updateZetaChange(layer, localZetaID);// before updateZetas
                updateZetas(layer, localZetaID);
                innerIteration = cg_zetas.iterations();
                if (innerIteration == 0 and outerIteration == 0) {
                    LOG(numerics) << "Convergence criterion for zeta (=" << RCLOSE_ZETA << ") too small - no iterations";
                    break;
                }

                /**
                 * @brief zeta change convergence // todo make function of this
                 */
                zetaFail = isZetaChangeGreater(layer, localZetaID);
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
