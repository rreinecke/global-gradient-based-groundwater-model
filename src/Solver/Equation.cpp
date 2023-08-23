#include "Equation.hpp"

namespace GlobalFlow {
namespace Solver {

Equation::Equation(NodeVector nodes, Simulation::Options options) : options(options) {
    this->numberOfNodesPerLayer = options.getNumberOfNodesPerLayer();
    this->numberOfLayers = options.getNumberOfLayers();
    this->numberOfNodesTotal = numberOfNodesPerLayer * numberOfLayers;
    LOG(userinfo) << "Setting up Equation for " << numberOfNodesPerLayer << " nodes on " << numberOfLayers <<
    " layer(s) (in total " << numberOfNodesTotal << " nodes)" << std::endl;

    this->IITER = options.getMaxIterations();
    this->RCLOSE_HEAD = options.getConverganceCriteriaHead();
    this->RCLOSE_ZETA = options.getConverganceCriteriaZeta();
    this->initialHead = options.getInitialHead();
    this->maxHeadChange = options.getMaxHeadChange();
    this->isAdaptiveDamping = options.isDampingEnabled();
    this->dampMin = options.getMinDamp();
    this->dampMax = options.getMaxDamp();
    this->vdf = options.isDensityVariable();

    this->inner_iterations = options.getInnerItter();

    this->nodes = nodes;


    long_vector __x(numberOfNodesTotal);
    long_vector __b(numberOfNodesTotal);

    x = std::move(__x);
    b = std::move(__b);

    Eigen::SparseMatrix<pr_t> __A(numberOfNodesTotal, numberOfNodesTotal);

    A = std::move(__A);
    A.reserve(long_vector::Constant(numberOfNodesTotal, 7));

    //Init first result vector x by writing initial heads
    //Initial head should be positive
    //resulting head is the real hydraulic head
    double tmp = 0;
#pragma omp parallel for
    for (int i = 0; i < numberOfNodesTotal; ++i) {
        if (not simpleHead) {
            tmp = nodes->at(i)->calcInitialHead(initialHead * si::meter).value();
            nodes->at(i)->setHead_direct(tmp);
            NANChecker(tmp, "Initial Head");
            x[nodes->at(i)->getProperties().get<large_num, Model::ID>()] = std::move(tmp);
        } else {
            x[nodes->at(i)->getProperties().get<large_num, Model::ID>()] =
                    nodes->at(i)->getProperties().get<quantity<Model::Meter>, Model::Head>().value();
            nodes->at(i)->initHead_t0();
        }
    }

    //set inner iterations
    cg.setMaxIterations(inner_iterations);
    cg.setTolerance(RCLOSE_HEAD);
    //cg.preconditioner().setInitialShift(1e-8);


    if(vdf) {
        LOG(userinfo) << "Simulating variable density flow" << std::endl;
        this->numberOfZones = options.getDensityZones().size();
        this->maxZetaChange = options.getMaxZetaChange();

#pragma omp parallel for
        for (int i = 0; i < numberOfNodesTotal; ++i) {
            nodes->at(i)->initZetasTZero();
        }

        //set inner iterations
        cg_zetas.setMaxIterations(inner_iterations);
        cg_zetas.setTolerance(RCLOSE_ZETA);
    }
}

Equation::~Equation() {
    LOG(debug) << "Destroying equation\n" << std::endl;
}

void inline
Equation::addToA(std::unique_ptr<Model::NodeInterface> const &node) {

    std::unordered_map<large_num, quantity<Model::MeterSquaredPerTime>> map;
    map = node->getMatrixEntries();

    auto nodeID = node->getProperties().get<large_num, Model::ID>();

    for (const auto &conductance : map) {
        // conductance.first is the nodeID of the respective neighbour node
        NANChecker(conductance.second.value(), "Matrix entry");
        if (isCached) {
            A.coeffRef(nodeID, conductance.first) = conductance.second.value();
        } else {
            A.insert(nodeID, conductance.first) = conductance.second.value();
        }
    }
}

void inline
Equation::addToA_zeta(large_num nodeIter, large_num iterOffset, int localZetaID) {
    std::unordered_map<large_num, quantity<Model::MeterSquaredPerTime>> map;
    large_num rowID = index_mapping[nodeIter];
    large_num colID;
    quantity<Model::MeterSquaredPerTime> zoneConductance;
    map = nodes->at(nodeIter + iterOffset)->getMatrixEntries(localZetaID); // gets matrix entries (zone conductances and porosity term)

    for (const auto &entry : map) { // entry: [1] node id of horizontal neighbours or this node, [2] conductance to neighbours in zeta zone
        colID = index_mapping[entry.first - iterOffset];
        if (colID != -1) {
            NANChecker(entry.second.value(), "Matrix entry (zetas)");
            zoneConductance = entry.second;
            //LOG(debug) << "colID = " << colID << ", rowID = " << rowID << ", iterOffset = " << iterOffset;
            //if (isCached_zetas) {
                A_zetas.coeffRef(rowID, colID) = zoneConductance.value();
            //} else {
            //    A_zetas.insert(rowID, colID) = zoneConductance.value();
            //}
        }
    }
}

void inline
Equation::updateMatrix() {
    Index n = A.outerSize();

#ifdef EIGEN_HAS_OPENMP
    Eigen::initParallel();
    Index threads = Eigen::nbThreads();
#endif
#pragma omp parallel for schedule(dynamic,(n+threads*4-1)/(threads*4)) num_threads(threads)
    for (large_num j = 0; j < numberOfNodesTotal; ++j) {
        //LOG(userinfo) << "node: " << j;
        large_num id = nodes->at(j)->getProperties().get<large_num, Model::ID>();
        //---------------------Left
        addToA(nodes->at(id));
        //---------------------Right
        b(id) = nodes->at(id)->getRHS().value();
        NANChecker(b[id], "Right hand side");
    }

    //Check if after iteration former 0 values turned to non-zero
    if ((not A.isCompressed()) and isCached) { // Question: is this necessary?
        LOG(numerics) << "Recompressing Matrix";
        A.makeCompressed();
    }

    //LOG(debug) << "A (updateMatrix):\n" << A << std::endl;
    //LOG(debug) << "x (updateMatrix):\n" << x << std::endl;
    //LOG(debug) << "b (updateMatrix):\n" << b << std::endl;
}

void inline
Equation::updateMatrix_zetas(large_num iterOffset, int localZetaID) {
    Index n = A_zetas.outerSize();
    large_num numInactive{0};
    index_mapping.clear();

    // finding inactive nodes
//#pragma omp parallel for
    for (large_num i = 0; i < numberOfNodesPerLayer; ++i) {
        if ( nodes->at(i + iterOffset)->isZetaActive(localZetaID) ) {
            index_mapping[i] = i - numInactive;
        } else {
            ++numInactive; // tracking how many have been set inactive
            index_mapping[i] = -1; // these entries will be ignored ( e.g. in loop filling A_zeta, x_zeta and b_zeta)
        }
    }

    const long numActive = numberOfNodesPerLayer - numInactive;
    Eigen::SparseMatrix<pr_t> __A_zetas(numActive, numActive);
    A_zetas = std::move(__A_zetas);
    A_zetas.reserve(long_vector::Constant(numActive, 5));
    long_vector __b_zetas(numActive);
    b_zetas = std::move(__b_zetas);
    long_vector __x_zetas(numActive);
    x_zetas = std::move(__x_zetas);
    
//#pragma omp parallel for
    for (large_num j = 0; j < numberOfNodesPerLayer; ++j) {
        auto id = index_mapping[j];
        if (id != -1) {
            //---------------------Left: fill matrix A_zeta and initiate x_zetas
            addToA_zeta(j, iterOffset, localZetaID);
            x_zetas(id) = nodes->at(j + iterOffset)->getZeta(localZetaID).value();
            //---------------------Right
            b_zetas(id) = nodes->at(j + iterOffset)->getRHS(localZetaID).value();
        }
    }

    //Check if after iteration former 0 values turned to non-zero
    if ((not A_zetas.isCompressed()) and isCached_zetas) {
        LOG(numerics) << "Recompressing Matrix (zetas)";
        A_zetas.makeCompressed();
    }
}

void inline
Equation::preconditioner() {
    LOG(numerics) << "Decomposing Matrix";
    cg.compute(A);
    if (cg.info() != Success) {
        LOG(numerics) << "Fail in decomposing matrix";
        throw "Fail in decomposing matrix";
    }
}

void inline
Equation::preconditioner_zetas() {
    Eigen::SelfAdjointEigenSolver<SparseMatrix<pr_t>> selfAdjointSolver(A_zetas);
    if ( selfAdjointSolver.info() != Success) {
        LOG(numerics) << "A_zetas is not self-adjoint!";
        throw "Fail before decomposing";
    }

    LOG(numerics) << "Decomposing Matrix (zetas)";
    cg_zetas.compute(A_zetas);
    if (cg_zetas.info() != Success) {
        LOG(userinfo) << "Fail in decomposing matrix (zetas)";
        throw "Fail in decomposing matrix (zetas)";
    }
}

void inline
Equation::updateIntermediateHeads() {
    long_vector changes = adaptiveDamping.getDamping(getResiduals(), x, isAdaptiveDamping);

#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        large_num id = nodes->at(k)->getProperties().get<large_num, Model::ID>();
        //nodes->at(k)->setHead((double) heads[id] * si::meter);
        nodes->at(id)->setHeadChange((double) changes[id] * si::meter);
        //LOG(debug) << "head (updateIntermediateHeads): " << nodes->at(id)->getHead().value();
    }
}

void inline
Equation::updateIntermediateZetas(large_num iterOffset, int localZetaID) {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
        auto id = index_mapping[k];
        if (id != -1) {
            // todo get nodeID
            nodes->at(k + iterOffset)->setZeta(localZetaID, (double) x_zetas[id] * si::meter);
            nodes->at(k + iterOffset)->setZetaChange(localZetaID, (double) x_zetas[id] * si::meter);
        }
    }
}

void inline
Equation::updateFinalHeads() {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesTotal; ++k) {
        nodes->at(k)->updateHeadChange();
    }
}

void inline
Equation::updateTopZetasToHeads() {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->setTopZetaToHead();
        }
}

/*void inline // todo: in SWI2: required for zeta time step adjustment
Equation::checkAllZetaSlopes(int localZetaID) {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
            nodes->at(k)->checkZetaSlopes();
        }
}*/

void inline
Equation::updateZetasAfterEquation() {
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
            nodes->at(k)->verticalZetaMovement();
            //if (k >= 264 && k <= 266){
            //    LOG(debug) << "zeta[" << 1 << "] at node " << k << " after verticalZetaMovement: " << nodes->at(k)->getZeta(1).value();
            //}
        }

    LOG(debug) << "Calculating horizontal zeta movement";
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->horizontalZetaMovement();
        }

        LOG(debug) << "Clipping inner zetas";
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->clipInnerZetas();
        }
        LOG(debug) << "Correcting crossing zetas";
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->correctCrossingZetas();
        }

        LOG(debug) << "Preventing zeta locking";
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            nodes->at(k)->preventZetaLocking();
        }

/*#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            for (int localZetaID = 1; localZetaID < numberOfZones; localZetaID++) {
                LOG(debug) << "zeta surface " << localZetaID << " at node " << k << ": " << nodes->at(k)->getZeta(localZetaID).value();
            }
        }*/
    }

void inline
Equation::updateVDFBudget(large_num iterOffset) {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
        nodes->at(k + iterOffset)->saveVDFMassBalance();
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
    LOG(numerics) << "Updating Matrix";
    updateMatrix();

    if (!isCached) {
        LOG(numerics) << "Compressing matrix";
        A.makeCompressed();

        LOG(numerics) << "Cached Matrix";
        isCached = true;
    }

    preconditioner(); // decomposing matrix
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxHeadChange, x);

    LOG(numerics) << "Running Time Step";

    double maxHead{0};
    double oldMaxHead{0};
    int itterScale{0};

    // Returns true if max headchange is greater than defined val
    auto isHeadChangeGreater = [this,&maxHead]() -> bool {
        double lowerBound = maxHeadChange;
        double changeMax = 0;
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            double val = std::abs(
                    nodes->at(k)->getProperties().get<quantity<Model::Meter>, Model::HeadChange>().value());
            changeMax = (val > changeMax) ? val : changeMax;
        }
	    maxHead = changeMax;
        LOG(numerics) << "MAX Head Change: " << changeMax;
        return changeMax > lowerBound;
    };

    long iterations{0};
    bool headFail{false};
    char smallHeadChanges{0};
    bool headConverged{false};

    while (iterations < IITER) {
        LOG(numerics) << "Outer iteration: " << iterations;
        //LOG(debug) << "x (before outer iteration " << iterations << "):\n" << x << std::endl;

        //Solve inner iterations
        x = cg.solveWithGuess(b, x);
        //LOG(debug) << "b (Outer iteration " << iterations << "):\n" << b << std::endl;
        //LOG(debug) << "x (Outer iteration " << iterations << "):\n" << x << std::endl;

        updateIntermediateHeads();

        int innerItter{0};
        innerItter = cg.iterations();
        if (innerItter == 0 and iterations == 0) {
            LOG(numerics) << "convergence criterion to small - no iterations";
            break;
        }

        /**
         * @brief head change convergence
         */
        headFail = isHeadChangeGreater();
        if (headFail) {
            //convergence is not reached
            //reset counter
            smallHeadChanges = 0;
        } else {
            //itter converged with head criterion
            if (headConverged and smallHeadChanges == 0) {
                LOG(numerics) << "Conditional convergence - check mass balance";
                break;
            } else {
                headConverged = true;

                smallHeadChanges++;
                if (smallHeadChanges >= 2) {
                    LOG(numerics) << "Reached head change convergence";
                    //LOG(debug) << "x (converged):\n" << x << std::endl;
                    break;
                }
            }
        }

        if(maxHead == oldMaxHead){
            //The head change is really the same -> increase inner iterations
            itterScale = itterScale + 10;
        }else{
            itterScale = 0;
        }
        cg.setMaxIterations(inner_iterations + itterScale);
        oldMaxHead = maxHead;

        /**
         * @brief residual norm convergence
         */
        LOG(numerics) << "Inner iterations: " << innerItter;
        if (cg.info() == Success and iterations != 0) {
            LOG(numerics) << "cg solver success"; // this is always reached in outer iteration 1 (starts at 0)

            break;
        }

        LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
        LOG(numerics) << "|Residual|_l2: " << cg.error();
        LOG(numerics) << "Head change bigger: " << headFail;

        updateMatrix();
        preconditioner();

        iterations++;
    }
    //LOG(debug) << "A:\n" << A << std::endl;
    //LOG(debug) << "x:\n" << x << std::endl;
    //LOG(debug) << "b (= rhs):\n" << b << std::endl;

    if (iterations == IITER) {
        std::cerr << "Fail in solving matrix with max iterations\n";
        LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
        LOG(numerics) << "|Residual|_l2: " << cg.error();
    }

    updateFinalHeads();
    updateBudget();

    __itter = iterations;
    __error = cg.error_inf();


    /**
     * ###############################
     * # Solve Zeta Surface Equation #
     * ###############################
     */
     if(vdf) {
         solve_zetas();
     }
}

/**
 * Solve Zeta Surface Equation
 */
void
Equation::solve_zetas(){
    LOG(numerics) << "If unconfined: clipping top zeta to new surface heights";
    updateTopZetasToHeads();

    for (large_num layer = 0; layer < numberOfLayers; layer++) {
        large_num iterOffset = layer * numberOfNodesPerLayer;
        LOG(debug) << "Finding zeta surface heights in layer " << layer;
#pragma omp parallel for
        for (int localZetaID = 1; localZetaID < numberOfZones; localZetaID++) {
            LOG(debug) << "Solving for zeta surface " << localZetaID << "";
            LOG(numerics) << "Updating Matrix (zetas)";
            updateMatrix_zetas(iterOffset, localZetaID);

            if (A_zetas.size() == 0){ // if matrix empty continue with next iteration
                continue;
            }

            if (!isCached_zetas) {
                LOG(numerics) << "Compressing matrix (zetas)";
                A_zetas.makeCompressed();

                LOG(numerics) << "Cached Matrix (zetas)";
                isCached_zetas = true;
            }

            //LOG(debug) << "A_zetas (before preconditioner):\n" << A_zetas << std::endl;
            preconditioner_zetas();

            double maxZeta{0};
            double oldMaxZeta{0};
            int itterScale{0};

            // Returns true if max zetachange is greater than defined val
            auto isZetaChangeGreater = [this, &maxZeta, &iterOffset]() -> bool {
                double lowerBound = maxZetaChange;
                double changeMax = 0;

#pragma omp parallel for
                for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
                    double val;
                    for (int l = 1; l < numberOfZones; l++) { // localZetaID needs to be defined within "isZetaChangeGreater"
                        val = std::abs(nodes->at(k + iterOffset)->getZetaChange(l).value()); // todo improve for loop (by getting rid of it)
                        //LOG(debug) << "val (zetas change) (in solve_zeta): " << val << std::endl;
                        changeMax = (val > changeMax) ? val : changeMax;
                    }
                }
                maxZeta = changeMax;
                LOG(numerics) << "MAX Zeta Change: " << changeMax;
                return changeMax > lowerBound;
            };

            long iterations{0};
            bool zetaFail{false};
            char smallZetaChanges{0};
            bool zetaConverged{false};

            //LOG(debug) << "A_zetas (before iteration):\n" << A_zetas << std::endl;
            //LOG(debug) << "b_zetas (before iteration):\n" << b_zetas << std::endl;
            //LOG(debug) << "x_zetas (before iteration):\n" << x_zetas << std::endl;

            while (iterations < IITER) {
                LOG(numerics) << "Outer iteration (zetas): " << iterations;

                //Solve inner iterations
                x_zetas = cg_zetas.solveWithGuess(b_zetas, x_zetas);
                //LOG(debug) << "x_zetas of layer " << layer << " (after outer iteration " << iterations << "):\n" << x_zetas << std::endl;

                updateIntermediateZetas(iterOffset, localZetaID);

                int innerItter = cg_zetas.iterations();

                if (innerItter == 0 and iterations == 0) {
                    LOG(numerics) << "Zeta surfaces: convergence criterion to small - no iterations";
                    break;
                }

                /**
                 * @brief zeta change convergence // todo make function of this
                 */
                zetaFail = isZetaChangeGreater();
                if (zetaFail) {
                    //convergence is not reached
                    //reset counter
                    smallZetaChanges = 0;
                } else {
                    //itter converged with head criterion
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

                if (maxZeta == oldMaxZeta) {
                    //The zeta change is really the same -> increase inner iterations
                    itterScale = itterScale + 10;
                } else {
                    itterScale = 0;
                }
                cg_zetas.setMaxIterations(inner_iterations + itterScale);
                oldMaxZeta = maxZeta;

                /**
                 * @brief residual norm convergence // todo make function of this
                 */

                LOG(numerics) << "Inner iterations (zetas): " << innerItter;
                if (cg_zetas.info() == Success and iterations != 0) {
                    LOG(numerics) << "cg_zetas success";
                    break;
                }

                LOG(numerics) << "|Residual|_inf / |RHS|_inf (zetas): " << cg_zetas.error_inf();
                LOG(numerics) << "|Residual|_l2 (zetas): " << cg_zetas.error();
                LOG(numerics) << "Zeta change bigger: " << zetaFail;


                updateMatrix_zetas(layer, localZetaID);
                if (A_zetas.size() == 0) { // if matrix is empty, go to next iteration
                    continue;
                }
                //LOG(debug) << "A_zetas (after outer iteration " << iterations << "):\n" << A_zetas << std::endl;
                //LOG(debug) << "b_zetas (after outer iteration " << iterations << "):\n" << b_zetas << std::endl;
                preconditioner_zetas();

                iterations++;
            }

            if (iterations == IITER) {
                std::cerr << "Fail in solving matrix with max iterations (zetas)\n";
                LOG(numerics) << "|Residual|_inf / |RHS|_inf (zetas): " << cg_zetas.error_inf();
                LOG(numerics) << "|Residual|_l2 (zetas): " << cg_zetas.error();
            }

            //__itter_zetas = iterations; // Question: add this to output?
            //__error_zetas = cg_zetas.error_inf();

        }
        updateVDFBudget(iterOffset);
        //LOG(numerics) << "Checking zeta slopes (after zeta height convergence)";
        // checkAllZetaSlopes(); todo remove if not required (in SWI2 used for time-step adjustment)
    }
    updateZetasAfterEquation();

    LOG(numerics) << "Adjusting zeta heights (after zeta height convergence)";
    adjustZetaHeights();
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
