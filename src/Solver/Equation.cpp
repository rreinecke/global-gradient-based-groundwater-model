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
    this->maxAllowedHeadChange = options.getMaxHeadChange();
    this->isAdaptiveDamping = options.isDampingEnabled();
    this->dampMin = options.getMinDamp();
    this->dampMax = options.getMaxDamp();
    this->vdf = options.isDensityVariable(); // vdf = Variable Density Flow
    this->gnc = options.isGridRefined(); // gnc = Ghost Node Correction

    this->inner_iterations = options.getInnerItter();
    this->maxRefinement = options.getMaxRefinement();
    this->nodes = nodes;


    long_vector __x(numberOfNodesTotal);
    long_vector __b(numberOfNodesTotal);

    x = std::move(__x);
    b = std::move(__b);

    Eigen::SparseMatrix<pr_t> __A(numberOfNodesTotal, numberOfNodesTotal);

    A = std::move(__A);
    int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 2 + 1; // +2 for top/down, + 1 for this node
    A.reserve(long_vector::Constant(numberOfNodesTotal, numberOfEntries));

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
        this->maxAllowedZetaChange = options.getMaxZetaChange();

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
//#pragma omp parallel for
    for (const auto &conductance : map) {
        // conductance.first is the nodeID of the respective neighbour node
        NANChecker(conductance.second.value(), "Matrix entry");
        if (A.isCompressed()) {
            A.coeffRef(nodeID, conductance.first) = conductance.second.value();
        } else { // if not compressed: only during matrix initialization (first call of "updateMatrix()")
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
            A_zetas.coeffRef(rowID, colID) = zoneConductance.value();
        }
    }
}

void inline
Equation::updateMatrix() {
    LOG(numerics) << "Updating matrix";
    Index n = A.outerSize();
#ifdef EIGEN_HAS_OPENMP
    Eigen::initParallel(); // initialize Eigen (https://eigen.tuxfamily.org/dox/TopicMultiThreading.html)
    Index threads = Eigen::nbThreads(); // get number of threads specified in config file, and set in Simulation.cpp
#endif
#pragma omp parallel for schedule(dynamic,(n+threads*4-1)/(threads*4)) num_threads(threads)
    for (large_num id = 0; id < numberOfNodesTotal; ++id) {
        //if (std::div(j,100000).rem == 0) {LOG(numerics) << "... reached nodeID " << j;}
        //---------------------Left
        addToA(nodes->at(id));
        //---------------------Right
        b(long(id)) = nodes->at(id)->getRHS().value();
        NANChecker(b[long(id)], "Right hand side");
    }
    if (!A.isCompressed()) {
        LOG(numerics) << "Compressing Matrix";
        A.makeCompressed();
    }
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
    int numberOfEntries = (int) (std::sqrt(maxRefinement) * 4) + 1; // +2 for top/down, + 1 for this node
    A_zetas.reserve(long_vector::Constant(numActive, numberOfEntries));
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
    if (A_zetas.size() != 0 and !A_zetas.isCompressed()) {
        LOG(numerics) << "Compressing Matrix (zetas)";
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
    LOG(numerics) << "Decomposing Matrix (zetas)";
    cg_zetas.compute(A_zetas);
    if (cg_zetas.info() != Success) {
        LOG(userinfo) << "Fail in decomposing matrix (zetas)";
        throw "Fail in decomposing matrix (zetas)";
    }
}

bool inline
Equation::nanInHeadChanges() {
    long_vector changes = adaptiveDamping.getDamping(getResiduals(), x, isAdaptiveDamping);

    bool nanInChanges{false};
#pragma omp parallel for
    for (large_num id = 0; id < numberOfNodesTotal; ++id) {
        if (std::isnan(changes[id])){
            nanInChanges = true;
        }
    }
    return nanInChanges;
}

void inline
Equation::updateIntermediateHeads() {
    long_vector changes = adaptiveDamping.getDamping(getResiduals(), x, isAdaptiveDamping);

#pragma omp parallel for
    for (large_num id = 0; id < numberOfNodesTotal; ++id) {
        // set new head (= old head + change) and headChange
        nodes->at(id)->setHeadChange(changes[id] * si::meter);
    }
}

void inline
Equation::updateIntermediateZetas(large_num iterOffset, int localZetaID) {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
        auto id = index_mapping[k];
        if (id != -1) {
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
        nodes->at(k)->verticalZetaMovement();
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
    updateMatrix(); // updating matrix before iteration
    preconditioner(); // decomposing matrix before iteration
    adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxAllowedHeadChange, x);

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
        double countHeadAbove10k{0.0};
        double changeAtNode{0.0};
        double headAtNode{0.0};
        double nodeID_headMax{0};
        double nodeID_changeMax{0};
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodesTotal; ++k) {
            changeAtNode = std::abs(nodes->at(k)->getProperties().get<quantity<Model::Meter>, Model::HeadChange>().value());
            nodeID_changeMax = (changeAtNode > changeMax) ? k : nodeID_changeMax;
            changeMax = (changeAtNode > changeMax) ? changeAtNode : changeMax;
            headAtNode = std::abs(nodes->at(k)->getProperties().get<quantity<Model::Meter>, Model::Head>().value());
            nodeID_headMax = (headAtNode > headMax) ? k : nodeID_headMax;
            headMax = (headAtNode > headMax) ? headAtNode : headMax;
            headSum += headAtNode;
            if (headAtNode > 10000) {countHeadAbove10k++;}
        }
        headMean = headSum / numberOfNodesTotal;
        LOG(numerics) << "MAX Head: " << headMax << " (nodeID: " << nodeID_headMax << ")";
        LOG(numerics) << "MEAN Head: " << headMean << ", Head > 10,000 at " << countHeadAbove10k << " nodes";
        maxHeadChange = changeMax;
        LOG(numerics) << "MAX Head Change: " << changeMax << " (nodeID: " << nodeID_changeMax << ")";
        return changeMax > maxAllowedHeadChange;
    };

    long iterations{0};
    bool headFail{false};
    char smallHeadChanges{0};
    bool headConverged{false};
    int nanInHeadChangeCounter{0};
    while (iterations < IITER) {
        LOG(numerics) << "Outer iteration: " << iterations;
        auto x_t0 = x; // save heads of previous outer iteration
        x = cg.solveWithGuess(b, x); // solving inner iterations
        //LOG(debug) << "A:\n" << A << std::endl;
        //LOG(debug) << "x:\n" << x << std::endl;
        //LOG(debug) << "b (= rhs):\n" << b << std::endl;
        if (nanInHeadChanges()){  // if there is a nan value in calculated head change of any node
            LOG(numerics) << "Nan in head changes. Trying to solve the same outer iteration again.";
            x = x_t0; // reset heads back to heads of previous outer iteration
            nanInHeadChangeCounter++;
            // if there were nans in head change in several (e.g., 5) outer iteration trials
            if (nanInHeadChangeCounter >= 5) { // todo move to config
                throw "Fail: nan in head change solution is persistent.";
            }
            continue;
        }

        updateIntermediateHeads();

        int innerItter{0};
        innerItter = cg.iterations();
        if (innerItter == 0 and iterations == 0) {
            LOG(numerics) << "Convergence criterion met without iterations. ";
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
                    break;
                }
            }
        }

        if(maxHeadChange == oldMaxHeadChange){
            //The head change is really the same -> increase inner iterations
            itterScale = itterScale + 10;
        }else{
            itterScale = 0;
        }
        cg.setMaxIterations(inner_iterations + itterScale);
        oldMaxHeadChange = maxHeadChange;

        /**
         * @brief residual norm convergence
         */
        LOG(numerics) << "Inner iterations: " << innerItter;
        if (cg.info() == Success and iterations != 0) {
            LOG(numerics) << "cg solver success"; // this is always reached in outer iteration 1 (starts at 0)
            break;
        }
        LOG(numerics) << "Head change bigger: " << headFail;

        updateMatrix();
        preconditioner(); // decomposing matrix during iteration

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

    __itter = iterations;
    __error = cg.error_inf();

    LOG(numerics) << "Updating heads";
    updateFinalHeads();

    /**
     * ###############################
     * # Solve Zeta Surface Equation #
     * ###############################
     */
     if(vdf) {
         solve_zetas();
     }

     /*if(vdf) {
         LOG(numerics) << "Adjusting zeta heights (after zeta height convergence)";
         adjustZetaHeights();
         LOG(numerics) << "Updating zetasTZero";
         updateZetas_TZero();
     }*/
    /**
    * ###############################
    * # Update budgets #
    * ###############################
    */
     if(gnc) {
         updateGNCBudget();
     }
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
        large_num iterOffset = layer * numberOfNodesPerLayer;
        LOG(numerics) << "Finding zeta surface heights in layer " << layer;
        for (int localZetaID = 1; localZetaID < numberOfZones; localZetaID++) {
            LOG(numerics) << "Solving for zeta surface " << localZetaID << "";
            LOG(numerics) << "Updating Matrix (zetas)";
            updateMatrix_zetas(iterOffset, localZetaID);

            if (A_zetas.size() == 0){ // if matrix empty continue with next iteration
                continue;
            }

            //LOG(debug) << "A_zetas (before preconditioner):\n" << A_zetas << std::endl;
            preconditioner_zetas();

            double maxZetaChange{0};
            double oldMaxZetaChange{0};
            int itterScale{0};

            // Returns true if max zetachange is greater than defined val
            auto isZetaChangeGreater = [this, &maxZetaChange, &iterOffset, &localZetaID]() -> bool {
                double zetaChangeMax{0.0};
                double zetaChangeNode{0.0};
                double zetaAtNode{0.0};
                double zetaMax{0.0};
#pragma omp parallel for
                for (large_num k = 0; k < numberOfNodesPerLayer; ++k) {
                    zetaChangeNode = std::abs(nodes->at(k + iterOffset)->getZetaChange(localZetaID).value());
                    zetaChangeMax = (zetaChangeNode > zetaChangeMax) ? zetaChangeNode : zetaChangeMax;
                    zetaAtNode = std::abs(nodes->at(k)->getProperties().get<quantity<Model::Meter>, Model::Head>().value());
                    zetaMax = (zetaAtNode > zetaMax) ? zetaAtNode : zetaMax;
                }
                LOG(numerics) << "MAX Zeta: " << zetaMax;
                maxZetaChange = zetaChangeMax;
                LOG(numerics) << "MAX Zeta Change: " << zetaChangeMax;
                return zetaChangeMax > maxAllowedZetaChange;
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
                updateIntermediateZetas(iterOffset, localZetaID);

                int innerItter = (int) cg_zetas.iterations();

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

                if (maxZetaChange == oldMaxZetaChange) {
                    //The zeta change is really the same -> increase inner iterations
                    itterScale = itterScale + 10;
                } else {
                    itterScale = 0;
                }
                cg_zetas.setMaxIterations(inner_iterations + itterScale);
                oldMaxZetaChange = maxZetaChange;

                /**
                 * @brief residual norm convergence // Question: make function of this?
                 */
                LOG(numerics) << "Inner iterations (zetas): " << innerItter;
                if (cg_zetas.info() == Success and iterations != 0) {
                    LOG(numerics) << "cg_zetas success";
                    break;
                }
                //LOG(numerics) << "|Residual|_inf / |RHS|_inf (zetas): " << cg_zetas.error_inf();
                //LOG(numerics) << "|Residual|_l2 (zetas): " << cg_zetas.error();
                LOG(numerics) << "Zeta change bigger: " << zetaFail;

                updateMatrix_zetas(layer, localZetaID);

                // if matrix is empty, go to next iteration
                if (A_zetas.size() == 0) { continue; }

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
            //LOG(debug) << "x_zetas[" << localZetaID << "] on node layer " << layer << ":\n" << x_zetas << std::endl;

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
