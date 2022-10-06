#include "Equation.hpp"

namespace GlobalFlow {
namespace Solver {

Equation::Equation(large_num numberOfNodes, NodeVector nodes, Simulation::Options options) : options(options) {
    LOG(userinfo) << "Setting up Equation for " << numberOfNodes << std::endl;

    this->numberOfNodes = numberOfNodes;
    this->IITER = options.getMaxIterations();
    this->RCLOSE = options.getConverganceCriteria();
    this->initialHead = options.getInitialHead();
    this->maxHeadChange = options.getMaxHeadChange();
    this->isAdaptiveDamping = options.isDampingEnabled();
    this->dampMin = options.getMinDamp(); // Question: introduce new min/max damping values for zeta surfaces?
    this->dampMax = options.getMaxDamp();
    this->disable_dry_cells = options.disableDryCells();
    this->nwt = (options.getSolverName().compare("NWT") == 0);
    if (nwt)
        LOG(userinfo) << "Running with NWT solver" << std::endl;
    this->vdf = options.isDensityVariable();

    this->inner_iterations = options.getInnerItter();

    this->nodes = nodes;


    long_vector __x(numberOfNodes);
    long_vector __b(numberOfNodes);

    x = std::move(__x);
    b = std::move(__b);

    Eigen::SparseMatrix<pr_t> __A(numberOfNodes, numberOfNodes);

    A = std::move(__A);
    A.reserve(long_vector::Constant(numberOfNodes, 7));

    //Init first result vector x by writing initial heads
    //Initial head should be positive
    //resulting head is the real hydraulic head
    double tmp = 0;
#pragma omp parallel for
    for (int i = 0; i < numberOfNodes; ++i) {
        if (nwt) {
            nodes->at(i)->enableNWT();
        }
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

    if (nwt) {
        bicgstab.setMaxIterations(inner_iterations);
        bicgstab.setTolerance(RCLOSE);
    } else {
        //set inner iterations
        cg.setMaxIterations(inner_iterations);
        cg.setTolerance(RCLOSE);
        //cg.preconditioner().setInitialShift(1e-8);
    }

    if(vdf) {
        LOG(userinfo) << "Simulating variable density flow" << std::endl;
        this->numberOfZones = options.getNumberOfDensityZones();
        this->totalNumberOfZetas = numberOfNodes * (numberOfZones - 1);
        this->maxZetaChange = options.getMaxZetaChange();
        long_vector __x_zetas(totalNumberOfZetas);
        long_vector __b_zetas(totalNumberOfZetas);

        x_zetas = std::move(__x_zetas);
        b_zetas = std::move(__b_zetas);

        Eigen::SparseMatrix<pr_t> __A_zetas(totalNumberOfZetas, totalNumberOfZetas);

        A_zetas = std::move(__A_zetas);
        A_zetas.reserve(long_vector::Constant(totalNumberOfZetas, 5));

        if (nwt) {
            bicgstab_zetas.setMaxIterations(inner_iterations);
            bicgstab_zetas.setTolerance(RCLOSE);
        } else {
            //set inner iterations
            cg_zetas.setMaxIterations(inner_iterations);
            cg_zetas.setTolerance(RCLOSE);
        }
    }
}

Equation::~Equation() {
    LOG(debug) << "Destroying equation\n" << std::endl;
}

void inline
Equation::addToA(std::unique_ptr<Model::NodeInterface> const &node, bool cached) {
    std::unordered_map<large_num, quantity<Model::MeterSquaredPerTime>> map;
    if (nwt) {
        map = node->getJacobian();
    } else {
        map = node->getConductance();
    }
    auto nodeID = node->getProperties().get<large_num, Model::ID>();

    if (disable_dry_cells) {
        if (not map[nodeID].value()) {
            disabled_nodes.insert(nodeID);
        }
    }

    for (const auto &conductance : map) {
        // conductance.first is the nodeID of the respective neighbour node
        if (cached) {
            A.coeffRef(nodeID, conductance.first) = conductance.second.value();
        } else {
            A.insert(nodeID, conductance.first) = conductance.second.value();
        }
    }
}

void inline
Equation::addToA_zeta(std::unique_ptr<Model::NodeInterface> const &node, int localZetaID, large_num globalZetaID, bool cached) {
    std::unordered_map<large_num, quantity<Model::MeterSquaredPerTime>> map;
    quantity<Model::MeterSquaredPerTime> zoneConductance;

    map = node->getConductance_ZetaMovement(localZetaID);

    for (const auto &entry : map) { // entry contains: [1] node id of the horizontal neighbours, [2] conductance of zone n
        zoneConductance = entry.second;
        //LOG(debug) << "globalZetaID row: " + std::to_string(globalZetaID) + "| col: " + std::to_string(entry.first) + " (in addToA_zetas)" << std::endl;
        // conductance.first is the zetaID of the zeta surface in the respective neighbour node
        if (cached) {
            A_zetas.coeffRef(globalZetaID - numberOfNodes, entry.first - numberOfNodes) = zoneConductance.value();
        } else {
            A_zetas.insert(globalZetaID - numberOfNodes, entry.first - numberOfNodes) = zoneConductance.value();
        }
    }
}

void inline Equation::reallocateMatrix() {
    if (not disable_dry_cells) {
        return;
    }

    if (disabled_nodes.empty()) {
        _A_ = A;
        _b_ = b;
        _x_ = x;
        return;
    }

    large_num __missing{0};
    long size = A.rows() - disabled_nodes.size();
    Matrix<pr_t, 2, Dynamic> new_matrix(size, size);


    if (dry_have_changed) {
        index_mapping.clear();
    }

    for (large_num i = 0; i < A.rows(); i++) {
        if (not dry_have_changed) {
            auto m = index_mapping[i];
            if (m != -1) {
                new_matrix.row(i) = A.row(i);
            }
        } else {
            const bool is_in = disabled_nodes.find(i) != disabled_nodes.end();
            if (not is_in) {
                //Save all non-zero
                new_matrix.row(i) = A.row(i);
                index_mapping[i] = i - __missing;
            } else {
                __missing++;
                index_mapping[i] = -1;
            }
        }
    }
    _A_ = new_matrix.sparseView();

    long_vector buffer_b;
    long_vector buffer_x;
    for (int j = 0; j < b.size(); ++j) {
        auto m = index_mapping[j];
        if (m != -1) {
            buffer_b[m] = b[j];
            if (dry_have_changed) {
                buffer_x[m] = x[j];
            }
        }
    }
    _b_ = buffer_b;
    if (dry_have_changed) {
        _x_ = buffer_x;
    }
}

void inline Equation::reallocateMatrix_zetas() {
    if (not disable_dry_cells) {
        return;
    }

    if (disabled_nodes.empty()) { // todo disabled zetas?
        _A__zetas = A_zetas;
        _b__zetas = b_zetas;
        _x__zetas = x_zetas;
        return;
    }

    large_num __missing{0};
    long size = A_zetas.rows() - disabled_nodes.size();
    Matrix<pr_t, 2, Dynamic> new_matrix_zetas(size, size);


    if (dry_have_changed) {
        index_mapping.clear(); // todo index mappping zetas?
    }

    for (large_num i = 0; i < A_zetas.rows(); i++) {
        if (not dry_have_changed) {
            auto m = index_mapping[i];
            if (m != -1) {
                new_matrix_zetas.row(i) = A_zetas.row(i);
            }
        } else {
            const bool is_in = disabled_nodes.find(i) != disabled_nodes.end();
            if (not is_in) {
                //Save all non-zero
                new_matrix_zetas.row(i) = A_zetas.row(i);
                index_mapping[i] = i - __missing;
            } else {
                __missing++;
                index_mapping[i] = -1;
            }
        }
    }
    _A__zetas = new_matrix_zetas.sparseView();

    long_vector buffer_b_zetas;
    long_vector buffer_x_zetas;
    for (int j = 0; j < b_zetas.size(); ++j) {
        auto m = index_mapping[j];
        if (m != -1) {
            buffer_b_zetas[m] = b_zetas[j];
            if (dry_have_changed) {
                buffer_x_zetas[m] = x_zetas[j];
            }
        }
    }
    _b__zetas = buffer_b_zetas;
    if (dry_have_changed) {
        _x__zetas = buffer_x_zetas;
    }
}

void inline
Equation::updateMatrix() {
    Index n = A.outerSize();
    std::unordered_set<large_num> tmp_disabled_nodes = disabled_nodes;
    disabled_nodes.clear();

#ifdef EIGEN_HAS_OPENMP
    Eigen::initParallel();
    Index threads = Eigen::nbThreads();
#endif
#pragma omp parallel for schedule(dynamic,(n+threads*4-1)/(threads*4)) num_threads(threads)
    for (large_num j = 0; j < numberOfNodes; ++j) {
        //---------------------Left
        addToA(nodes->at(j), isCached);
        //---------------------Right
        double rhs{0.0};
        if (nwt) {
            rhs = nodes->at(j)->getRHS__NWT();
            b[nodes->at(j)->getProperties().get<large_num, Model::ID>()] = std::move(rhs);
            NANChecker(rhs, "Right hand side");
        } else {
            rhs = nodes->at(j)->getRHS().value();
            b(nodes->at(j)->getProperties().get<large_num, Model::ID>()) = rhs;
            NANChecker(rhs, "Right hand side");
        }
    }

    if (not disable_dry_cells) {
        dry_have_changed = not set_compare(tmp_disabled_nodes, disabled_nodes);
    }

    //some nodes are in dried condition
    //reallocate matrix, and a,b
    reallocateMatrix();

    //A is up to date and b contains only rhs term
    //Update with current heads
    if (nwt) {
        if (disable_dry_cells) {
            _b_ = _b_ + _A_ * _x_;
        } else {
            b = b + A * x;
        }
    }

    //Check if after iteration former 0 values turned to non-zero
    if (disable_dry_cells) {
        if ((not _A_.isCompressed()) and isCached) {
            LOG(numerics) << "Recompressing Matrix";
            _A_.makeCompressed();
        }
    } else {
        if ((not A.isCompressed()) and isCached) {
            LOG(numerics) << "Recompressing Matrix";
            A.makeCompressed();
        }
    }
    //LOG(debug) << "A:\n" << A << std::endl;
    //LOG(debug) << "x_zetas:\n" << x_zetas << std::endl;
    //LOG(debug) << "b (= rhs):\n" << b << std::endl;
}

void inline
Equation::updateMatrix_zetas(int localZetaID) { // todo: adapt for multiple layers and more than one moving zeta surface (each has one equation system)
        Index n = A_zetas.outerSize();
        std::unordered_set<large_num> tmp_disabled_nodes = disabled_nodes;
        disabled_nodes.clear();

        large_num globalZetaID = numberOfNodes; // because we start with localZetaID = 1 (and not 0) todo this needs to be documented well (implications for initial_zetas.csv)
#ifdef EIGEN_HAS_OPENMP
        Eigen::initParallel();
        Index threads = Eigen::nbThreads();
#endif
#pragma omp parallel for schedule(dynamic,(n+threads*4-1)/(threads*4)) num_threads(threads)
        for (large_num j = 0; j < numberOfNodes; ++j) {
            //---------------------Left
            addToA_zeta(nodes->at(j), localZetaID, globalZetaID, isCached);
            //---------------------Right
            b_zetas(globalZetaID - numberOfNodes) = nodes->at(j)->getRHS_zeta(localZetaID).value();
            globalZetaID++;
        }

        if (not disable_dry_cells) {
            dry_have_changed = not set_compare(tmp_disabled_nodes, disabled_nodes);
        }

        //some nodes are in dried condition
        //reallocate matrix, and a,b
        reallocateMatrix_zetas(); // Question: is this needed with zetas?

        //A is up to date and b contains only rhs term
        //Update with current heads
        if (nwt) {
            if (disable_dry_cells) {
                _b__zetas = _b__zetas + _A__zetas * _x__zetas;
            } else {
                b_zetas = b_zetas + A_zetas * x_zetas;
            }
        }

        //Check if after iteration former 0 values turned to non-zero
        if (disable_dry_cells) {
            if ((not _A__zetas.isCompressed()) and isCached) {
                LOG(numerics) << "Recompressing Matrix";
                _A__zetas.makeCompressed();
            }
        } else {
            if ((not A_zetas.isCompressed()) and isCached) {
                LOG(numerics) << "Recompressing Matrix";
                A_zetas.makeCompressed();
            }
        }
    LOG(debug) << "A_zetas:\n" << A_zetas << std::endl;
    LOG(debug) << "b_zetas (= rhs_zeta):\n" << b_zetas << std::endl;
    }

void inline
Equation::preconditioner() {
    LOG(numerics) << "Decomposing Matrix";
    if (nwt) {
        if (disable_dry_cells) {
            bicgstab.compute(_A_);
        } else {
            bicgstab.compute(A);
        }
        if (bicgstab.info() != Success) {
            LOG(numerics) << "Fail in decomposing matrix";
            throw "Fail in decomposing matrix";
        }
    } else {
        if (disable_dry_cells) {
            cg.compute(_A_);
        } else {
            cg.compute(A);
        }
        if (cg.info() != Success) {
            LOG(numerics) << "Fail in decomposing matrix";
            throw "Fail in decomposing matrix";
        }
    }
}

void inline
Equation::preconditioner_zetas() {
    LOG(numerics) << "Decomposing Matrix (zetas)";
    if (nwt) {
        if (disable_dry_cells) {
            bicgstab_zetas.compute(_A__zetas);
        } else {
            bicgstab_zetas.compute(A_zetas);
        }
        if (bicgstab_zetas.info() != Success) {
            LOG(numerics) << "Fail in decomposing matrix (zetas)";
            throw "Fail in decomposing matrix (zetas)";
        }
    } else {
        if (disable_dry_cells) {
            cg_zetas.compute(_A__zetas);
        } else {
            cg_zetas.compute(A_zetas);
        }
        if (cg_zetas.info() != Success) {
            LOG(numerics) << "Fail in decomposing matrix (zetas)";
            throw "Fail in decomposing matrix (zetas)";
        }
    }
}

void inline
Equation::updateIntermediateHeads() {
    long_vector changes;
    if (disable_dry_cells) {
        changes = adaptiveDamping.getDamping(getResiduals(), _x_, isAdaptiveDamping);
    } else {
        changes = adaptiveDamping.getDamping(getResiduals(), x, isAdaptiveDamping);
    }
    LOG(debug) << "Head changes (updateIntermediateHeads):\n" << changes << std::endl;

    bool reduced = disabled_nodes.empty();
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodes; ++k) {
        large_num id = nodes->at(k)->getProperties().get<large_num, Model::ID>();

        if (reduced) {
            nodes->at(k)->setHead((double) changes[id] * si::meter);
        } else {
            auto m = index_mapping[id];
            if (m != -1) {
                nodes->at(k)->setHead((double) changes[m] * si::meter);
            }
        }
    }
}

void inline
Equation::updateIntermediateZetas(int localZetaID) {
    long_vector changes;
    if (disable_dry_cells) {
        changes = adaptiveDamping.getDamping(getResiduals_zetas(), _x__zetas, isAdaptiveDamping);
    } else {
        changes = adaptiveDamping.getDamping(getResiduals_zetas(), x_zetas, isAdaptiveDamping);
    }
    LOG(debug) << "Zeta changes (updateIntermediateZetas):\n" << changes << std::endl;
    bool reduced = disabled_nodes.empty();

#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodes; ++k) {
        large_num id = nodes->at(k)->getProperties().get<large_num, Model::ID>();
        large_num globalZetaID = nodes->at(k)->getGlobalZetaID(localZetaID);

        if (reduced) {
            nodes->at(k)->addDeltaToZeta(localZetaID, (double) changes[globalZetaID-numberOfNodes] * si::meter); // setZeta
        } else {
            auto m = index_mapping[id];
            if (m != -1) {
                nodes->at(k)->addDeltaToZeta(localZetaID, (double) changes[globalZetaID-numberOfNodes] * si::meter); // setZeta
            }
        }
    }
}

void inline
Equation::updateFinalHeads() {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodes; ++k) {
        nodes->at(k)->updateHeadChange();
    }
}

void inline
Equation::updateFinalZetaChange(int localZetaID) {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodes; ++k) {
            nodes->at(k)->updateZetaChange(localZetaID);
        }
}

void inline
Equation::updateTopZetasToHeads() {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodes; ++k) {
            nodes->at(k)->clipTopZetasToHeads();
        }
}

/*void inline // todo: remove if not required (in SWI2: required for time step adjustment)
Equation::checkAllZetaSlopes(int localZetaID) {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodes; ++k) {
            nodes->at(k)->checkZetaSlopes();
        }
}*/

void inline
Equation::adjustAllZetaHeights(int localZetaID) {
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodes; ++k) {
            nodes->at(k)->adjustZetaHeights(localZetaID);
        }
}

void inline
Equation::updateBudget() {
#pragma omp parallel for
    for (large_num k = 0; k < numberOfNodes; ++k) {
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
        if (disable_dry_cells) {
            _A_.makeCompressed();
        } else {
            A.makeCompressed();
        }
        LOG(numerics) << "Cached Matrix";
        isCached = true;
    }

    preconditioner();
    if (disable_dry_cells) {
        adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxHeadChange, _x_);
    } else {
        adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxHeadChange, x);
    }

    LOG(numerics) << "Running Time Step";

    double maxHead{0};
    double oldMaxHead{0};
    int itterScale{0};

    // Returns true if max headchange is greater than defined val
    auto isHeadChangeGreater = [this,&maxHead]() -> bool {
        double lowerBound = maxHeadChange;
        double changeMax = 0;
#pragma omp parallel for
        for (large_num k = 0; k < numberOfNodes; ++k) {
            double
                    val = std::abs(
                    nodes->at(k)->getProperties().get<quantity<Model::Meter>, Model::HeadChange>().value());
            //LOG(debug) << "val (in solve): " << val << std::endl;
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

        //Solve inner iterations
        if (nwt) {
            if (disable_dry_cells) {
                _x_ = bicgstab.solveWithGuess(_b_, _x_);
            } else {
                x = bicgstab.solveWithGuess(b, x);
            }
        } else {
            if (disable_dry_cells) {
                _x_ = cg.solveWithGuess(_b_, _x_);
            } else {
                x = cg.solveWithGuess(b, x);
            }
        }
        LOG(debug) << "x (potential new hydraulic head) (outer iteration " << iterations << "):\n" << x << std::endl;

        updateIntermediateHeads();
        int innerItter{0};
        if (nwt) {
            innerItter = bicgstab.iterations();
        } else {
            innerItter = cg.iterations();
        }

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
        if (nwt) {
            LOG(numerics) << "Inner iterations: " << innerItter;
            if (bicgstab.info() == Success) {
                LOG(numerics) << "bicgstab solver success";
                break;
            }
        } else {
            LOG(numerics) << "Inner iterations: " << innerItter;
            if (cg.info() == Success and iterations != 0) {
                LOG(numerics) << "cg solver success"; // this is always reached in outer iteration 1 (starts at 0)

                // Question: why dont we iterate until convergence?? Or: what does cg.info() == Success mean?
                break;
            }
        }

        if (nwt) {
            LOG(numerics) << "Residual squared norm error: " << bicgstab.error();
            LOG(numerics) << "Head change bigger: " << headFail;
        } else {
            LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
	        LOG(numerics) << "|Residual|_l2: " << cg.error();
            LOG(numerics) << "Head change bigger: " << headFail;
        }

        updateMatrix();
        preconditioner();

        iterations++;
    }

    if (iterations == IITER) {
        std::cerr << "Fail in solving matrix with max iterations\n";
        if (nwt) {
            LOG(numerics) << "Residual squared norm error: " << bicgstab.error();
        } else {
            LOG(numerics) << "|Residual|_inf / |RHS|_inf: " << cg.error_inf();
            LOG(numerics) << "|Residual|_l2: " << cg.error();
        }
    }

    updateFinalHeads();
    updateBudget();

    __itter = iterations;
    __error = nwt ? bicgstab.error() : cg.error_inf();


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

    for (int localZetaID = 1; localZetaID < numberOfZones; localZetaID++) {
        //Init first result vector x_zetas by writing initial zetas
#pragma omp parallel for
        for (int i = 0; i < numberOfNodes; ++i) {
            x_zetas[nodes->at(i)->getProperties().get<large_num, Model::ID>()] =
                    nodes->at(i)->getZetas()[localZetaID].value();
            //nodes->at(i)->initHead_t0();
        }

        LOG(numerics) << "Updating Matrix (zeta surface " << localZetaID <<")";
        updateMatrix_zetas(localZetaID);

        if (!isCached_zetas) {
            LOG(numerics) << "Compressing matrix (zetas)";
            if (disable_dry_cells) {
                _A__zetas.makeCompressed();
            } else {
                A_zetas.makeCompressed();
            }
            LOG(numerics) << "Cached Matrix (zetas)";
            isCached_zetas = true;
        }

        preconditioner_zetas();
        if (disable_dry_cells) {
            adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxZetaChange, _x__zetas);
        } else {
            adaptiveDamping = AdaptiveDamping(dampMin, dampMax, maxZetaChange, x_zetas);
        }

        LOG(numerics) << "Running Time Step (zetas)";

        double maxZeta{0};
        double oldMaxZeta{0};
        int itterScale{0};

        // Returns true if max zetachange is greater than defined val
        auto isZetaChangeGreater = [this, &maxZeta]() -> bool {
            double lowerBound = maxZetaChange;
            double changeMax = 0;
#pragma omp parallel for
            for (large_num k = 0; k < totalNumberOfZetas; ++k) {
                double val;
                for (int localZetaID = 1; localZetaID < numberOfZones; localZetaID++) { // need to define localZetaID inside "isZetaChangeGreater"
                    val = std::abs(nodes->at(k)->getZetasChange()[localZetaID].value()); // todo improve for loop (by getting rid of it)
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

        while (iterations < IITER) {

            LOG(numerics) << "Outer iteration (zetas): " << iterations;

            //Solve inner iterations
            if (nwt) {
                if (disable_dry_cells) {
                    _x__zetas = bicgstab_zetas.solveWithGuess(_b__zetas, _x__zetas);
                } else {
                    x_zetas = bicgstab_zetas.solveWithGuess(b_zetas, x_zetas);
                }
            } else {
                if (disable_dry_cells) {
                    _x__zetas = cg_zetas.solveWithGuess(_b__zetas, _x__zetas);
                } else {
                    x_zetas = cg_zetas.solveWithGuess(b_zetas, x_zetas);
                }
            }
            LOG(debug) << "x_zetas (potential new zeta heights) before updating and convergence check (outer iteration "
                       << iterations << "):\n" << x_zetas << std::endl;

            updateIntermediateZetas(localZetaID);
            int innerItter{0};
            if (nwt) {
                innerItter = bicgstab_zetas.iterations();
            } else {
                innerItter = cg_zetas.iterations();
            }

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
             * @brief residual norm convergence // todo make function of this (need to
             */
            if (nwt) {
                LOG(numerics) << "Inner iterations (zetas): " << innerItter;
                if (bicgstab_zetas.info() == Success) {
                    LOG(numerics) << "bicgstab_zetas success";
                    break;
                }
            } else {
                LOG(numerics) << "Inner iterations (zetas): " << innerItter;
                if (cg_zetas.info() == Success and iterations != 0) {
                    LOG(numerics) << "cg_zetas success";
                    break;
                }
            }

            if (nwt) {
                LOG(numerics) << "Residual squared norm error (zetas)" << bicgstab_zetas.error();
                LOG(numerics) << "Zeta change bigger: " << zetaFail;
            } else {
                LOG(numerics) << "|Residual|_inf / |RHS|_inf (zetas): " << cg_zetas.error_inf();
                LOG(numerics) << "|Residual|_l2 (zetas): " << cg_zetas.error();
                LOG(numerics) << "Zeta change bigger: " << zetaFail;
            }

            updateMatrix_zetas(localZetaID);
            preconditioner_zetas();

            iterations++;

        }

        if (iterations == IITER) {
            std::cerr << "Fail in solving matrix with max iterations (zetas)\n";
            if (nwt) {
                LOG(numerics) << "Residual squared norm error (zetas): " << bicgstab_zetas.error();
            } else {
                LOG(numerics) << "|Residual|_inf / |RHS|_inf (zetas): " << cg_zetas.error_inf();
                LOG(numerics) << "|Residual|_l2 (zetas): " << cg_zetas.error();
            }
        }

        //LOG(numerics) << "Checking zeta slopes (after zeta height convergence)";
        // checkAllZetaSlopes(); todo remove if not required (in SWI2 used for time-step adjustment)

        LOG(numerics) << "Adjusting zeta heights (after zeta height convergence)";
        adjustAllZetaHeights(localZetaID);
        // updateZetaBudget(); // Question: calculate zeta budgets?

        updateFinalZetaChange(localZetaID);

        __itter = iterations;
        __error = nwt ? bicgstab_zetas.error() : cg_zetas.error_inf();
    }
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