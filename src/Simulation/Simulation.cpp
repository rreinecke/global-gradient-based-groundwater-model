#include "Simulation.hpp"

namespace GlobalFlow {
namespace Simulation {

Simulation::Simulation(Options op, DataReader *reader) : op(op), reader(reader) {
    Eigen::setNbThreads(op.getThreads());

    if (op.cacheEnabled()) {
        serialize = true;
        loadNodes = true;
    }

    //FIXME not pretty could be changed to https://jonasdevlieghere.com/containers-of-unique-pointers/
    //This might be a hughe memory leak at the end :/
    NodeVector ptr(new vector<unique_ptr<Model::NodeInterface>>);
    nodes = std::move(ptr);

    int numOfStaticNodes = initNodes();
    LOG(userinfo) << "Creating Equation..";
    eq = make_unique<GlobalFlow::Solver::Equation>(numOfStaticNodes, nodes, op);
}
}
}
