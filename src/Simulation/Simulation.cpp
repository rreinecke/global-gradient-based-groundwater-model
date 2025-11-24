#include "Simulation.hpp"

BOOST_CLASS_IMPLEMENTATION(GlobalFlow::Model::NodeInterface, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(GlobalFlow::Model::StandardNode, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(GlobalFlow::Model::StaticHeadNode, boost::serialization::object_serializable);

namespace GlobalFlow {
namespace Simulation {

Simulation::Simulation(Options op, DataReader *reader) : op(op), reader(reader) {
    Eigen::setNbThreads(op.getThreads());

    if (op.cacheEnabled()) {
        serialize = true;
        loadNodes = true;
    }

    //FIXME not pretty could be changed to https://jonasdevlieghere.com/containers-of-unique-pointers/
    //This might be a huge memory leak at the end :/
    NodeVector ptr(new  std::vector< std::unique_ptr<Model::NodeInterface>>);
    nodes = std::move(ptr);
    int numOfStaticNodes{0};

    if(loadNodes){
        LOG(stateinfo) << "Attempting to load old state";
        if(boost::filesystem::exists(saveName)){
            restore();
            succefullyRestored = true;
        }else{
            LOG(userinfo) << "No existing state to load";
            LOG(userinfo) << "Starting new run";
            initNodes();
        }
    } else{
         initNodes();
    }
    LOG(userinfo) << "Creating Equation..";
    eq = make_unique<GlobalFlow::Solver::Equation>(nodes, op);
}
}
}
