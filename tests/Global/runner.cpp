#include "runner.hpp"

namespace GlobalFlow {

    using namespace std;

    //int sim_id{0};

    void GlobalRunner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_global.json");
    }

    void GlobalRunner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        LOG(debug) << "reader initiated" << std::endl;

        sim = Simulation::Simulation(op, reader);
        LOG(debug) << "sim initiated" << std::endl;

        LOG(debug) << "NodeID 1: " << sim.getNodes()->at(1);

        // For node infos:
        ofstream myfile;
        myfile.open ("node_attributes_output.csv");
        myfile << "nodeID,lon,lat,neighbour_count,elevation,hyd_cond,recharge" << std::endl;

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            sim.getNodes()->at(j)->setSimpleK();

            myfile << j << "," <<
                sim.getNodes()->at(j)->getLon() << "," <<
                sim.getNodes()->at(j)->getLat() << "," <<
                sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                sim.getNodes()->at(j)->getElevation().value() << "," <<
                sim.getNodes()->at(j)->getK().value() << "," <<
                sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() <<
                std::endl;
        }
        LOG(debug) << "simple k set for all nodes" << std::endl;
        myfile.close();
        _eq = sim.getEquation();
    }

    void GlobalRunner::simulate() {

        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::MONTH, 1);
        for (Simulation::step step : stepper) {
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved step with " << step.first->getItter() << " iterations and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);
            step.first->toggleSteadyState();
        }
        sim.save();
    }

    void GlobalRunner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out_global.json", sim).write();
    }

    GlobalRunner::GlobalRunner() {}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::GlobalRunner runner;
    std::cout << FBLU("- loading settings...") << std::endl;
    runner.loadSettings();
    std::cout << FBLU("- simulation setup...") << std::endl;
    runner.setupSimulation();
    std::cout << FBLU("- simulating...") << std::endl;
    runner.simulate();
    std::cout << FBLU("- writing output...") << std::endl;
    runner.writeData();
    return EXIT_SUCCESS;
}
