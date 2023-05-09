#include "runner.hpp"

namespace GlobalFlow {

    using namespace std;

    //int sim_id{0};

    void NZRunner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_nzVDF.json");
    }

    void NZRunner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader(op.getStepSizeModifier());
        sim = Simulation::Simulation(op, reader);

        LOG(debug) << "NodeID 1: " << sim.getNodes()->at(1);

        // For neighbour counting:
        ofstream myfile;
        myfile.open ("node_attributes_output.csv");
        myfile << "nodeID,lon,lat,neighbour_count" << std::endl;
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            sim.getNodes()->at(j)->setSimpleK();

            double lon = sim.getNodes()->at(j)->getLon();
            double lat = sim.getNodes()->at(j)->getLat();
            int neig_count = sim.getNodes()->at(j)->getListOfNeighbours().size();

            myfile << j << "," << lon << "," << lat << "," << neig_count << std::endl;
        }
        LOG(debug) << "simple k set for all nodes" << std::endl;
        myfile.close();
        _eq = sim.getEquation();
    }

    void NZRunner::simulate() {

        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::DAY, 1);
        for (Simulation::step step : stepper) {
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved step with " << step.first->getItter() << " iterations and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);
            step.first->toggleSteadyState();
        }
        sim.save();
    }

    void NZRunner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out_nzVDF.json", sim).write();
    }

    NZRunner::NZRunner() {}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::NZRunner runner;
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
