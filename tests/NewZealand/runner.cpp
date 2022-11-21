#include "runner.hpp"

namespace GlobalFlow {

    using namespace std;

    //int sim_id{0};

    void NZRunner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_nz.json");
    }

    void NZRunner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader(op.getStepsizeModifier());
        LOG(debug) << "reader initiated" << std::endl;

        sim = Simulation::Simulation(op, reader);
        LOG(debug) << "sim initiated" << std::endl;

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            sim.getNodes()->at(j)->setSimpleK();
        }
        LOG(debug) << "simple k set for all nodes" << std::endl;

        _eq = sim.getEquation();
    }

    void NZRunner::simulate() {

        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::DAY, 1);
        for (Simulation::step step : stepper) {
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved step with" << step.first->getItter() << "iterations and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);
            step.first->toggleSteadyState();
        }
        sim.save();
    }

    void NZRunner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out_nz.json", sim).write();
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
