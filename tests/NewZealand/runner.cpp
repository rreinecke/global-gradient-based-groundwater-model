#include "runner.hpp"

namespace GlobalFlow {

    using namespace std;

    //int sim_id{0};

    void NZRunner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config.json");
    }

    void NZRunner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader(op.getStepsizeModifier());
        sim = Simulation::Simulation(op, reader);
        //disabling e-folding
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            sim.getNodes()->at(j)->setSimpleK();
        }
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
	DataProcessing::DataOutput::OutputManager("data/out_NewZealand.json", sim).write();
    }

    NZRunner::NZRunner() {}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::NZRunner runner;
    runner.loadSettings();
    runner.setupSimulation();
    runner.simulate();
    runner.writeData();
    return EXIT_SUCCESS;
}
