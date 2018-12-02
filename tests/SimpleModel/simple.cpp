#include "simple.hpp"


namespace GlobalFlow {

using namespace std;

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simple.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleDataReader(op.getStepsizeModifier());
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::DAY, 1);
    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running a steady state step";
        step.first->toggleSteadyState();
        step.first->solve();
        sim.printMassBalances();
        step.first->toggleSteadyState();
    }

    LOG(userinfo) << "Running a transient steps";
    Simulation::Stepper transientStepper = Simulation::Stepper(_eq, Simulation::DAY, 10);
    for (Simulation::step step : transientStepper) {
        step.first->solve();
        sim.printMassBalances();
    }

    DataProcessing::DataOutput::OutputManager("data/out_simple.json", sim).write();
    delete reader;
}

void StandaloneRunner::getResults() {

}

void StandaloneRunner::writeData(std::string) {

}

StandaloneRunner::StandaloneRunner() {

}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::StandaloneRunner runner;
    runner.loadSettings();
    runner.setupSimulation();
    runner.simulate();
    return EXIT_SUCCESS;
}
