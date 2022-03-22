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
        step.first->toogleSteadyState();
        step.first->solve();
        sim.printMassBalances(userinfo);
    }
    DataProcessing::DataOutput::OutputManager("data/out_simple.json", sim).write();
    //sim.save();
}

void StandaloneRunner::getResults() {

}

void StandaloneRunner::writeData(std::string s) {

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
