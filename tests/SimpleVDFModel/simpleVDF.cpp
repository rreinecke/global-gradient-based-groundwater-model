#include "simpleVDF.hpp"


namespace GlobalFlow {

using namespace std;

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleVDF.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleVDFDataReader(op.getStepsizeModifier());
    sim = Simulation::Simulation(op, reader);
    //disabling e-folding
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        sim.getNodes()->at(j)->setSimpleK();
    }
    _eq = sim.getEquation();
}

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_DAYS, 200);
    int stepNumber = 1;

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);

        if (stepNumber == 1) {
            step.first->toggleSteadyState();
        }
        step.first->solve();
        sim.printMassBalances(debug);
        stepNumber++;
    }

    DataProcessing::DataOutput::OutputManager("data/out_simpleVDF.json", sim).write();
    //sim.save();
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