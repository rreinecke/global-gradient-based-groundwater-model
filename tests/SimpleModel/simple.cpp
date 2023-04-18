#include "simple.hpp"


namespace GlobalFlow {

using namespace std;

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simple.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleDataReader(op.getStepSizeModifier());
    sim = Simulation::Simulation(op, reader);
    //disabling e-folding
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        sim.getNodes()->at(j)->setSimpleK();
    }
    _eq = sim.getEquation();
}

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::DAY, 1);
    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running a steady state step";
        step.first->toggleSteadyState();
        step.first->solve();
        sim.printMassBalances(debug);
        step.first->toggleSteadyState();
    }

    DataProcessing::DataOutput::OutputManager("data/out_simple.json", sim).write();

    LOG(userinfo) << "Running transient steps";
    Simulation::Stepper transientStepper = Simulation::Stepper(_eq, Simulation::DAY, 10);
    for (Simulation::step step : transientStepper) {
        step.first->solve();
        sim.printMassBalances(debug);
    }

    //Changing stresses
    std::vector<int> ids = {8,
            9,
            18,
            19,
            28,
            29,
            38,
            39,
            48,
            49,
            58,
            59,
            68,
            69,
            78,
            79,
            88,
            89,
            98,
            99};
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        if(std::find(ids.begin(), ids.end(), j) != ids.end())
            sim.getNodes()->at(j)->updateUniqueFlow(0.5, Model::RECHARGE, false);
    }
  
    LOG(userinfo) << "Running transient steps with changed stresses";
    Simulation::Stepper transientStepper2 = Simulation::Stepper(_eq, Simulation::DAY, 10);
    for (Simulation::step step : transientStepper2) {
        step.first->solve();
        sim.printMassBalances(debug);
    }

    DataProcessing::DataOutput::OutputManager("data/out_simple.json", sim).write();
    delete reader;
}

void StandaloneRunner::getResults() {

}

void StandaloneRunner::writeData() {

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