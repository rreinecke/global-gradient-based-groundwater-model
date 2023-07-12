#include "simpleVDF2.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleVDF2.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleVDF2DataReader();
    sim = Simulation::Simulation(op, reader);

    LOG(debug) << sim.getNodes()->at(1);

    //disabling e-folding
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        sim.getNodes()->at(j)->setSimpleK();
    }
    _eq = sim.getEquation();
}

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_DAYS, 1000);
    int stepNumber = 1;

    // for saving zetas in a csv
    std::ofstream myfile;
    myfile.open ("zetas.csv");
    myfile << "timestep,nodeID,zetaID,zeta" << std::endl;

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);

        if (stepNumber == 1) {
            step.first->toggleSteadyState();
        }
        step.first->solve();
        sim.printMassBalances(debug);

        // for saving zetas in a csv
        int zetaID;
        double zeta;
        zetaID = 1;
        for (int nodeID = 0; nodeID < sim.getNodes()->size(); ++nodeID) {
            zeta = sim.getNodes()->at(nodeID)->getZeta(zetaID).value();
            myfile << stepNumber << "," << nodeID << "," << zetaID << "," << zeta << std::endl;
        }
        zetaID = 2;
        for (int nodeID = 0; nodeID < sim.getNodes()->size(); ++nodeID) {
            zeta = sim.getNodes()->at(nodeID)->getZeta(zetaID).value();
            myfile << stepNumber << "," << nodeID << "," << zetaID << "," << zeta << std::endl;
        }

        stepNumber++;
    }
    myfile.close(); // for saving zetas in a csv
    //sim.save();
}

void StandaloneRunner::getResults() {}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simpleVDF2.json", sim).write();
}

StandaloneRunner::StandaloneRunner() {}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::StandaloneRunner runner;
    runner.loadSettings();
    runner.setupSimulation();
    runner.simulate();
    runner.writeData();
    return EXIT_SUCCESS;
}
