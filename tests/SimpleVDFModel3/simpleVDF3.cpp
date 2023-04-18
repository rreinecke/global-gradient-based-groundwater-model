#include "simpleVDF3.hpp"


namespace GlobalFlow {

using namespace std;

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleVDF3.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleVDF3DataReader(op.getStepSizeModifier());
    sim = Simulation::Simulation(op, reader);
    //disabling e-folding
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        sim.getNodes()->at(j)->setSimpleK();
    }
    _eq = sim.getEquation();
}

void StandaloneRunner::simulate() {
    LOG(userinfo) << "Running stress period 1";
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::YEAR, 1);
    int stepNumber = 1;

    // for saving zetas in a csv
    ofstream myfile;
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

    /*
    LOG(userinfo) << "Running stress period 2";
    Simulation::Stepper stepper2 = Simulation::Stepper(_eq, Simulation::YEAR, 1);
    for (Simulation::step step : stepper2) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);
        step.first->solve();
        sim.printMassBalances(debug);
        stepNumber++;
    }*/
}

void StandaloneRunner::getResults() {}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simpleVDF3.json", sim).write();
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
