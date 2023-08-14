#include "swi2_ex2.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_swi2_ex2.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SWI2_ex2_DataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}

void StandaloneRunner::writeNodeInfosToCSV(){
    // empty
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
    DataProcessing::DataOutput::OutputManager("data/out_swi2_ex2.json", sim).write();
}

StandaloneRunner::StandaloneRunner() = default;

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::StandaloneRunner runner;
    runner.loadSettings();
    runner.setupSimulation();
    runner.writeNodeInfosToCSV();
    runner.simulate();
    runner.writeData();
    return EXIT_SUCCESS;
}
