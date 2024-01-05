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

    void StandaloneRunner::writeNodeInfosToCSV() {
        // For node infos:
        std::ofstream myfile;
        myfile.open ("node_attributes_swi2.csv");
        myfile << "nodeID,spatID,lon,lat,neighbour_count,elevation,bottom,hyd_cond,zeta[1]" << std::endl;

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            myfile << j << "," <<
                   sim.getNodes()->at(j)->getSpatID() << "," <<
                   sim.getNodes()->at(j)->getLon() << "," <<
                   sim.getNodes()->at(j)->getLat() << "," <<
                   sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                   sim.getNodes()->at(j)->getElevation().value() << "," <<
                   sim.getNodes()->at(j)->getBottom().value() << "," <<
                   sim.getNodes()->at(j)->getK().value() << "," <<
                   sim.getNodes()->at(j)->getZeta(1).value() <<
                   std::endl;
        }
    }

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_DAYS, 1000);
    int stepNumber{1};

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);
        step.first->toggleSteadyState();
        step.first->solve();
        sim.printMassBalances(debug);
        step.first->toggleSteadyState();
        ++stepNumber;
    }
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
