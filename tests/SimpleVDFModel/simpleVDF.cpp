#include "simpleVDF.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleVDF.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleVDFDataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}


void StandaloneRunner::writeNodeInfosToCSV() {
    // For node infos:
    std::ofstream myfile;
    myfile.open ("node_attributes_simpleVDF.csv");
    myfile << "nodeID,lon,lat,neighbour_count,elevation,bottom,hyd_cond,zeta[1]" << std::endl;

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        myfile << j << "," <<
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
    LOG(userinfo) << "Running stress period 1";
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::YEAR, 100);
    int stepNumber = 1;
    LOG(debug) << sim.getNodes()->at(0); // printing node properties in debug file

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);
        step.first->toggleSteadyState();
        step.first->solve();
        sim.printMassBalances(debug);
        step.first->toggleSteadyState();
        stepNumber++;
    }

    // saving zetas in a csv
    std::ofstream myfile;
    myfile.open ("zetas.csv");
    myfile << "nodeID,zetaID,lat,lon,zeta" << std::endl;
    int zetaID = 1;
    for (int nodeID = 0; nodeID < sim.getNodes()->size(); ++nodeID) {
        myfile << nodeID << "," <<
               zetaID << "," <<
               sim.getNodes()->at(nodeID)->getLat() << "," <<
               sim.getNodes()->at(nodeID)->getLon() << "," <<
               sim.getNodes()->at(nodeID)->getZeta(zetaID).value() <<
               std::endl;
    }

    LOG(userinfo) << "Running stress period 2";
    Simulation::Stepper stepper2 = Simulation::Stepper(_eq, Simulation::YEAR, 2);
    for (Simulation::step step : stepper2) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);
        step.first->solve();
        sim.printMassBalances(debug);

        stepNumber++;
    }

    //sim.save();
}

void StandaloneRunner::getResults() {}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simpleVDF.json", sim).write();
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
