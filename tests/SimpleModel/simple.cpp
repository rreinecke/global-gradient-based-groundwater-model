#include "simple.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simple.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleDataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}

void StandaloneRunner::writeNodeInfosToCSV(){
    std::ofstream myfile;
    myfile.open ("node_attributes_simple.csv");
    myfile << "nodeID,spatID,refID,lon,lat,neighbour_count,neighbours,elevation,bottom,hyd_cond,hasGHB,recharge" << std::endl;

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        std::string neighbours;
        for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
            neighbours += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
        }
        myfile << j << "," <<
               sim.getNodes()->at(j)->getSpatID() << "," <<
               sim.getNodes()->at(j)->getRefID() << "," <<
               sim.getNodes()->at(j)->getLon() << "," <<
               sim.getNodes()->at(j)->getLat() << "," <<
               sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
               neighbours << "," <<
               sim.getNodes()->at(j)->getElevation().value() << "," <<
               sim.getNodes()->at(j)->getBottom().value() << "," <<
               sim.getNodes()->at(j)->getK().value() << "," <<
               sim.getNodes()->at(j)->hasGHB() << "," <<
               sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() <<
               std::endl;
    }
    myfile.close();
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

    int stepNumber = 1;
    Simulation::Stepper transientStepper = Simulation::Stepper(_eq, Simulation::DAY, 10);
    for (Simulation::step step : transientStepper) {
        LOG(userinfo) << "Running transient step " << stepNumber;
        step.first->solve();
        sim.printMassBalances(debug);
        stepNumber++;
    }

    //Changing stresses
    std::vector<int> ids = {90,
                            91,
                            92,
                            93,
                            94,
                            95,
                            96,
                            97,
                            98,
                            99};
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        if(std::find(ids.begin(), ids.end(), j) != ids.end())
            sim.getNodes()->at(j)->updateUniqueFlow(0.5, Model::RECHARGE, false);
    }
  
    LOG(userinfo) << "Running transient steps with changed stresses";
    Simulation::Stepper transientStepper2 = Simulation::Stepper(_eq, Simulation::DAY, 10);
    for (Simulation::step step : transientStepper2) {
        LOG(userinfo) << "Running transient step " << stepNumber;
        step.first->solve();
        sim.printMassBalances(debug);
        stepNumber++;
    }

    delete reader;
}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simple.json", sim).write();
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