#include "swi2_ex3.hpp"

namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_swi2_ex3.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SWI2_ex3_DataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}

void StandaloneRunner::writeNodeInfosToCSV(){
    std::ofstream myfile;
    myfile.open ("swi2_ex3_node_attributes.csv");
    myfile << "nodeID,spatID,refID,lon,lat,area,edgeLR,edgeFB,neighbour_count,neighbours,elevation,bottom,hyd_cond,hasGHB,zeta[1],recharge" << std::endl;

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
               sim.getNodes()->at(j)->getArea().value() << "," <<
               sim.getNodes()->at(j)->getEdgeLengthLeftRight().value() << "," <<
               sim.getNodes()->at(j)->getEdgeLengthFrontBack().value() << "," <<
               sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
               neighbours << "," <<
               sim.getNodes()->at(j)->getElevation().value() << "," <<
               sim.getNodes()->at(j)->getBottom().value() << "," <<
               sim.getNodes()->at(j)->getK().value() << "," <<
               sim.getNodes()->at(j)->hasGHB() << "," <<
               sim.getNodes()->at(j)->getZeta(1).value() << "," <<
               sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() <<
               std::endl;
    }
    myfile.close();
}

void StandaloneRunner::simulate() {
    LOG(userinfo) << "Running stress period 1";
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_YEARS, 500);
    int stepNumber{1};

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);
        step.first->toggleSteadyState(); // turn steady state on
        step.first->solve(); // solve equations
        sim.printMassBalances(debug);
        step.first->toggleSteadyState(); // turn steady state off
        ++stepNumber;
    }

    //Changing recharge at nodes 199 and 599
    std::unordered_map<int, double> recharge_new = {{199, 0.00025},{599, 0.0005}}; // new recharge (per unit area [m^2])
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        if(recharge_new.find(j) != recharge_new.end()) {
            // multiply new recharge with area (here: 20 m^2), then update recharge
            double recharge = recharge_new.at(j) * sim.getNodes()->at(j)->getProperties().get<Model::quantity<Model::SquareMeter>,Model::Area>().value();
            sim.getNodes()->at(j)->updateUniqueFlow(recharge, Model::RECHARGE, false);
        }
    }

    // for saving zetas in a csv
    std::ofstream myfile;
    myfile.open ("swi2_ex3_zeta1_timestep_500.csv");
    myfile << "timestep,nodeID,zeta1" << std::endl;
    // for saving zetas in a csv
    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        double zeta = sim.getNodes()->at(j)->getZeta(1).value();
        myfile << stepNumber << "," << j << "," << zeta << std::endl;
    }
    myfile.close();

    LOG(userinfo) << "Running stress period 2";
    Simulation::Stepper stepper2 = Simulation::Stepper(_eq, Simulation::TWO_YEARS, 500);
    for (Simulation::step step : stepper2) {
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
    DataProcessing::DataOutput::OutputManager("data/out_swi2_ex3.json", sim).write();
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
