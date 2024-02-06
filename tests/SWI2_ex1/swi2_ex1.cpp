#include "swi2_ex1.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_swi2_ex1.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SWI2_ex1_DataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}

void StandaloneRunner::writeNodeInfosToCSV(){
        std::ofstream myfile("swi2_ex1_node_attributes.csv");
        myfile << "nodeID,lon,lat,neig_count,neigs,EL,bottom,K,Por_eff,zeta1,GWR, C_GHB, EL_GHB" << std::endl;

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }
            myfile << sim.getNodes()->at(j)->getID() << "," <<
                   sim.getNodes()->at(j)->getLon() << "," <<
                   sim.getNodes()->at(j)->getLat() << "," <<
                   sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                   neighboursStr << "," <<
                   sim.getNodes()->at(j)->getElevation().value() << "," <<
                   sim.getNodes()->at(j)->getBottom().value() << "," <<
                   sim.getNodes()->at(j)->getK().value() << "," <<
                   sim.getNodes()->at(j)->getEffectivePorosity().value() << "," <<
                   sim.getNodes()->at(j)->getZeta(1).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowConductance(Model::GENERAL_HEAD_BOUNDARY) << "," <<
                   sim.getNodes()->at(j)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY) <<
                   std::endl;
        }
        myfile.close();
    }

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_DAYS, 200);
    int stepNumber = 1;

    // for saving zetas in a csv
    std::ofstream myfile("swi2_ex1_timesteps.csv");
    myfile << "timestep,nodeID,zeta1,head" << std::endl;

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " << stepNumber;
        step.first->toggleSteadyState();
        step.first->solve();
        sim.printMassBalances(debug);

        // for saving zetas in a csv
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            myfile << stepNumber
            << "," << sim.getNodes()->at(j)->getID()
            << "," << sim.getNodes()->at(j)->getZeta(1).value()
            << "," << sim.getNodes()->at(j)->getHead().value()
            << std::endl;
        }
        step.first->toggleSteadyState();
        stepNumber++;
    }
    myfile.close(); // for saving zetas in a csv
}

void StandaloneRunner::getResults() {}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_swi2_ex1.json", sim).write();
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
