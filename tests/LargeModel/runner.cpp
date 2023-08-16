#include "runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config.json");
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::MONTH, 1);
        //LOG(debug) << "NodeID 1: " << sim.getNodes()->at(1);

        int stepNumber{1};
        for (Simulation::step step : stepper) {
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Step " << stepNumber << ":\n";
            LOG(userinfo) << "- Solved with " << step.first->getItter() << " iterations and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);
            step.first->toggleSteadyState();
            ++stepNumber;
        }
        sim.save();
    }

    void Runner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out.json", sim).write();
    }

    void Runner::writeNodeInfosToCSV(){
        // For node infos:
        std::ofstream myfile;
        myfile.open ("node_attributes_output.csv");
        myfile << "nodeID,spatID,lon,lat,neighbour_count,neighbours,hyd_cond,hasGHB,elevation,initial_head" << std::endl;
        // spatID,area,neighbour_count,lake,global_lake,wetland,global_wetland,recharge
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            std::string neighbours{""};
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighbours += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }
            myfile <<
                   sim.getNodes()->at(j)->getID() << "," <<
                   std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision) << "," <<
                   sim.getNodes()->at(j)->getLon() << "," <<
                   sim.getNodes()->at(j)->getLat() << "," <<
                   //sim.getNodes()->at(j)->getArea().value() << "," <<
                   sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                   neighbours << "," <<
                   sim.getNodes()->at(j)->getK().value() << "," <<
                   sim.getNodes()->at(j)->hasGHB() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::LAKE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_LAKE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::WETLAND).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_WETLAND).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getConductance().value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getBottom().value() << "," <<
                   sim.getNodes()->at(j)->getElevation().value() << "," <<
                   sim.getNodes()->at(j)->getHead().value() <<
                   //sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getFlowHead().value() << "," <<
                   //sim.getNodes()->at(j)->getZeta(1).value() <<
                   std::endl;

        }
        myfile.close();
    }

    Runner::Runner() {}

}//ns



int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::Runner runner;
    std::cout << FBLU("- loading settings...") << std::endl;
    runner.loadSettings();
    std::cout << FBLU("- simulation setup...") << std::endl;
    runner.setupSimulation();
    runner.writeNodeInfosToCSV();
    std::cout << FBLU("- simulating...") << std::endl;
    runner.simulate();
    std::cout << FBLU("- writing output...") << std::endl;
    runner.writeData();
    return EXIT_SUCCESS;
}
