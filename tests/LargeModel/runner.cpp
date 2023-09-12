#include "runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_two_layers.json");
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::MONTH, 10);
        //LOG(debug) << "NodeID 1: " << sim.getNodes()->at(1);

        // for saving zetas in a csv
        std::ofstream myfile;
        myfile.open ("timestep_results.csv");
        //myfile << "timestep,spatID,lon,lat,head,zeta1,zeta2,zeta3,zeta4" << std::endl;
        myfile << "timestep,spatID,lon,lat,head" << std::endl;

        int stepNumber{1};
        for (Simulation::step step : stepper) {
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved step " << stepNumber << " with " << step.first->getItter() << " iterations"; // and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);

            // for saving zetas in a csv
            for (int j = 0; j < sim.getNodes()->size(); ++j) {
                myfile << stepNumber << "," <<
                          sim.getNodes()->at(j)->getSpatID() << "," <<
                          sim.getNodes()->at(j)->getLon() << "," <<
                          sim.getNodes()->at(j)->getLat() << "," <<
                          sim.getNodes()->at(j)->getHead().value() <<
                          //sim.getNodes()->at(j)->getZetaIfActive(1).value() << "," <<
                          //sim.getNodes()->at(j)->getZetaIfActive(2).value() << "," <<
                          //sim.getNodes()->at(j)->getZetaIfActive(3).value() << "," <<
                          //sim.getNodes()->at(j)->getZetaIfActive(4).value() <<
                          std::endl;
            }

            step.first->toggleSteadyState();
            ++stepNumber;
        }
        sim.save();
        myfile.close(); // for saving zetas in a csv
    }

    void Runner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out.json", sim).write();
    }

    void Runner::writeNodeInfosToCSV(){
        // For node infos:
        std::ofstream myfile;
        myfile.open ("node_attributes_output.csv");
        myfile << "nodeID,spatID,lon,lat,area,neighbour_count,hyd_cond,hasGHB,effPor,recharge,lake,wetland,river,river_cond,river_bottom,river_head,elevation,initial_head" << std::endl;
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }
            myfile <<
                   sim.getNodes()->at(j)->getID() << "," <<
                   std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision) << "," <<
                   sim.getNodes()->at(j)->getLon() << "," <<
                   sim.getNodes()->at(j)->getLat() << "," <<
                   sim.getNodes()->at(j)->getArea().value() << "," <<
                   sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                   //neighboursStr << "," <<
                   sim.getNodes()->at(j)->getK().value() << "," <<
                   sim.getNodes()->at(j)->hasGHB() << "," <<
                   sim.getNodes()->at(j)->getEffectivePorosity() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::LAKE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_LAKE).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::WETLAND).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_WETLAND).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getConductance().value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getBottom().value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getFlowHead().value() << "," <<
                   sim.getNodes()->at(j)->getElevation().value() << "," <<
                   sim.getNodes()->at(j)->getHead().value() <<
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
