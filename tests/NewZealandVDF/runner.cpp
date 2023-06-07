#include "runner.hpp"

namespace GlobalFlow {

    using namespace std;

    //int sim_id{0};

    void NZRunner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_nzVDF.json");
    }

    void NZRunner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);

        LOG(debug) << "NodeID 1: " << sim.getNodes()->at(1);

        // For node infos:
        ofstream myfile;
        myfile.open ("node_attributes_output.csv");
        myfile << "nodeID,spatID,lon,lat,neighbour_count,elevation,hyd_cond,head,hasGHB,recharge,lake,global_lake,wetland,global_wetland,river_mm,river_mm_flowhead" << std::endl;

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            sim.getNodes()->at(j)->setSimpleK();

            myfile <<
                sim.getNodes()->at(j)->getID() << "," <<
                sim.getNodes()->at(j)->getSpatID() << "," <<
                sim.getNodes()->at(j)->getLon() << "," <<
                sim.getNodes()->at(j)->getLat() << "," <<
                sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                sim.getNodes()->at(j)->getElevation().value() << "," <<
                sim.getNodes()->at(j)->getK().value() << "," <<
                sim.getNodes()->at(j)->getHead().value() << "," <<
                sim.getNodes()->at(j)->hasGHB() << "," <<
                sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() << "," <<
                sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::LAKE).value() << "," <<
                sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_LAKE).value() << "," <<
                sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::WETLAND).value() << "," <<
                sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_WETLAND).value() << "," <<
                //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value() << "," <<
                //sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getFlowHead().value() <<
                std::endl;
        }
        LOG(debug) << "simple k set for all nodes" << std::endl;
        myfile.close();
        _eq = sim.getEquation();
    }

    void NZRunner::simulate() {

        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::MONTH, 60);
        for (Simulation::step step : stepper) {
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved step with " << step.first->getItter() << " iterations and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);
            step.first->toggleSteadyState();
        }
        sim.save();
    }

    void NZRunner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out_nzVDF.json", sim).write();
    }

    NZRunner::NZRunner() {}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::NZRunner runner;
    std::cout << FBLU("- loading settings...") << std::endl;
    runner.loadSettings();
    std::cout << FBLU("- simulation setup...") << std::endl;
    runner.setupSimulation();
    std::cout << FBLU("- simulating...") << std::endl;
    runner.simulate();
    std::cout << FBLU("- writing output...") << std::endl;
    runner.writeData();
    return EXIT_SUCCESS;
}
