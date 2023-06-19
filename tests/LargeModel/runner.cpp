#include "runner.hpp"

namespace GlobalFlow {

    using namespace std;

    //int sim_id{0};

    void NZRunner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_nz.json");
    }

    void NZRunner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            sim.getNodes()->at(j)->setSimpleK();
        }
        LOG(debug) << "simple k set for all nodes" << std::endl;
        _eq = sim.getEquation();
    }

    void NZRunner::simulate() {
        Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::MONTH, 1);
        LOG(debug) << "NodeID 1: " << sim.getNodes()->at(1);

        // For node infos:
        ofstream myfile;
        myfile.open ("node_attributes_output.csv");
        myfile << "nodeID,spatID,lon,lat,hyd_cond,hasGHB,river_flow,river_conduct,river_bottom,elevation,initial_head,river_head" << std::endl;
        // spatID,area,neighbour_count,lake,global_lake,wetland,global_wetland,recharge
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            myfile <<
                   sim.getNodes()->at(j)->getID() << "," <<
                   setprecision(7) << sim.getNodes()->at(j)->getSpatID() << setprecision(default_precision) << "," <<
                   sim.getNodes()->at(j)->getLon() << "," <<
                   sim.getNodes()->at(j)->getLat() << "," <<
                   //sim.getNodes()->at(j)->getArea().value() << "," <<
                   //sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                   sim.getNodes()->at(j)->getK().value() << "," <<
                   sim.getNodes()->at(j)->hasGHB() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::LAKE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_LAKE).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::WETLAND).value() << "," <<
                   //sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_WETLAND).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getConductance().value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getBottom().value() << "," <<
                   sim.getNodes()->at(j)->getElevation().value() << "," <<
                   sim.getNodes()->at(j)->getHead().value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowByName(Model::RIVER_MM).getFlowHead().value() <<
                   std::endl;
        }
        myfile.close();

        for (Simulation::step step : stepper) {
            //step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved step with " << step.first->getItter() << " iterations and error of: " << step.first->getError() << std::endl;
            sim.printMassBalances(debug);
            //step.first->toggleSteadyState();
        }
        sim.save();
    }

    void NZRunner::writeData() {
	DataProcessing::DataOutput::OutputManager("data/out_nz.json", sim).write();
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
