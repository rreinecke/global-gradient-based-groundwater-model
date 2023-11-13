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
        std::ofstream myfile;
        myfile.open ("node_attributes_vdf1.csv");
        myfile << "nodeID,lon,lat,neighbour_count,elevation,bottom,hyd_cond,zeta[1],recharge" << std::endl;

        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            myfile << j << "," <<
                   sim.getNodes()->at(j)->getLon() << "," <<
                   sim.getNodes()->at(j)->getLat() << "," <<
                   sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
                   sim.getNodes()->at(j)->getElevation().value() << "," <<
                   sim.getNodes()->at(j)->getBottom().value() << "," <<
                   sim.getNodes()->at(j)->getK().value() << "," <<
                   sim.getNodes()->at(j)->getZeta(1).value() << "," <<
                   sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() <<
                   std::endl;
        }
        myfile.close();
    }

void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_DAYS, 200);
    int stepNumber = 1;

    //LOG(debug) << sim.getNodes()->at(1);

    // for saving zetas in a csv
    std::ofstream myfile;
    myfile.open ("swi2_ex1_zeta1_timesteps.csv");
    myfile << "timestep,nodeID,zeta1" << std::endl;

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " << stepNumber;

        if (stepNumber == 1) {
            step.first->toggleSteadyState();
        }
        step.first->solve();
        sim.printMassBalances(debug);

        // for saving zetas in a csv
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            double zeta = sim.getNodes()->at(j)->getZeta(1).value();
            myfile << stepNumber << "," << j << "," << zeta << std::endl;
        }

        stepNumber++;
    }
    myfile.close(); // for saving zetas in a csv
    //sim.save();
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
