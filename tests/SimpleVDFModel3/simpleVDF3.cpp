#include "simpleVDF3.hpp"


namespace GlobalFlow {

using namespace std;

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleVDF3.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleVDF3DataReader();
    sim = Simulation::Simulation(op, reader);


    // For node infos:
    ofstream myfile;
    myfile.open ("node_attributes_output.csv");
    myfile << "nodeID,lon,lat,neighbour_count,elevation,bottom,hyd_cond,zeta[1]" << std::endl;

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        sim.getNodes()->at(j)->setSimpleK();

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
    _eq = sim.getEquation();
}

void StandaloneRunner::simulate() {
    LOG(userinfo) << "Running stress period 1";
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::TWO_YEARS, 500);
    int stepNumber = 1;

    // for saving zetas in a csv
    ofstream myfile;
    myfile.open ("zetas.csv");
    myfile << "timestep,nodeID,zetaID,zeta" << std::endl;

    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);

        if (stepNumber == 1) {
            step.first->toggleSteadyState();
        }
        step.first->solve();
        sim.printMassBalances(debug);

        // for saving zetas in a csv
        int zetaID = 1;
        double zeta;
        for (int nodeID = 0; nodeID < sim.getNodes()->size(); ++nodeID) {
            zeta = sim.getNodes()->at(nodeID)->getZeta(zetaID).value();
            myfile << stepNumber << "," << nodeID << "," << zetaID << "," << zeta << std::endl;
        }

        stepNumber++;
    }

    //Changing stresses
    std::vector<int> ids = {199, 599};

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        if(std::find(ids.begin(), ids.end(), j) != ids.end()) {
            double recharge = 0.5 * sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value();
            sim.getNodes()->at(j)->updateUniqueFlow(recharge, Model::RECHARGE, false);
        }
    }
    LOG(userinfo) << "Running stress period 2";
    Simulation::Stepper stepper2 = Simulation::Stepper(_eq, Simulation::TWO_YEARS, 500);
    for (Simulation::step step : stepper2) {
        LOG(userinfo) << "Running steady state step " + std::to_string(stepNumber);
        step.first->solve();
        sim.printMassBalances(debug);

        stepNumber++;
    }

    myfile.close(); // for saving zetas in a csv
    //sim.save();
}

void StandaloneRunner::getResults() {}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simpleVDF3.json", sim).write();
}

StandaloneRunner::StandaloneRunner() {}

}//ns

int main() {
    std::cout << FBLU("Starting Simulation") << std::endl;
    GlobalFlow::StandaloneRunner runner;
    runner.loadSettings();
    runner.setupSimulation();
    runner.simulate();
    runner.writeData();
    return EXIT_SUCCESS;
}
