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
    std::ofstream myfile("swi2_ex3_node_attributes.csv");
    myfile << "nodeID,spatID,lon,lat,area,edgeLR,edgeFB,neig_count,neighbours,EL,bottom,K,hasGHB,C_GHB,EL_GHB,Por_eff,zeta1_ini,GWR" << std::endl;

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        std::string neighbours;
        for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
            neighbours += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
        }
        myfile << j << "," <<
               sim.getNodes()->at(j)->getSpatID() << "," <<
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
               sim.getNodes()->at(j)->getExternalFlowConductance(Model::GENERAL_HEAD_BOUNDARY) << "," <<
               sim.getNodes()->at(j)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY) << "," <<
               sim.getNodes()->at(j)->getEffectivePorosity().value() << "," <<
               sim.getNodes()->at(j)->getZeta(1).value() << "," <<
               sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value() <<
               std::endl;
    }
    myfile.close();
}

void StandaloneRunner::simulate() {
    std::vector<bool> isSteadyState = op.getStressPeriodSteadyState();
    std::vector<int> numberOfSteps = op.getStressPeriodSteps();
    std::vector<std::string> stepSizes = op.getStressPeriodStepSizes();
    std::vector<bool> isDensityVariable = op.getStressPeriodVariableDensity();

    int stepNumber{1};

    std::ofstream myfile("swi2_ex3_timesteps.csv");
    myfile << "timestep,nodeID,zeta1,head" << std::endl;

    for (int strssPrd = 0; strssPrd < isSteadyState.size(); ++strssPrd) {
        LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd] << " step(s), with stepsize " <<
                      stepSizes[strssPrd];
        // set zetas if previous stress period had no variable density simulation
        if (strssPrd == 1) {
            // Updating recharge at nodes 199 and 599
            std::unordered_map<large_num, double> nodeID_to_newRecharge = {{199, 0.00025},{599, 0.0005}}; // new recharge (per unit area [m^2])
            for (const auto &[nodeID, newRecharge]: nodeID_to_newRecharge ) {
                sim.getNodes()->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                            0 * Model::si::meter,
                                                            newRecharge * sim.getNodes()->at(nodeID)->getArea().value(),
                                                            0 * Model::si::meter);
            }
        }

        Simulation::Stepper stepper = Simulation::Stepper(_eq, stepSizes[strssPrd], isSteadyState[strssPrd],
                                                          isDensityVariable[strssPrd], numberOfSteps[strssPrd]);
        for (Simulation::step step : stepper) {
            step.first->solve();
            LOG(userinfo) << "Step " << stepNumber << ": ";
            LOG(userinfo) << " - Groundwater flow solved with " << step.first->getItter() << " iteration(s)";
            if (isDensityVariable[strssPrd]) {
                LOG(userinfo) << " - Variable density solved with " << step.first->getItter_zetas() << " iteration(s)";
            }
            sim.printMassBalances(debug, isDensityVariable[strssPrd]);
            // for saving zetas in a csv
            for (int j = 0; j < sim.getNodes()->size(); ++j) {
                myfile << stepNumber
                       << "," << sim.getNodes()->at(j)->getID()
                       << "," << sim.getNodes()->at(j)->getZeta(1).value()
                       << "," << sim.getNodes()->at(j)->getHead().value()
                       << std::endl;
            }
            ++stepNumber;
        }
    }
    myfile.close();
    sim.saveNodeState();
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
