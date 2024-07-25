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
    std::vector<bool> isSteadyState = op.getStressPeriodSteadyState();
    std::vector<int> numberOfSteps = op.getStressPeriodSteps();
    std::vector<std::string> stepSizes = op.getStressPeriodStepSizes();
    std::vector<bool> isDensityVariable = op.getStressPeriodVariableDensity();

    int stepNumber{1};

    for (int strssPrd = 0; strssPrd < isSteadyState.size(); ++strssPrd) {
        LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd] << " step(s), with stepsize " <<
                      stepSizes[strssPrd];

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
            ++stepNumber;
        }
    }
    sim.saveNodeState();
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
