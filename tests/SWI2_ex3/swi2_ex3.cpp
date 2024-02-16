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
    myfile << "nodeID,spatID,refID,lon,lat,area,edgeLR,edgeFB,neig_count,neighbours,EL,bottom,K,hasGHB,C_GHB,EL_GHB,Por_eff,zeta1_ini,GWR" << std::endl;

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
    boost::gregorian::date date = boost::gregorian::day_clock::universal_day();
    std::stringstream ss;
    ss << date.day() << date.month() << date.year();
    std::string simDate = ss.str();
    std::string pathToOutput = "/mnt/storage/output_" + simDate + "/";
    std::vector<std::string> variablesToSave = {"head", "zeta0", "zeta1", "zeta2", "ghb", "sum_neig"};

    for (int strssPrd = 0; strssPrd < isSteadyState.size(); ++strssPrd) {
        LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd] << " step(s), with stepsize " <<
                      stepSizes[strssPrd];
        // set zetas if previous stress period had no variable density simulation
        if (strssPrd == 1) {
            //Changing recharge at nodes 199 and 599
            std::unordered_map<int, double> recharge_new = {{199, 0.00025},{599, 0.0005}}; // new recharge (per unit area [m^2])
            for (int j = 0; j < sim.getNodes()->size(); ++j) {
                if(recharge_new.find(j) != recharge_new.end()) {
                    // multiply new recharge with area (here: 20 m^2), then update recharge
                    double recharge = recharge_new.at(j) * sim.getNodes()->at(j)->getProperties().get<Model::quantity<Model::SquareMeter>,Model::Area>().value();
                    sim.getNodes()->at(j)->updateUniqueFlow(recharge, Model::RECHARGE, false);
                }
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
            sim.saveStepResults(pathToOutput, stepNumber, variablesToSave, isDensityVariable[strssPrd]);
            ++stepNumber;
        }
    }
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
