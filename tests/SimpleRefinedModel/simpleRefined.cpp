#include "simpleRefined.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleRefined.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleRefinedDataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}

void StandaloneRunner::writeNodeInfosToCSV(){
    std::ofstream myfile;
    myfile.open ("node_attributes_refined.csv");
    myfile << "nodeID,spatID,refID,refinedInto,lon,lat,neighbour_count,neighbours,elevation,bottom,hyd_cond,hasGHB,recharge" << std::endl;

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
        std::string neighbours;
        for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
            neighbours += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
        }
        myfile << j << "," <<
               sim.getNodes()->at(j)->getSpatID() << "," <<
               sim.getNodes()->at(j)->getRefID() << "," <<
               sim.getNodes()->at(j)->getRefinedInto() << "," <<
               sim.getNodes()->at(j)->getLon() << "," <<
               sim.getNodes()->at(j)->getLat() << "," <<
               sim.getNodes()->at(j)->getListOfNeighbours().size() << "," <<
               neighbours << "," <<
               sim.getNodes()->at(j)->getElevation().value() << "," <<
               sim.getNodes()->at(j)->getBottom().value() << "," <<
               sim.getNodes()->at(j)->getK().value() << "," <<
               sim.getNodes()->at(j)->hasGHB() << "," <<
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
        if (strssPrd == 2) {
            //Changing recharge
            for (int j = 0; j < sim.getNodes()->size(); ++j) {
                if (sim.getNodes()->at(j)->hasTypeOfExternalFlow(Model::RECHARGE)){
                    sim.getNodes()->at(j)->updateUniqueFlow(0.5, Model::RECHARGE, false);
                }
            }
        }

        Simulation::Stepper stepper = Simulation::Stepper(_eq, stepSizes[strssPrd], isSteadyState[strssPrd],
                                                          isDensityVariable[strssPrd], numberOfSteps[strssPrd]);
        for (Simulation::step step : stepper) {
            step.first->solve();
            sim.printMassBalances(debug, isDensityVariable[strssPrd]);
            sim.saveStepResults(pathToOutput, stepNumber, variablesToSave, isDensityVariable[strssPrd]);
            LOG(userinfo) << "Step " << stepNumber << ": ";
            LOG(userinfo) << " - Groundwater flow solved with " << step.first->getItter() << " iteration(s)";
            if (isDensityVariable[strssPrd]) {
                LOG(userinfo) << " - Variable density solved with " << step.first->getItter_zetas() << " iteration(s)";
            }
            ++stepNumber;
        }
    }
    sim.saveNodeState();
    delete reader;
}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simpleRefined.json", sim).write();
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