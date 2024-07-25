#include "simpleVDF.hpp"


namespace GlobalFlow {

void StandaloneRunner::loadSettings() {
    op = Simulation::Options();
    op.load("data/config_simpleVDF.json");
}

void StandaloneRunner::setupSimulation() {
    reader = new DataProcessing::SimpleVDFDataReader();
    sim = Simulation::Simulation(op, reader);
    _eq = sim.getEquation();
}


void StandaloneRunner::writeNodeInfosToCSV() {
    // For node infos:
    std::ofstream myfile;
    myfile.open ("node_attributes_simpleVDF.csv");
    myfile << "nodeID,lon,lat,neighbour_count,elevation,bottom,hyd_cond,zeta[1]" << std::endl;

    for (int j = 0; j < sim.getNodes()->size(); ++j) {
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

void StandaloneRunner::getResults() {}

void StandaloneRunner::writeData() {
    DataProcessing::DataOutput::OutputManager("data/out_simpleVDF.json", sim).write();
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
