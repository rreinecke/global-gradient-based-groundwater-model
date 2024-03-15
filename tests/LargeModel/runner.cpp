#include "runner.hpp"

namespace GlobalFlow {

    void Runner::loadSettings() {
        pathToConfig = "data/config_na.json"; // nodes per layer: grid_na_dk: 452736
        op = Simulation::Options();
        op.load(pathToConfig);
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader(); // sets data reader
        sim = Simulation::Simulation(op, reader); // calls data reader initializing nodes
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
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
            LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd]
                          << " step(s), with stepsize " << stepSizes[strssPrd];
            // set zetas if previous stress period had no variable density simulation
            if (strssPrd > 0) {
                if (isDensityVariable[strssPrd] and !isDensityVariable[strssPrd-1]) {
                    LOG(userinfo) << "Setting initial zetas using Ghyben-Herzberg";
                    reader->setZetasGhybenHerzberg(op.getNumberOfLayers(), op.getNumberOfNodesPerLayer(), 10, // todo add to config
                                                   op.getDensityZones());
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
    }

    void Runner::writeData() {
	    DataProcessing::DataOutput::OutputManager("data/out.json", sim).write();
    }

    void Runner::writeNodeInfosToCSV(){
        std::ofstream myfile("node_attributes_large.csv");
        myfile << "nodeID,spatID,layer,lon,lat,area,neighbour_count,K,hasGHB,C$_{GHB}$,EL$_{GHB}$,Por$_{eff}$,EL,GWR,"
               << "C$_{river}$,EL$_{river}$,H$_{ini}$" << std::endl;
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            /*std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }*/ // keeping this for debug
            myfile << sim.getNodes()->at(j)->getID()
                   << "," << std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision)
                   << "," << sim.getNodes()->at(j)->getLayer()
                   << "," << sim.getNodes()->at(j)->getLon()
                   << "," << sim.getNodes()->at(j)->getLat()
                   << "," << sim.getNodes()->at(j)->getArea().value()
                   << "," << sim.getNodes()->at(j)->getListOfNeighbours().size()
                   << "," << sim.getNodes()->at(j)->getK().value()
                   << "," << sim.getNodes()->at(j)->hasGHB()
                   << "," << sim.getNodes()->at(j)->getExternalFlowConductance(Model::GENERAL_HEAD_BOUNDARY)
                   << "," << sim.getNodes()->at(j)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY)
                   << "," << sim.getNodes()->at(j)->getEffectivePorosity()
                   << "," << sim.getNodes()->at(j)->getElevation().value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowConductance(Model::RIVER_MM)
                   << "," << sim.getNodes()->at(j)->getExternalFlowElevation(Model::RIVER_MM)
                   << "," << sim.getNodes()->at(j)->getHead().value()
                   << std::endl;
        }
        myfile.close();
    }
    Runner::Runner() = default;
}


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
