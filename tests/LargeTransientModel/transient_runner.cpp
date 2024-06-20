#include "transient_runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        pathToConfig = "data/config_transient.json"; // nodes per layer: grid_na: 396787, grid_na_dk: 452736
        op = Simulation::Options();
        op.load(pathToConfig);
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::TransientDataReader();
        sim = Simulation::Simulation(op, reader);
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
        int stepNumber{1};
        boost::gregorian::date date = boost::gregorian::day_clock::universal_day();
        std::stringstream ss;
        ss << date.day() << date.month() << date.year();
        std::string simDate = ss.str();

        std::string pathToOutput = "/mnt/storage/output_transient_" + simDate + "/";
        std::string pathToRecharge = "";
        int year = 1900;
        int month = 12;

        std::vector<bool> isSteadyState = op.getStressPeriodSteadyState();
        std::vector<int> numberOfSteps = op.getStressPeriodSteps();
        std::vector<std::string> stepSizes = op.getStressPeriodStepSizes();
        std::vector<bool> isDensityVariable = op.getStressPeriodVariableDensity();

        std::vector<std::string> variablesToSave = {"head", "zeta1"};

        for (int strssPrd = 0; strssPrd < isSteadyState.size(); ++strssPrd) {
            LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd] << " step(s), with stepsize " <<
                          stepSizes[strssPrd];

            Simulation::Stepper stepper = Simulation::Stepper(_eq, stepSizes[strssPrd], isSteadyState[strssPrd],
                                                              isDensityVariable[strssPrd], numberOfSteps[strssPrd]);
            for (Simulation::step step : stepper) {
                if (month == 12) {
                    year++;
                    month = 1;
                } else {
                    month++;
                }
                pathToRecharge = "../../data/recharge/recharge_WaterGAP_" + std::to_string(year) +
                               "-" + std::to_string(month) + ".csv";
                LOG(debug) << "Reading current data: " << pathToRecharge;
                reader->readGWRecharge(pathToRecharge);
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
        // For node infos:
        std::ofstream myfile("node_attributes_transient.csv");
        myfile << "nodeID,spatID,lon,lat,area,neighbour_count,neighbours,K,hasGHB,ghb,ghb_conductance,ghb_elevation,"
               << "effPor,elevation,Qriver,river_conductance,river_elevation,initial_head" << std::endl;
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }
            myfile << sim.getNodes()->at(j)->getID()
                   << "," << std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision)
                   << "," << sim.getNodes()->at(j)->getLon()
                   << "," << sim.getNodes()->at(j)->getLat()
                   << "," << sim.getNodes()->at(j)->getArea().value()
                   << "," << sim.getNodes()->at(j)->getListOfNeighbours().size()
                   << "," << neighboursStr
                   << "," << sim.getNodes()->at(j)->getK().value()
                   << "," << sim.getNodes()->at(j)->hasGHB()
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GENERAL_HEAD_BOUNDARY).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowConductance(Model::GENERAL_HEAD_BOUNDARY)
                   << "," << sim.getNodes()->at(j)->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY)
                   << "," << sim.getNodes()->at(j)->getEffectivePorosity()
                   << "," << sim.getNodes()->at(j)->getElevation().value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowConductance(Model::RIVER_MM)
                   << "," << sim.getNodes()->at(j)->getExternalFlowElevation(Model::RIVER_MM)
                   << "," << sim.getNodes()->at(j)->getHead().value()
                   << std::endl;
        }
        myfile.close();
    }
    Runner::Runner() = default;
}//ns


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
