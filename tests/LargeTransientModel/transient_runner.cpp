#include "transient_runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        pathToConfig = "data/config_nz.json"; // nodes per layer: grid_na: 396787, grid_na_dk: 452736
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

        std::string pathToOutput;
        if (pathToConfig == "data/config_nz.json") {
            pathToOutput = "output/";
        } else {
            pathToOutput = "/mnt/storage/output_transient_" + simDate + "/";
        }

        std::vector<std::string> rechargeFiles = {"recharge/recharge_WaterGAP_1901-1.csv", // todo try to read netCDF: https://gerasimosmichalitsianos.wordpress.com/2017/12/13/usingcppwithnetcdf/
                                                  "recharge/recharge_WaterGAP_1901-2.csv",
                                                  "recharge/recharge_WaterGAP_1901-3.csv",
                                                  "recharge/recharge_WaterGAP_1901-4.csv",
                                                  "recharge/recharge_WaterGAP_1901-5.csv",
                                                  "recharge/recharge_WaterGAP_1901-6.csv",
                                                  "recharge/recharge_WaterGAP_1901-7.csv",
                                                  "recharge/recharge_WaterGAP_1901-8.csv",
                                                  "recharge/recharge_WaterGAP_1901-9.csv",
                                                  "recharge/recharge_WaterGAP_1901-10.csv",
                                                  "recharge/recharge_WaterGAP_1901-11.csv",
                                                  "recharge/recharge_WaterGAP_1901-12.csv"};

        std::vector<bool> isSteadyState = op.getStressPeriodSteadyState();
        std::vector<int> numberOfSteps = op.getStressPeriodSteps();
        std::vector<std::string> stepSizes = op.getStressPeriodStepSizes();
        std::vector<bool> isDensityVariable = op.getStressPeriodVariableDensity();

        std::vector<std::string> variablesToSave = {"head", "zeta0", "zeta1", "zeta2", "ghb", "sum_neig"};

        for (int strssPrd = 0; strssPrd < isSteadyState.size(); ++strssPrd) {
            LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd] << " step(s), with stepsize " <<
                          stepSizes[strssPrd];
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
        // For node infos:
        std::ofstream myfile("node_attributes_large.csv");
        myfile << "nodeID,spatID,lon,lat,area,neighbour_count,neighbours,K,hasGHB,ghb,ghb_conductance,ghb_elevation,"
               << "effPor,elevation,recharge,Qriver,river_conductance,river_elevation,initial_head,"
               << "zeta0,zeta1active,zeta1,zeta2active,zeta2,zeta3" << std::endl;
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
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowConductance(Model::RIVER_MM)
                   << "," << sim.getNodes()->at(j)->getExternalFlowElevation(Model::RIVER_MM)
                   << "," << sim.getNodes()->at(j)->getHead().value()
                   << "," << sim.getNodes()->at(j)->getZeta(0).value()
                   << "," << sim.getNodes()->at(j)->isZetaActive(1)
                   << "," << sim.getNodes()->at(j)->getZeta(1).value()
                   << "," << sim.getNodes()->at(j)->isZetaActive(2)
                   << "," << sim.getNodes()->at(j)->getZeta(2).value()
                   << "," << sim.getNodes()->at(j)->getZeta(3).value()
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
