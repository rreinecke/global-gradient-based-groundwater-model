#include "runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        pathToConfig = "data/config_na.json";
        op = Simulation::Options();
        op.load(pathToConfig);
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
        Simulation::TimeFrame stepSize = Simulation:: YEAR;
        int stepCount = 2500;
        std::string pathToOutput;

        LOG(userinfo) << "Stepsize is " << stepSize << " day(s)";
        Simulation::Stepper stepper = Simulation::Stepper(_eq, stepSize, 1);
        for (Simulation::step step : stepper) {
            LOG(userinfo) << "Running a steady state step";
            step.first->toggleSteadyState();
            step.first->solve();
            LOG(userinfo) << "Solved steady state step with " << step.first->getItter() << " iteration(s)";
            sim.printMassBalances(debug);
            step.first->toggleSteadyState();
        }

        Simulation::Stepper transientStepper = Simulation::Stepper(_eq, stepSize, stepCount);
        LOG(userinfo) << "Runnning " << stepCount << " transient step(s)";

        int stepNumber{1};

        for (Simulation::step step : transientStepper) {
            if (pathToConfig == "data/config_nz.json") {
                pathToOutput = "";
            } else {
                pathToOutput = "/mnt/storage/";
            }

            step.first->solve();
            LOG(userinfo) << "Solved step " << stepNumber << " with " << step.first->getItter() << " iteration(s)";
            sim.printMassBalances(debug);
            // saving timestep results in CSVs
            std::ofstream zetasFile(pathToOutput + "output/zetas_timestep_" + std::to_string(stepNumber) + "_of_" +
                                    std::to_string(stepCount) + ".csv");
            if (zetasFile.is_open()) {
                zetasFile << "timestep,nodeID,head,zeta0,zeta1active,zeta1,zeta2active,zeta2,zeta3" << std::endl;
                for (int j = 0; j < sim.getNodes()->size(); ++j) {
                    zetasFile << stepNumber
                              << "," << sim.getNodes()->at(j)->getID()
                              << "," << sim.getNodes()->at(j)->getHead().value()
                              << "," << sim.getNodes()->at(j)->getZeta(0).value()
                              << "," << sim.getNodes()->at(j)->isZetaActive(1)
                              << "," << sim.getNodes()->at(j)->getZeta(1).value()
                              << "," << sim.getNodes()->at(j)->isZetaActive(2)
                              << "," << sim.getNodes()->at(j)->getZeta(2).value()
                              << "," << sim.getNodes()->at(j)->getZeta(3).value()
                              << std::endl;
                }
                zetasFile.close();
            }
            std::ofstream flowsFile(pathToOutput + "output/flows_timestep_" + std::to_string(stepNumber) + "_of_" + std::to_string(stepCount) + ".csv");
            if (flowsFile.is_open()) {
                flowsFile << "timestep,nodeID,ghb,front,back,left,right" << std::endl;
                std::unordered_map<Model::NeighbourPosition, double> flowMap;
                for (int j = 0; j < sim.getNodes()->size(); ++j) {
                    flowMap = sim.getNodes()->at(j)->getFlowToOrFromNeighbours();
                    flowsFile << stepNumber
                           << "," << sim.getNodes()->at(j)->getID()
                           << ","
                           << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GENERAL_HEAD_BOUNDARY).value()
                           << "," << flowMap[Model::NeighbourPosition::FRONT]
                           << "," << flowMap[Model::NeighbourPosition::BACK]
                           << "," << flowMap[Model::NeighbourPosition::LEFT]
                           << "," << flowMap[Model::NeighbourPosition::RIGHT]
                           << std::endl;
                }
                flowsFile.close();
            }
            ++stepNumber;
        }
        sim.save();
    }

    void Runner::writeData() {
	    DataProcessing::DataOutput::OutputManager("data/out.json", sim).write();
    }

    void Runner::writeNodeInfosToCSV(){
        // For node infos:
        std::ofstream myfile("node_attributes_large.csv");
        myfile << "nodeID,spatID,lon,lat,area,K,hasGHB,effPor,elevation,recharge,Qriver,initial_head,zeta0,zeta1active,zeta1,zeta2active,zeta2,zeta3" << std::endl; // zeta0,zeta1Active,zeta1,zeta2active,zeta2
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            /*std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }*/
            myfile << sim.getNodes()->at(j)->getID()
                   << "," << std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision)
                   << "," << sim.getNodes()->at(j)->getLon()
                   << "," << sim.getNodes()->at(j)->getLat()
                   << "," << sim.getNodes()->at(j)->getArea().value()
                   //<< "," << sim.getNodes()->at(j)->getListOfNeighbours().size()
                   //<< "," << neighboursStr
                   << "," << sim.getNodes()->at(j)->getK().value()
                   << "," << sim.getNodes()->at(j)->hasGHB()
                   << "," << sim.getNodes()->at(j)->getEffectivePorosity()
                   << "," << sim.getNodes()->at(j)->getElevation().value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value()
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
