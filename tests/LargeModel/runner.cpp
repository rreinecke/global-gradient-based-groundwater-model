#include "runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        op = Simulation::Options();
        op.load("data/config_na.json");
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
        Simulation::TimeFrame stepSize = Simulation:: YEAR;
        int stepCount = 100;

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
        LOG(userinfo) << "Runnning " << stepCount << " transient step(s). Stepsize is " << stepSize << " day(s)";

        // for saving zetas in a csv
        std::ofstream myfile;
        myfile.open ("timestep_results.csv");
        myfile << "timestep,nodeID,lon,lat,head,zeta0,zeta1Active,zeta1,zeta2active,zeta2,zeta3active,zeta3,zeta4" << std::endl;

        int stepNumber{1};
        for (Simulation::step step : transientStepper) {
            step.first->solve();
            LOG(userinfo) << "Solved step " << stepNumber << " with " << step.first->getItter() << " iteration(s)";
            sim.printMassBalances(debug);

            // for saving zetas in a csv
            for (int j = 0; j < sim.getNodes()->size(); ++j) {
                myfile << stepNumber
                       << "," << sim.getNodes()->at(j)->getID()
                       << "," << sim.getNodes()->at(j)->getLon()
                       << "," << sim.getNodes()->at(j)->getLat()
                       << "," << sim.getNodes()->at(j)->getHead().value()
                       << "," << sim.getNodes()->at(j)->getZeta(0).value()
                       << "," << sim.getNodes()->at(j)->isZetaActive(1)
                       << "," << sim.getNodes()->at(j)->getZeta(1).value()
                       << "," << sim.getNodes()->at(j)->isZetaActive(2)
                       << "," << sim.getNodes()->at(j)->getZeta(2).value()
                        << "," << sim.getNodes()->at(j)->isZetaActive(3)
                        << "," << sim.getNodes()->at(j)->getZeta(3).value()
                       << "," << sim.getNodes()->at(j)->getZeta(4).value()
                       << std::endl;
            }
            ++stepNumber;
        }
        sim.save();
        myfile.close(); // for saving zetas in a csv
    }

    void Runner::writeData() {
	    DataProcessing::DataOutput::OutputManager("data/out.json", sim).write();
    }

    void Runner::writeNodeInfosToCSV(){
        // For node infos:
        std::ofstream myfile;
        myfile.open ("node_attributes_large.csv");
        myfile << "nodeID,spatID,numNeig,neig,hasGHB,riverCond,riverBottom,riverElevation,elevation,initial_head" << std::endl; // zeta0,zeta1Active,zeta1,zeta2active,zeta2
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }
            myfile << sim.getNodes()->at(j)->getID()
                   << "," << std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision)
                   //<< "," << sim.getNodes()->at(j)->getLon()
                   //<< "," << sim.getNodes()->at(j)->getLat()
                   //<< "," << sim.getNodes()->at(j)->getArea().value()
                   << "," << sim.getNodes()->at(j)->getListOfNeighbours().size()
                   << "," << neighboursStr
                   //<< "," << sim.getNodes()->at(j)->getK().value()
                   << "," << sim.getNodes()->at(j)->hasGHB()
                   //<< "," << sim.getNodes()->at(j)->getEffectivePorosity()
                   //<< "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RECHARGE).value()
                   //<< "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::LAKE).value()
                   //<< "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_LAKE).value()
                   //<< "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::WETLAND).value()
                   //<< "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::GLOBAL_WETLAND).value()
                   //<< "," << sim.getNodes()->at(j)->getExternalFlowVolumeByName(Model::RIVER_MM).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowConductance(Model::RIVER_MM).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowBottom(Model::RIVER_MM).value()
                   << "," << sim.getNodes()->at(j)->getExternalFlowElevation(Model::RIVER_MM).value()
                   << "," << sim.getNodes()->at(j)->getElevation().value()
                   << "," << sim.getNodes()->at(j)->getHead().value()
                   //<< "," << sim.getNodes()->at(j)->getZeta(0).value()
                   //<< "," << sim.getNodes()->at(j)->isZetaActive(1)
                   //<< "," << sim.getNodes()->at(j)->getZeta(1).value()
                   //<< "," << sim.getNodes()->at(j)->isZetaActive(2)
                   //<< "," << sim.getNodes()->at(j)->getZeta(2).value()
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
