#include "runner.hpp"

namespace GlobalFlow {

    //int sim_id{0};

    void Runner::loadSettings() {
        pathToConfig = "data/config_na.json"; // nodes per layer: grid_na: 396787, grid_na_dk: 452736
        op = Simulation::Options();
        op.load(pathToConfig);
    }

    void Runner::setupSimulation() {
        reader = new DataProcessing::GlobalDataReader();
        sim = Simulation::Simulation(op, reader);
        _eq = sim.getEquation();
    }

    void Runner::simulate() {
        std::vector<int> steadyStateStressPeriodSteps = op.getSteadyStateStressPeriodSteps();
        std::vector<int> transientStressPeriodSteps = op.getTransientStressPeriodSteps();
        std::vector<std::string> steadyStateStressPeriodStepsizes = op.getSteadyStateStressPeriodStepsizes();
        std::vector<std::string> transientStressPeriodStepsizes = op.getTransientStressPeriodStepsizes();

        int stepNumber{1};
        boost::gregorian::date date = boost::gregorian::day_clock::universal_day();
        std::stringstream ss;
        ss << date.day() << date.month() << date.year();
        std::string simDate = ss.str();

        std::string pathToOutput = "/mnt/storage/output_" + simDate + "/";
        for (int i = 0; i < steadyStateStressPeriodSteps.size(); ++i) {
            if (steadyStateStressPeriodSteps[i] == 0) { continue; }
            LOG(userinfo) << "Steady state stress period " << i+1 << ": " <<
                          steadyStateStressPeriodSteps[i] << " step(s), with stepsize " <<
                          steadyStateStressPeriodStepsizes[i];
            Simulation::Stepper steadyStepper = Simulation::Stepper(_eq, steadyStateStressPeriodStepsizes[i],
                                                                    steadyStateStressPeriodSteps[i]);
            for (Simulation::step step : steadyStepper) {
                step.first->toggleSteadyState();
                step.first->solve();
                LOG(userinfo) << "Solved step " << stepNumber << " (steady state) with " << step.first->getItter()
                              << " iteration(s)";
                sim.printMassBalances(debug);
                sim.saveStepResults(pathToOutput, stepNumber);
                step.first->toggleSteadyState();
                ++stepNumber;
            }
        }

        for (int i = 0; i < transientStressPeriodSteps.size(); ++i) {
            if (transientStressPeriodSteps[i] == 0) { continue; }
            LOG(userinfo) << "Transient stress period " << i+1 << ": " <<
                          transientStressPeriodSteps[i] << " step(s), with stepsize " <<
                          transientStressPeriodStepsizes[i];
            Simulation::Stepper transientStepper = Simulation::Stepper(_eq, transientStressPeriodStepsizes[i],
                                                                       transientStressPeriodSteps[i]);
            for (Simulation::step step: transientStepper) {
                step.first->solve();
                LOG(userinfo) << "Solved step " << stepNumber << " (transient) with " << step.first->getItter()
                              << " iteration(s)";
                sim.printMassBalances(debug);
                sim.saveStepResults(pathToOutput, stepNumber);
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
        myfile << "nodeID,spatID,lon,lat,area,neighbour_count,K,hasGHB,C_ghb,EL_ghb,"
               << "Por_eff,EL,GWR,C_river,EL_river,H_ini," << std::endl;
               //<< "zeta0,zeta1active,zeta1,zeta2active,zeta2,zeta3" << std::endl; // keeping this for debug
        for (int j = 0; j < sim.getNodes()->size(); ++j) {
            const auto default_precision = (int) std::cout.precision();
            /*std::string neighboursStr;
            for (auto neighbour : sim.getNodes()->at(j)->getListOfNeighbours()){
                neighboursStr += "N:" + std::to_string(neighbour.first) + " ID:" + std::to_string(neighbour.second) + "; ";
            }*/ // keeping this for debug
            myfile << sim.getNodes()->at(j)->getID()
                   << "," << std::setprecision(7) << sim.getNodes()->at(j)->getSpatID() << std::setprecision(default_precision)
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
                   /*<< "," << sim.getNodes()->at(j)->getZeta(0).value()
                   << "," << sim.getNodes()->at(j)->isZetaActive(1)
                   << "," << sim.getNodes()->at(j)->getZeta(1).value()
                   << "," << sim.getNodes()->at(j)->isZetaActive(2)
                   << "," << sim.getNodes()->at(j)->getZeta(2).value()
                   << "," << sim.getNodes()->at(j)->getZeta(3).value()*/ // keeping this for debug
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
