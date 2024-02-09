#include "runner.hpp"

namespace GlobalFlow {

    void Runner::loadSettings() {
        pathToConfig = "data/config_na.json"; // nodes per layer: grid_na_dk: 452736, filtered: 381205
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
                step.first->setSteadyState();
                step.first->solve();
                LOG(userinfo) << "Solved step " << stepNumber << " (steady state) with " << step.first->getItter()
                              << " iteration(s)";
                sim.printMassBalances(debug);
                sim.saveStepResults(pathToOutput, stepNumber);
                step.first->setTransient();
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
        myfile << "nodeID,spatID,layer,lon,lat,area,neighbour_count,K,headActive,hasGHB,C$_{GHB}$,EL$_{GHB}$,Por$_{eff}$,EL,GWR,"
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
                   << "," << sim.getNodes()->at(j)->getHeadActive()
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
