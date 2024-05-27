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
        std::vector<std::string> variablesToSave = {"head", "zeta1"};

        for (int strssPrd = 0; strssPrd < isSteadyState.size(); ++strssPrd) {
            LOG(userinfo) << "Stress period " << strssPrd+1 << ": " << numberOfSteps[strssPrd]
                          << " step(s), with stepsize " << stepSizes[strssPrd];
            // set zetas if previous stress period had no variable density simulation
            if (strssPrd > 0) {
                if (isDensityVariable[strssPrd] and !isDensityVariable[strssPrd-1]) {
                    reader->setInitialZetas(1);
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
        myfile << "nodeID,spatID,layer,lon,lat,front,back,left,right,area,neighbour_count,K,hasGHB,C$_{GHB}$,EL$_{GHB}$,Por$_{eff}$,EL,GWR,"
               << "C$_{river}$,EL$_{river}$,H$_{ini}$" << std::endl;
        int front;
        int back;
        int left;
        int right;

        for (auto & node : *sim.getNodes()) {
            try { front = node->getNeighbour(Model::FRONT)->getSpatID(); } catch (...) { front = -1;}
            try { back = node->getNeighbour(Model::BACK)->getSpatID(); } catch (...) { back = -1;}
            try { left = node->getNeighbour(Model::LEFT)->getSpatID(); } catch (...) { left = -1;}
            try { right = node->getNeighbour(Model::RIGHT)->getSpatID(); } catch (...) { right = -1;}

            const auto default_precision = (int) std::cout.precision();
            myfile << node->getID()
                   << "," << std::setprecision(7) << node->getSpatID() << std::setprecision(default_precision)
                   << "," << node->getLayer()
                   << "," << node->getLon()
                   << "," << node->getLat()
                   << "," << front << "," << back << "," << left << "," << right
                   << "," << node->getArea().value()
                   << "," << node->getListOfNeighbours().size()
                   << "," << node->getK().value()
                   << "," << node->hasGHB()
                   << "," << node->getExternalFlowConductance(Model::GENERAL_HEAD_BOUNDARY)
                   << "," << node->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY)
                   << "," << node->getEffectivePorosity()
                   << "," << node->getElevation().value()
                   << "," << node->getExternalFlowVolumeByName(Model::RECHARGE).value()
                   << "," << node->getExternalFlowConductance(Model::RIVER_MM)
                   << "," << node->getExternalFlowElevation(Model::RIVER_MM)
                   << "," << node->getHead().value()
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
