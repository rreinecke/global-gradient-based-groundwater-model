#ifndef TESTING_SWI2_EX1_DATAREADER_HPP
#define TESTING_SWI2_EX1_DATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"
#include "../../src/Model/Node.hpp"
#include "../../src/Model/Units.hpp"

namespace GlobalFlow {
    namespace DataProcessing {

        class SWI2_ex1_DataReader : public DataReader {
        public:
            SWI2_ex1_DataReader() { }

            void readData(Simulation::Options op) override {
                LOG(userinfo) << "Building the initial model layer";
                std::vector<std::vector<int>> grid;
                grid = readGrid(nodes,
                                buildDir(op.getNodesDir()),
                                op.getNumberOfNodesPerLayer(),
                                op.getNumberOfLayers(),
                                op.getNumberOfRows(),
                                op.getNumberOfCols(),
                                op.getInitialK()[0],
                                op.getInitialHead(),
                                op.getAquiferDepth()[0],
                                op.getAnisotropy()[0],
                                op.getSpecificYield(),
                                op.getSpecificStorage(),
                                op.getEdgeLengthLeftRight(),
                                op.getEdgeLengthFrontBack(),
                                op.useEfolding(),
                                op.isConfined(0),
                                op.isDensityVariable(),
                                op.getEffectivePorosity(),
                                op.getMaxTipSlope(),
                                op.getMaxToeSlope(),
                                op.getMinDepthFactor(),
                                op.getSlopeAdjFactor(),
                                op.getVDFLock(),
                                op.getDensityZones());

                LOG(userinfo) << "Building grid by rows and columns (boundaries need to be specified in with a file)";
                DataProcessing::buildByGrid(nodes, grid, op.getNumberOfNodesPerLayer(), op.getNumberOfLayers());

                LOG(userinfo) << "Reading elevation";
                readElevation(buildDir(op.getElevation()));

                if (op.isKFromFile()) {
                    LOG(userinfo) << "Reading hydraulic conductivity";
                    readConduct(buildDir(op.getLithology()));
                }

                if (op.isInitialHeadFromFile()){
                    LOG(userinfo) << "Initializing head";
                    readInitialHeads((buildDir(op.getInitialHeadsDir())));
                }

                LOG(userinfo) << "Reading the groundwater recharge";
                readGWRecharge(buildDir(op.getRecharge()));

                LOG(userinfo) << "Reading the boundary condition";
                readHeadBoundary(buildDir(op.getKGHBDir()));

                if (op.isDensityVariable()) {
                    LOG(userinfo) << "Reading zetas";
                    readInitialZetas(op.getNumberOfNodesPerLayer(), op.getNumberOfLayers(),
                                     buildDir(op.getInitialZetasDir()));
                    if (op.isEffectivePorosityFromFile()){
                        LOG(userinfo) << "Reading effective porosity";
                        readEffectivePorosity(buildDir(op.getEffectivePorosityDir()));
                    }
                    if (op.isZonesSourcesSinksFromFile()){
                        LOG(userinfo) << "Reading zones of sinks and sources";
                        readZonesSourcesSinks(buildDir(op.getZonesOfSourcesAndSinksDir()), op.getDensityZones());
                    }
                }
            }

        private:

            void readGWRecharge(std::string path) override {
                io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
                in.read_header(io::ignore_no_column, "spatID", "layer", "recharge");
                large_num spatID{0};
                int refID{0};
                int layer{0};
                double recharge{0};
                large_num nodeID{0};

                while (in.read_row(spatID, layer, recharge)) {
                    try {
                        nodeID = this->lookupSpatIDtoNodeIDs.at(spatID).at(layer).at(refID);
                    }
                    catch (const std::out_of_range &ex) {
                        //if Node does not exist ignore entry
                        continue;
                    }

                    nodes->at(nodeID)->addExternalFlow(Model::RECHARGE,
                                                       0 * Model::si::meter,
                                                       recharge * nodes->at(nodeID)->getProperties().get<Model::quantity<Model::SquareMeter>,Model::Area>().value(),
                                                       0 * Model::si::meter);
                }
            }
        };
    }
}
#endif //TESTING_SWI2_EX1_DATAREADER_HPP
