/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef GLOBAL_FLOW_SIMULATION_H
#define GLOBAL_FLOW_SIMULATION_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <type_traits>
#include <boost/filesystem.hpp>

#include "../Solver/Equation.hpp"
#include "../DataProcessing/DataReader.hpp"
#include "../Logging/Logging.hpp"
#include "../Model/Node.hpp"
#include "../DataProcessing/Neighbouring.hpp"
#include "../Misc/Helpers.hpp"
//#include "sensitivity.hpp"


namespace GlobalFlow {
    namespace Simulation {
        using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<Model::NodeInterface>>>;
        using eq_ptr = std::unique_ptr<GlobalFlow::Solver::Equation>;
        using namespace boost::multiprecision;
        using namespace boost::units;
        namespace fs = boost::filesystem;

        /**
         * @class MassError
         * Simple container for the mass error calulations
         */
        class MassError {
        public:
            MassError(mpf_float_1000 OUT, mpf_float_1000 IN, mpf_float ERR) : OUT(OUT), IN(IN), ERR(ERR) {}

            mpf_float_1000 OUT;
            mpf_float_1000 IN;
            mpf_float ERR = 0;
        };

        /**
         * @class Simulation
         * The simulation class which holds the equation, options and data instance
         * Further contains methods for calulating the mass balance and sensitivity methods
         */
        class Simulation {
            eq_ptr eq;
            DataReader *reader;
            Options op;

        public:
            Simulation() {}
            Simulation(Options op, DataReader *reader); //impl in cpp - not pretty but works

            Solver::Equation *getEquation() { return eq.get(); };

            void save() {
                //Serialize current node state
                if (serialize) {
                    cout << "Saving nodes for faster reboot .. \n";
                    {
                        std::ofstream ofs("savedNodes", ios::out | ios::binary);
                        boost::archive::binary_oarchive outStream(ofs);
                        // write class instance to archive
                        outStream << nodes;
                    }
                    cout << "Nodes saved\n";
                }
            };

            /**
             * Get basic node information by its id
             * @param nodeID
             * @return A string of information
             */
            std::string NodeInfosByID(unsigned long nodeID) {
                Model::NodeInterface *nodeInterface = nodes->at(nodeID).get();
                string out("\n");
                out += "IN: ";
                out += to_string(nodeInterface->getCurrentIN().value());
                out += "\nOUT: ";
                out += to_string(nodeInterface->getCurrentOUT().value());
                out += "\nElevation: ";
                out +=
                        to_string(nodeInterface->getProperties().get<quantity<Model::Meter>,
                                Model::Elevation>().value());
                out += "\nRHS: ";
                out += to_string(nodeInterface->getRHS().value());
                out += "\nHEAD: ";
                out += to_string(nodeInterface->getProperties().get<quantity<Model::Meter>, Model::Head>().value());
                out += "\nArea: ";
                out +=
                        to_string(nodeInterface->getProperties().get<quantity<Model::SquareMeter>,
                                Model::Area>().value());
                out += "\nStorageFlow: ";
                out += to_string(nodeInterface->getTotalStorageFlow().value());
                out += "\nNONStorageFlowIN: ";
                out += to_string(nodeInterface->getNonStorageFlow([](double a) -> bool {
                    return a > 0;
                }).value());
                out += "\nNONStorageFlowOUT: ";
                out += to_string(nodeInterface->getNonStorageFlow([](double a) -> bool {
                    return a < 0;
                }).value());
                out += "\n";
                std::ostringstream strs;
                strs << nodeInterface;
                out += strs.str();

                return out;
            }

            /**
             * Get budget per node
             * @param ids
             * @return
             */
            template<int FieldNum>
            std::string getFlowSumByIDs(std::array<int, FieldNum> ids) {
                double in{0};
                double out{0};
                for (int j = 0; j < nodes->size(); ++j) {
                    auto id = std::find(std::begin(ids), std::end(ids), nodes->at(j)->getID());
                    if (id != std::end(ids)) {
                        in += nodes->at(j)->getIN().value();
                        out += nodes->at(j)->getOUT().value();
                    }
                }
                ostringstream output;
                output << to_string(in) << "\n";
                output << to_string(out);
                return output.str();
            }

            /**
             * Return all external flows seperatly
             */
            std::string NodeFlowsByID(unsigned long nodeID) {
                long id{0};
                for (int j = 0; j < nodes->size(); ++j) {
                    if (nodes->at(j)->getID() == nodeID) {
                        id = j;
                    }
                }
                Model::NodeInterface *nodeInterface = nodes->at(id).get();
                string out("");
                //Flows Budget for
                //ID, Elevation, Head, IN, OUT, Recharge, River_MM, Lake, Wetland
                out += to_string(nodeID);
                out += ",";
                out += to_string(nodeInterface->getProperties().get<quantity<Model::Meter>,
                        Model::Elevation>().value());
                out += ",";
                out += to_string(nodeInterface->getProperties().get<quantity<Model::Meter>,
                        Model::Head>().value());
                out += ",";
                out += to_string(nodeInterface->getIN().value());
                out += ",";
                out += to_string(nodeInterface->getOUT().value());
                out += ",";
                out += to_string(nodeInterface->getExternalFlowVolumeByName(Model::RECHARGE).value());
                out += ",";
                out += to_string(nodeInterface->getExternalFlowVolumeByName(Model::RIVER_MM).value());
                out += ",";
                out += to_string(nodeInterface->getExternalFlowVolumeByName(Model::LAKE).value());
                out += ",";
                out += to_string(nodeInterface->getExternalFlowVolumeByName(Model::WETLAND).value());
                return out;
            }

            /**
             * Calulate the mass error
             * @param fun1 Function to get OutFlow
             * @param fun2 Function to get InFlow
             * @return
             */
            template<class FunOut, class FunIn>
            MassError getError(FunOut fun1, FunIn fun2) {
                mpf_float_1000 out = 0;
                mpf_float_1000 in = 0;
                mpf_float error = 0;
                for (int j = 0; j < nodes->size(); ++j) {
                    out = out + fun1(j);
                    in = in + fun2(j);
                }
                if (abs(in - out) > 0.00001) {
                    error = ((100 * (in - out)) / ((in + out) / 2));
                }
                MassError err(out, in, error);
                return err;
            }

            /**
             * Get the total mass balance
             * @return
             */
            MassError getMassError() {
                return getError([this](int pos) {
                                    return nodes->at(pos)->getOUT().value();
                                },
                                [this](int pos) {
                                    return nodes->at(pos)->getIN().value();
                                });
            }

            /**
             * Get the mass balance for the current step
             * @return
             */
            MassError getCurrentMassError() {
                return getError([this](int pos) {
                                    return nodes->at(pos)->getCurrentOUT().value();
                                },
                                [this](int pos) {
                                    return nodes->at(pos)->getCurrentIN().value();
                                });
            }

            /**
             * Get the flow lost to external flows
             * @return
             */
            double getLossToRivers() {
                double out = 0.0;
                for (int j = 0; j < nodes->size(); ++j) {
                    out += nodes->at(j)->getNonStorageFlow([](double a) -> bool {
                        return a < 0;
                    }).value();
                }
                return out;
            }

            enum Flows {
                RIVERS = 1,
                DRAINS,
                RIVER_MM,
                LAKES,
                WETLANDS,
                GLOBAL_WETLANDS,
                RECHARGE,
                FASTSURFACE,
                NAG,
                STORAGE,
                GENERAL_HEAD_BOUNDARY
            };

            /**
             * Decide if its an In or Outflow
             * @param fun
             * @return
             */
            template<class Fun>
            MassError inline getError(Fun fun) {
                return getError([this, fun](int pos) {
                                    double flow = fun(pos);
                                    return flow < 0 ? flow : 0;
                                },
                                [this, fun](int pos) {
                                    double flow = fun(pos);
                                    return flow > 0 ? flow : 0;
                                });
            }

            /**
             * Helper function for printing the mass balance for each flow
             * @param flow
             * @return
             */
            string getFlowByName(Flows flow) {
                ostringstream stream;
                MassError tmp(0, 0, 0);
                switch (flow) {
                    case GENERAL_HEAD_BOUNDARY:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::GENERAL_HEAD_BOUNDARY).value();
                        });
                        break;
                    case RECHARGE:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::RECHARGE).value();
                        });
                        break;
                    case FASTSURFACE:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::FAST_SURFACE_RUNOFF).value();
                        });
                        break;
                    case NAG:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::NET_ABSTRACTION).value();
                        });
                        break;
                    case RIVERS:
                        tmp = getError(
                                [this](int i) {
                                    return nodes->at(i)->getExternalFlowVolumeByName(Model::RIVER).value();
                                });
                        break;
                    case DRAINS:
                        tmp = getError(
                                [this](int i) {
                                    return nodes->at(i)->getExternalFlowVolumeByName(Model::DRAIN).value();
                                });
                        break;
                    case RIVER_MM:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::RIVER_MM).value();
                        });
                        break;
                    case LAKES:
                        tmp = getError(
                                [this](int i) {
                                    return nodes->at(i)->getExternalFlowVolumeByName(Model::LAKE).value();
                                });
                        break;
                    case WETLANDS:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::WETLAND).value();
                        });
                        break;
                    case GLOBAL_WETLANDS:
                        tmp = getError([this](int i) {
                            return nodes->at(i)->getExternalFlowVolumeByName(Model::GLOBAL_WETLAND).value();
                        });
                        break;
                    case STORAGE:
                        tmp = getError([this](int i) { return nodes->at(i)->getTotalStorageFlow().value(); });
                        break;
                }

                stream << "IN :" << tmp.IN << "  OUT :" << tmp.OUT;
                return stream.str();
            }

            /**
             * Prints all mass balances
             */
            void printMassBalances() {
                MassError currentErr = getCurrentMassError();
                MassError totalErr = getMassError();
                LOG(stateinfo) << "All units in meter per stepsize";
                LOG(stateinfo) << "Step mass error: " << currentErr.ERR << "  IN: " << currentErr.IN << "  Out: "
                              << currentErr.OUT;
                LOG(stateinfo) << "Total mass error: " << totalErr.ERR << "IN: " << totalErr.IN << "Out: "
                              << totalErr.OUT;
                LOG(stateinfo) << "General Head Boundary" << getFlowByName(GENERAL_HEAD_BOUNDARY);
                LOG(stateinfo) << "Rivers: " << getFlowByName(RIVERS);
                //LOG(stateinfo) << "Drains: " << getFlowByName(DRAINS);
                LOG(stateinfo) << "Dynamic Rivers: " << getFlowByName(RIVER_MM);
                LOG(stateinfo) << "Lakes: " << getFlowByName(LAKES);
                LOG(stateinfo) << "Wetlands: " << getFlowByName(WETLANDS);
                LOG(stateinfo) << "Global wetlands: " << getFlowByName(GLOBAL_WETLANDS);
                LOG(stateinfo) << "Recharge: " << getFlowByName(RECHARGE);
                //LOG(userinfo) << "Fast Surface Runoff: " << getFlowByName(FASTSURFACE) << "\n";
                LOG(stateinfo) << "Net abstraction from groudnwater: " << getFlowByName(NAG);
                LOG(stateinfo) << "Storage (only valid if transient run): " << getFlowByName(STORAGE);
            }

            DataReader *getDataReader() {
                return this->reader;
            }

            NodeVector const &getNodes(){
                return nodes;
            }

            /**
             * Get the residuals of the current iteration
             * @param path
             */
            void writeResiduals(string path) {
                Eigen::VectorXd vector = eq->getResiduals();
                std::ofstream ofs;
                ofs.open(path, std::ofstream::out | std::ofstream::trunc);
                ofs << "X,Y,data" << "\n";
                for (int i = 0; i < vector.size(); ++i) {
                    //i is node position ->get X and Y
                    auto lat = nodes->at(i)->getProperties().get<double, Model::Lat>();
                    auto lon = nodes->at(i)->getProperties().get<double, Model::Lon>();
                    double val = vector[i];
                    ofs << lat << "," << lon << "," << val << "\n";
                }
                ofs.close();
            }


        private:
            NodeVector nodes;

            int initNodes() {
                LOG(userinfo) << FMAG(BOLD("Starting G³M"));
                nodes->reserve(op.getNumberOfNodes() * op.getNumberOfLayers());
                reader->initNodes(nodes);
                reader->readData(op);
                return op.getNumberOfNodes() * op.getNumberOfLayers();
            };


            bool serialize = false;
            bool loadNodes = false;
        };

    }
}//ns
#endif
