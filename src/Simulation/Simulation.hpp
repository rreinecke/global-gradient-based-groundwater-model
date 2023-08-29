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
#include <boost/multiprecision/gmp.hpp>
#include <type_traits>
#include <boost/filesystem.hpp>
#include "../Misc/serializer.hpp"
#include "../Solver/Equation.hpp"
#include "../DataProcessing/DataReader.hpp"
#include "../Logging/Logging.hpp"
#include "../Model/Node.hpp"
#include "../DataProcessing/Neighbouring.hpp"
#include "../Misc/Helpers.hpp"
//#include "sensitivity.hpp"

namespace boost { namespace serialization {

        using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

        template<class Archive>
        inline void serialize(Archive & ar, NodeVector& foo, const unsigned int file_version){
            LOG(GlobalFlow::debug) << "Serializing called";
        }

        template<class Archive>
        inline void save_construct_data(Archive & ar, const NodeVector* foo, const unsigned int file_version){
            LOG(GlobalFlow::debug) << "Serializing vector";
            ar & foo->get()->size();
            for(auto it=foo->get()->begin(), end=foo->get()->end(); it!=end; ++it){
                ar & *it;
            }
        }

        template<class Archive>
        inline void load_construct_data(Archive & ar, NodeVector* foo, const unsigned int file_version){
            size_t size{0};
            ar & size;
            foo->get()->clear();
            foo->get()->reserve(size);
            for (int i = 0; i <size ; ++i) {
                std::unique_ptr<GlobalFlow::Model::NodeInterface> p;
                ar & p;
                foo->get()->push_back(std::move(p));
            }
        }
    }}


namespace GlobalFlow {
    namespace Simulation {
        using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<Model::NodeInterface>>>;
        using eq_ptr = std::unique_ptr<GlobalFlow::Solver::Equation>;
        using namespace boost::multiprecision;
        using namespace boost::units;
        namespace fs = boost::filesystem;

        /**
         * @class MassError
         * Simple container for the mass error calculations
         */
        class MassError {
        public:
            MassError(mpf_float_1000 OUT, mpf_float_1000 IN, mpf_float_1000 ERR) : OUT(OUT), IN(IN), ERR(ERR) {}

            mpf_float_1000 OUT;
            mpf_float_1000 IN;
            mpf_float_1000 ERR = 0;
        };

        /**
         * @class Simulation
         * The simulation class which holds the equation, options and data instance
         * Further contains methods for calculating the mass balance and sensitivity methods
         */
        class Simulation {
            eq_ptr eq;
            DataReader *reader;
            Options op;

        public:
            Simulation() {}
            Simulation(Options op, DataReader *reader); //impl in cpp - not pretty but works

            Solver::Equation *getEquation() { return eq.get(); };

            /**
             * Serialize current node state
             */
            void save() {
                if (serialize) {
                    LOG(stateinfo) << "Saving state for faster reboot..";
                    {
                        std::ofstream ofs(saveName,  std::ios::out |  std::ios::binary);
                        boost::archive::binary_oarchive outStream(ofs);
                        outStream << nodes;
                    }
                    LOG(stateinfo) << "Nodes saved";
                }
            };

            void restore(){
                if(loadNodes) {
                    LOG(stateinfo) << "Restoring state..";
                    {
                        std::ifstream in(saveName,  std::ios::in |  std::ios::binary);
                        boost::archive::binary_iarchive inStream(in);
                        inStream >> nodes;
                    }
                    LOG(stateinfo) << "Restored state successfully..";
                }
            }

            /**
             * Get basic node information by its id
             * @param nodeID
             * @return A string of information
             */
            std::string NodeInfosByID(unsigned long nodeID) {
                Model::NodeInterface *nodeInterface = nodes->at(nodeID).get();
                std::string out("\n");
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
                std::ostringstream output;
                output << to_string(in) << "\n";
                output << to_string(out);
                return output.str();
            }

            /**
             * Return all external flows separately
             */
            std::string NodeFlowsByID(unsigned long nodeID) {
                long id{0};
                for (int j = 0; j < nodes->size(); ++j) {
                    if (nodes->at(j)->getID() == nodeID) {
                        id = j;
                    }
                }
                Model::NodeInterface *nodeInterface = nodes->at(id).get();
                std::string out("");
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
             * Decide if its an In or Outflow
             * @param fun
             * @return
             */
            template<class Fun>
            MassError inline getError(Fun fun) {
                return getError([this, fun](int pos) {
                                    double flow = fun(pos);
                                    return flow < 0 ? flow : 0;
                                }, // fun1 for getError(fun1, fun2): if flow at pos is negative return it, else return 0
                                [this, fun](int pos) {
                                    double flow = fun(pos);
                                    return flow > 0 ? flow : 0;
                                }); // fun2 for getError(fun1, fun2): if flow at pos is positive return it, else return 0
            }

            /**
             * Calculate the mass error
             * @param fun1 Function to get OutFlow
             * @param fun2 Function to get InFlow
             * @return
             */
            template<class FunOut, class FunIn>
            MassError getError(FunOut fun1, FunIn fun2) {
                mpf_float_1000 out = 0;
                mpf_float_1000 in = 0;
                mpf_float_1000 error = 0;

                for (int j = 0; j < nodes->size(); ++j) {
                    out = out + fun1(j);
                    in = in + fun2(j);
                }
                if (abs(in - out) > 0.00001) {
                    error = ((100 * (in - abs(out))) / ((in + abs(out)) / 2));
                }
                MassError err(out, in, error);
                return err;
            }

            /**
             * Get the total mass balance
             * @return
             */
            MassError getMassError() {
                return getError([this](int pos) { return nodes->at(pos)->getOUT().value(); },
                                [this](int pos) { return nodes->at(pos)->getIN().value(); } );
            }

            /**
             * Get the mass balance for the current step
             * @return
             */
            MassError getCurrentMassError() {
                return getError([this](int pos) { return nodes->at(pos)->getCurrentOUT().value(); },
                                [this](int pos) { return nodes->at(pos)->getCurrentIN().value(); } );
            }

            /**
             * Get the mass balance for the current step
             * @return
             */
            MassError getVDFMassError() {
                return getError([this](int pos) { return nodes->at(pos)->getVDFOUT().value(); },
                                [this](int pos) { return nodes->at(pos)->getVDFIN().value(); } );
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
                GLOBAL_LAKES, // Question: GLOBAL_LAKES?
                WETLANDS,
                GLOBAL_WETLANDS,
                RECHARGE,
                FASTSURFACE,
                NAG,
                STORAGE,
                GENERAL_HEAD_BOUNDARY
            };

            /**
             * Helper function for printing the mass balance for each flow
             * @param flow
             * @return
             */
            std::string getMassErrorByFlowName(Flows flow) {
                std::ostringstream stream;
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
                    case GLOBAL_LAKES:
                        tmp = getError(
                                [this](int i) {
                                    return nodes->at(i)->getExternalFlowVolumeByName(Model::GLOBAL_LAKE).value();
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

                stream << "In: " << tmp.IN << "  Out: " << tmp.OUT;
                return stream.str();
            }

            /**
             * Prints all mass balances
             */
            void printMassBalances(custom_severity_level level) {
                LOG(level) << "All units in cubic meter per step size";
                MassError currentErr = getCurrentMassError();
                LOG(level) << "Step mass error: " << currentErr.ERR << "  In: " << currentErr.IN << "  Out: "
                           << currentErr.OUT;
                MassError totalErr = getMassError();
                LOG(level) << "Total mass error: " << totalErr.ERR << "  In: " << totalErr.IN << "  Out: "
                           << totalErr.OUT;
                MassError vdfErr = getVDFMassError();
                LOG(level) << "VDF mass error: " << vdfErr.ERR << "  In (sum over all zones): " << vdfErr.IN <<
                "  Out (sum over all zones): " << vdfErr.OUT;
                LOG(level) << "General Head Boundary: " << getMassErrorByFlowName(GENERAL_HEAD_BOUNDARY);
                LOG(level) << "Rivers: " << getMassErrorByFlowName(RIVERS);
                //LOG(stateinfo) << "Drains: " << getMassErrorByFlowName(DRAINS);
                LOG(level) << "Rivers MM: " << getMassErrorByFlowName(RIVER_MM);
                LOG(level) << "Lakes: " << getMassErrorByFlowName(LAKES);
                LOG(level) << "Global lakes: " << getMassErrorByFlowName(GLOBAL_LAKES);
                LOG(level) << "Wetlands: " << getMassErrorByFlowName(WETLANDS);
                LOG(level) << "Global wetlands: " << getMassErrorByFlowName(GLOBAL_WETLANDS);
                LOG(level) << "Recharge: " << getMassErrorByFlowName(RECHARGE);
                //LOG(userinfo) << "Fast Surface Runoff: " << getMassErrorByFlowName(FASTSURFACE) << "\n";
                LOG(level) << "Net abstraction from groundwater: " << getMassErrorByFlowName(NAG);
                LOG(level) << "Storage (only valid if transient run): " << getMassErrorByFlowName(STORAGE);
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
            void writeResiduals( std::string path ) {
                Eigen::Matrix<double, Eigen::Dynamic, 1> vector = eq->getResiduals();
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

            bool isRestored(){ return succefullyRestored; }

        private:
            NodeVector nodes;

            void initNodes() {
                LOG(userinfo) << FMAG(BOLD("Starting G³M"));
                nodes->reserve(op.getNumberOfNodesPerLayer() * op.getNumberOfLayers());
                reader->initNodes(nodes);
                reader->readData(op);
            };


            bool serialize{false};
            bool loadNodes{false};
            std::string saveName{".cached_simulation"};
            bool succefullyRestored{false};
        };

    }
}//ns
#endif