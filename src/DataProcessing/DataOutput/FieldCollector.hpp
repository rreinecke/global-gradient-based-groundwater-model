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

#ifndef GLOBAL_FLOW_FIELDCOLLECTOR_HPP
#define GLOBAL_FLOW_FIELDCOLLECTOR_HPP

#include "Converter.hpp"

namespace GlobalFlow {
namespace DataProcessing {
namespace DataOutput {

using large_num = unsigned long int;

/**
 * Internal data fields that can be written out
 */
enum class FieldType {
        ID,
        ARCID,
        AREA,
        CONDUCT,
        ELEVATION,
        SLOPE,
        X, Y,
        HEAD,
        EQ_HEAD,
        IN, OUT,
        EQ_FLOW,
        LATERAL_FLOW,
        WETLANDS,
        LAKES,
        FLOW_HEAD,
        RECHARGE,
        DYN_RIVER,
        NODE_VELOCITY,
        RIVER_OUT,
        RIVER_IN,
        WTD,
        RIVER_CONDUCT,
        DRAIN_CONDUCT,
        WETLAND_CONDUCT,
        GL_WETLAND_CONDUCT,
        LAKE_CONDUCT,
        NON_VALID
};
const std::unordered_map<std::string, FieldType> fieldMapping{
        {"ID",                FieldType::ID},
        {"ArcID",             FieldType::ARCID},
        {"Area",              FieldType::AREA},
        {"Conductivity",      FieldType::CONDUCT},
        {"Elevation",         FieldType::ELEVATION},
        {"Slope",             FieldType::SLOPE},
        {"X",                 FieldType::X},
        {"Y",                 FieldType::Y},
        {"Head",              FieldType::HEAD},
        {"Eq_Head",           FieldType::EQ_HEAD},
        {"In",                FieldType::IN},
        {"Out",               FieldType::OUT},
        {"Eq_Flow",           FieldType::EQ_FLOW},
        {"Lateral_Flow",      FieldType::LATERAL_FLOW},
        {"Wetlands",          FieldType::WETLANDS},
        {"Lakes",             FieldType::LAKES},
        {"FlowHead",          FieldType::FLOW_HEAD},
        {"Recharge",          FieldType::RECHARGE},
        {"Dyn_River",         FieldType::DYN_RIVER},
        {"Velocity",          FieldType::NODE_VELOCITY},
        {"River_Out",         FieldType::RIVER_OUT},
        {"River_In",          FieldType::RIVER_IN},
        {"DepthToWaterTable", FieldType::WTD},
        {"River_Conduct",     FieldType::RIVER_CONDUCT},
        {"Drain_Conduct",     FieldType::DRAIN_CONDUCT},
        {"Wetland_Conduct",    FieldType::WETLAND_CONDUCT},
        {"Gl_Wetland_Conduct", FieldType::GL_WETLAND_CONDUCT},
        {"Lake_Conduct",       FieldType::LAKE_CONDUCT}
};


template<typename T> using data_vector = std::vector<T>;

/**
 * @class FieldCollector
 * Iterates over internal fields and searches for data to be written out
 * This is currently relatively inefficient
 */
class FieldCollector {
    private:
        const FieldType enumField;

        template<typename T, class Fun>
        data_vector<T> getData(Simulation::Simulation const &simulation, Fun fun) {
            data_vector<T> data;
            //TODO currently all layers
            for (int i = 0; i < simulation.nodes->size(); ++i) {
                data.push_back(static_cast<T>(fun(i)));
            }
            return data;
        }

    public:
        FieldCollector(FieldType enumField) : enumField(enumField) {}

        FieldCollector() : enumField(FieldType::NON_VALID) {}

        using pos_v = std::vector<std::pair<double, double>>;

        pos_v getPositions(Simulation::Simulation const &simulation) {
            return getData<std::pair<double, double>>(simulation, [simulation](int i) {
                double x{simulation.nodes->at(i)->getProperties().get<double, Model::Lat>()};
                double y{simulation.nodes->at(i)->getProperties().get<double, Model::Lon>()};
                return make_pair(x, y);
            });
        }

        std::vector<large_num> getIds(Simulation::Simulation const &simulation) {
            return getData<large_num>(simulation, [simulation](int i) {
                return simulation.nodes->at(i)->getID();
            });
        }

        /**
         * Collects the data from simulation nodes
         * @note Relatively inefficient currently
         * @param simulation The simulation
         * @return The collected data
         */
        template<typename T>
        data_vector<T> get(Simulation::Simulation const &simulation) {

            switch (enumField) {
                case FieldType::NON_VALID: {
                    throw std::bad_function_call();
                }
                case FieldType::ID : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties().get<large_num, Model::ID>());
                    });
                }
                case FieldType::ARCID : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties().get<large_num, Model::ArcID>());
                    });
                }
                case FieldType::AREA : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties()
                                                  .get < Model::quantity < Model::SquareMeter > ,
                                          Model::Area > ().value());
                    });
                }
                case FieldType::CONDUCT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getK().value());
                    });
                }
                case FieldType::ELEVATION : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties()
                                                  .get < Model::quantity < Model::Meter > ,
                                          Model::Elevation > ().value());
                    });
                }
                case FieldType::SLOPE : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties().get < Model::quantity <
                                          Model::Dimensionless > ,
                                          Model::Slope > ().value());
                    });
                }
                case FieldType::X : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties().get<double, Model::Lat>());
                    });
                }
                case FieldType::Y : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getProperties().get<double, Model::Lon>());
                    });
                }
                case FieldType::HEAD : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(
                                simulation.nodes->at(i)->getProperties().get < Model::quantity < Model::Meter > ,
                                Model::Head >
                                ().value());
                    });
                }
                case FieldType::EQ_HEAD : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(
                                simulation.nodes->at(i)->getProperties().get < Model::quantity < Model::Meter > ,
                                Model::EQHead >
                                ().value());
                    });
                }
                case FieldType::IN : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getIN().value());
                    });
                }
                case FieldType::OUT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getOUT().value());
                    });
                }
                case FieldType::EQ_FLOW : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getEqFlow().value());
                    });
                }
                case FieldType::LATERAL_FLOW : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>((simulation.nodes->at(i)->getLateralFlows().value() /
                                           simulation.nodes->at(
                                                   i)->getProperties().get < Model::quantity < Model::SquareMeter > ,
                                Model::Area >
                                ().value()) *
                                          1000);
                    });
                }
                case FieldType::WETLANDS : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->hasTypeOfExternalFlow(Model::WETLAND));
                    });
                }
                case FieldType::LAKES : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->hasTypeOfExternalFlow(Model::LAKE));
                    });
                }
                case FieldType::FLOW_HEAD: {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out = simulation.nodes->at(i)->getExternalFlowByName(Model::RIVER_MM).getFlowHead().value();
                        }
                        catch (exception &e) {
                        }
                        return convert<T>(out);
                    });
                }
                case FieldType::RECHARGE: {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out = simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RECHARGE).value();
                            out = (out / simulation.nodes->at(
                                    i)->getProperties().get < Model::quantity < Model::SquareMeter > , Model::Area >
                                                                                                       ().value()) *
                                  1000;
                        }
                        catch (exception &e) {
                        }
                        return convert<T>(out);
                    });
                }
                case FieldType::DYN_RIVER: {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out =
                                    simulation.nodes->at(i)->getExternalFlowByName(Model::RIVER_MM).getDyn(
                                            simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RECHARGE),
                                            simulation.nodes->at(i)->getProperties().get < Model::quantity <
                                            Model::Meter > ,
                                            Model::EQHead > (),
                                            simulation.nodes->at(
                                                    i)->getProperties().get < Model::quantity < Model::Meter > ,
                                            Model::Head > (),
                                            simulation.nodes->at(i)->getEqFlow()
                                    ).value();
                        }
                        catch (exception &e) {
                        }
                        return convert<T>(out);
                    });
                }
                case FieldType::NODE_VELOCITY: {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(simulation.nodes->at(i)->getVelocityVector());
                    });
                }
                case FieldType::RIVER_IN: {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{0};
                        try {
                            out += simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RIVER_MM).value();
                            out += simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RIVER).value();
                            out = (out / simulation.nodes->at(
                                    i)->getProperties().get < Model::quantity < Model::SquareMeter > , Model::Area >
                                                                                                       ().value()) *
                                  1000;
                        }
                        catch (exception &e) {
                        }
                        if (out < 0) {
                            out = 0;
                        }

                        return convert<T>(out);
                    });
                }
                case FieldType::RIVER_OUT: {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{0};
                        try {
                            out += simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RIVER_MM).value();
                            out += simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RIVER).value();
                            out = (out / simulation.nodes->at(
                                    i)->getProperties().get < Model::quantity < Model::SquareMeter > , Model::Area >
                                                                                                       ().value()) *
                                  1000;
                        }
                        catch (exception &e) {
                        }
                        if (out > 0) {
                            out = 0;
                        }

                        return convert<T>(out);
                    });
                }
                case FieldType::WTD : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        return convert<T>(
                                simulation.nodes->at(i)->getProperties().get < Model::quantity < Model::Meter > ,
                                Model::Elevation > ().value() -
                                                   simulation.nodes->at(i)->getProperties().get <
                                Model::quantity < Model::Meter > ,
                                Model::Head >
                                ().value());
                    });
                }
                case FieldType::RIVER_CONDUCT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out =
                                    simulation.nodes->at(i)->getExternalFlowByName(Model::RIVER_MM).getDyn(
                                            simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RECHARGE),
                                            simulation.nodes->at(i)->getProperties().get < Model::quantity <
                                            Model::Meter > , Model::EQHead > (),
                                            simulation.nodes->at(
                                                    i)->getProperties().get < Model::quantity < Model::Meter > ,
                                            Model::Head > (),
                                            simulation.nodes->at(i)->getEqFlow()
                                    ).value();
                        }
                        catch (exception &e) {}
                        return convert<T>(out);
                    });
                }
                case FieldType::DRAIN_CONDUCT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out =
                                    simulation.nodes->at(i)->getExternalFlowByName(Model::DRAIN).getDyn(
                                            simulation.nodes->at(i)->getExternalFlowVolumeByName(Model::RECHARGE),
                                            simulation.nodes->at(i)->getProperties().get < Model::quantity <
                                            Model::Meter > , Model::EQHead > (),
                                            simulation.nodes->at(
                                                    i)->getProperties().get < Model::quantity < Model::Meter > ,
                                            Model::Head > (),
                                            simulation.nodes->at(i)->getEqFlow()
                                    ).value();
                        }
                        catch (exception &e) {}
                        return convert<T>(out);
                    });
                }
                case FieldType::WETLAND_CONDUCT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out = simulation.nodes->at(i)->getExternalFlowByName(
                                    Model::WETLAND).getConductance().value();
                        }
                        catch (exception &e) {}
                        return convert<T>(out);
                    });
                }
                case FieldType::GL_WETLAND_CONDUCT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out = simulation.nodes->at(i)->getExternalFlowByName(
                                    Model::GLOBAL_WETLAND).getConductance().value();
                        }
                        catch (exception &e) {}
                        return convert<T>(out);
                    });
                }
                case FieldType::LAKE_CONDUCT : {
                    return getData<T>(simulation, [simulation, this](int i) {
                        double out{NAN};
                        try {
                            out = simulation.nodes->at(i)->getExternalFlowByName(Model::LAKE).getConductance().value();
                        }
                        catch (exception &e) {}
                        return convert<T>(out);
                    });
                }
                default:
                    throw std::bad_function_call();
            }
        }
};
}
}
}
#endif