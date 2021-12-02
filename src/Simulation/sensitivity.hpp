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
#ifndef GLOBAL_FLOW_SENSITIVITY_HPP
#define GLOBAL_FLOW_SENSITIVITY_HPP

#include "../Misc/Helpers.hpp"

void scaleByIds(std::vector<int> ids, string field, double mult) {
    for (auto id : ids) {
        int i{-1};
        try {
            i = reader->getGlobIDMapping().at(id);
        } catch (const std::out_of_range &ex) {
            continue;
        }
        if (field == "RECHARGE") {
            try {
                Model::ExternalFlow flow = nodes->at(i)->getExternalFlowByName(Model::RECHARGE);
                double val = flow.getRecharge().value() * mult;
                nodes->at(i)->updateUniqueFlow(val, Model::RECHARGE);
            }
            catch (const std::out_of_range &ex) {
                continue;
            }
        }
    }
}

template<class Fun, class ChangeFunction>
void scaleByFunction(Fun fun, ChangeFunction apply) {
    for (int i = 0; i < nodes->size(); ++i) {
        if (fun(nodes->at(i))) {
            apply(nodes->at(i));
        }
    }
};

/**
 * Helper function for sensitivity
 * @param fun
 * @param field
 * @param mult
 */
template<class Fun>
void scaleByFunction(Fun fun, string field, double mult) {
    for (int i = 0; i < nodes->size(); ++i) {
        if (fun(nodes->at(i))) {
            if (field == "K") {
                quantity < Model::Velocity > k = nodes->at(
                        i)->getProperties().get < Model::quantity < Model::Velocity >, Model::K > ();
                nodes->at(i)->setK(k * mult);
            }
            if (field == "Anisotropy") {
                quantity < Model::Dimensionless > k =
                        nodes->at(i)->getProperties().get < quantity < Model::Dimensionless >,
                        Model::Anisotropy > ();
                nodes->at(i)->getProperties().set < quantity < Model::Dimensionless > ,
                        Model::Anisotropy > (k * mult);
            }
            if (field == "Depth_0") {
                if (nodes->at(i)->getProperties().get<int, Model::Layer>() == 0) {
                    quantity < Model::Meter > d = nodes->at(i)->getProperties().get < quantity < Model::Meter >,
                            Model::VerticalSize > ();
                    nodes->at(i)->getProperties().set < quantity < Model::Meter > , Model::VerticalSize >
                                                                                    (d * mult);
                    //TODO fix elevation for nodes below..
                }
            }
            if (field == "RECHARGE") {
                try {
                    Model::ExternalFlow flow = nodes->at(i)->getExternalFlowByName(Model::RECHARGE);
                    //LOG(debug) << "Multiplier: " << mult;
                    double val = flow.getRecharge().value() * mult;
                    nodes->at(i)->updateUniqueFlow(val, Model::RECHARGE);
                }
                catch (const std::out_of_range &ex) {
                    continue;
                }
            }
            if (field == "K_Wetlands") {
                try {
                    nodes->at(i)->updateExternalFlowConduct(mult, Model::WETLAND);
                }
                catch (const std::out_of_range &ex) {
                    continue;
                }
            }
            if (field == "K_Lakes") {
                try {
                    nodes->at(i)->updateExternalFlowConduct(mult, Model::LAKE);
                }
                catch (const std::out_of_range &ex) {
                    continue;
                }
            }
            if (field == "Lakes_bottom") {
                try {
                    nodes->at(i)->updateLakeBottoms(mult);
                }
                catch (const std::out_of_range &ex) {
                    continue;
                }
            }
            if (field == "River_Mult") {
                try {
                    nodes->at(i)->scaleDynamicRivers(mult);
                }
                catch (const std::out_of_range &ex) {
                    continue;
                }
            }
            if (field == "HEAD") {
                quantity < Model::Meter > m = nodes->at(
                        i)->getProperties().get < Model::quantity < Model::Meter >, Model::Head > ();
                nodes->at(i)->getProperties().set < Model::quantity < Model::Meter > , Model::Head > (m * mult);
            }

        }
    }

}


/**
 * Expects a file with global_ID, parameter, multiplier
 * Multiplier is log scaled
 * @param path to csv file
 * @returns a vector of node-ids (used to write out heads in correct order)
 */
std::vector<int> readMultipliersPerID(string path) {
    std::vector<int> ids;
    io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
    in.read_header(io::ignore_no_column, "global_ID", "mult", "parameter");
    int id;
    string parameter;
    string prev_parameter;
    bool n_param{false};
    double mult;
    while (in.read_row(id, mult, parameter)) {
        int i = reader->getGlobIDMapping().at(id);

        //Only save head ids once
        if (not n_param) {
            prev_parameter = parameter;
            n_param = true;
            ids.push_back(i);
        } else if (prev_parameter == parameter) {
            ids.push_back(i);
        }

        mult = pow(10, mult);
        NANChecker(mult, "Sensitivity multiplier");

        if (parameter == "K") {
            quantity < Model::Velocity > k = nodes->at(
                    i)->getProperties().get < Model::quantity < Model::Velocity >, Model::K > ();
            nodes->at(i)->setK(k * mult);
        }
        if (parameter == "K_Lakes") {
            try {
                nodes->at(i)->updateExternalFlowConduct(mult, Model::LAKE);
            }
            catch (const std::out_of_range &ex) {}
        }
        if (parameter == "K_Wetlands") {
            try {
                nodes->at(i)->updateExternalFlowConduct(mult, Model::WETLAND);
            }
            catch (const std::out_of_range &ex) {}
        }
        if (parameter == "K_Global_Wet") {
            try {
                nodes->at(i)->updateExternalFlowConduct(mult, Model::GLOBAL_WETLAND);
            }
            catch (const std::out_of_range &ex) {}
        }
        if (parameter == "RE") {
            try {
                Model::ExternalFlow flow = nodes->at(i)->getExternalFlowByName(Model::RECHARGE);
                double val = flow.getRecharge().value() * mult;
                nodes->at(i)->updateUniqueFlow(val, Model::RECHARGE);
            }
            catch (const std::out_of_range &ex) {}
        }
        if (parameter == "River_Mult") {
            try {
                nodes->at(i)->scaleDynamicRivers(mult);
            }
            catch (const std::out_of_range &ex) {}
        }
        if (parameter == "AqThickness") {
            quantity < Model::Meter > d = nodes->at(i)->getProperties().get < quantity < Model::Meter >,
                    Model::VerticalSize > ();
            nodes->at(i)->getProperties().set < quantity < Model::Meter > , Model::VerticalSize > (d * mult);
        }
        if (parameter == "SWB_Elevation") {
            try {
                nodes->at(i)->updateExternalFlowFlowHead(mult, Model::GLOBAL_WETLAND);
            } catch (const std::out_of_range &ex) {}
            try {
                nodes->at(i)->updateExternalFlowFlowHead(mult, Model::WETLAND);
            } catch (const std::out_of_range &ex) {}
            try {
                nodes->at(i)->updateExternalFlowFlowHead(mult, Model::LAKE);
            } catch (const std::out_of_range &ex) {}
            try {
                nodes->at(i)->updateExternalFlowFlowHead(mult, Model::RIVER_MM);
            } catch (const std::out_of_range &ex) {}
        }
    }
    return ids;
}

/**
 * Scale values for sensitivity analysis
 * @param path
 */
void readMultipliers(string path) {
    io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
    in.read_header(io::ignore_no_column, "name", "val");
    string name;
    double value;
    while (in.read_row(name, value)) {
        value = pow(10, value);
        NANChecker(value, "Sensitivity multiplier");
        if (name == "K") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "K", value);
        } else if (name == "RE") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "RECHARGE", value);
        } else if (name == "K_Wetlands") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "K_Wetlands", value);
        } else if (name == "K_Lakes") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "K_Lakes", value);
        } else if (name == "Lakes_bottom") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "Lakes_Bottom", value);
        } else if (name == "River_Mult") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "Dyn_River_Mult", value);
        } else if (name == "AqThickness") {
            scaleByFunction([](const std::unique_ptr<Model::NodeInterface> &node) {
                return true;
            }, "Depth0", value);
        }
    }
}


#endif
