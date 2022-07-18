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

#ifndef GLOBAL_FLOW_NODE_HPP
#define GLOBAL_FLOW_NODE_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <numeric>
#include <future>
#include <forward_list>
#include "ExternalFlows.hpp"
#include "FluidMechanics.hpp"
#include "PhysicalProperties.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {
/**
 * Neighbouring positions for cells
 *     TOP
 * LEFT * RIGHT
 *     DOWN
 *
 * In Z (Top view):
 * FRONT (larger ID) * BACK (smaler ID)
 */
        enum NeighbourPosition {
            TOP = 1,
            DOWN,
            RIGHT,
            LEFT,
            FRONT,
            BACK
        };
    }
}

namespace std {
    template<>
    struct hash<GlobalFlow::Model::NeighbourPosition> {
        typedef GlobalFlow::Model::NeighbourPosition argument_type;
        typedef std::underlying_type<argument_type>::type underlying_type;
        typedef std::hash<underlying_type>::result_type result_type;

        result_type
        operator()(const argument_type &arg) const {
            std::hash<underlying_type> hasher;
            return hasher(static_cast< underlying_type >( arg ));
        }
    };
}

namespace GlobalFlow {
    namespace Model {

        using namespace std;
        using namespace boost::units;

        class NodeInterface;

        using p_node = unique_ptr<GlobalFlow::Model::NodeInterface>;
        using NodeVector = std::shared_ptr<vector<unique_ptr<GlobalFlow::Model::NodeInterface>>>;

/**
 * Interface defining required fields for a node.
 * A node is the central computational and spatial unit.
 * A simulated area is seperated into a discrete raster of cells or nodes
 * (separate computational units which stay in contact to ech other).
 * Is equal to 'cell'.
 *
 * Nodes can be of different physical property e.g. different size.
 */
        class NodeInterface {
        private:
            virtual void
            __setHead(t_meter head) = 0;

            virtual t_meter
            __calcInitialHead(t_meter initialParam) = 0;

            virtual bool
            __isStaticNode() = 0;

            /**
             * @brief Calculates e-folding depth from input data.
             * @param z The vertical size of the computing cell.
             * @return A dimensionless factor that can be used to modify
             * hydraulic conductance depending on depth.
             *
             * E-Folding function as defined by Ying Fan et al. e^(-(Depth - Factor)).
             */
            t_dim efoldingFromData(t_meter z) {
                t_meter folding = get<t_meter, EFolding>();
                if (folding == 0.0 * si::meter)
                    return 1 * si::si_dimensionless;
                //Alter if a different size should be used and not full vertical size
                //z = (z / (2 * si::si_dimensionless));
                t_dim out = exp(-z / folding);
                if (out == 0 * si::si_dimensionless)
                    return 1e-7 * si::si_dimensionless;
                return out;
            }

            /**
             * @brief Apply function to all layers of the model.
             * @param &&...p A list of functions forwarded to setAttribute.
             *
             * Apply an attribute to all layers at same position as the node.
             * A function passed in addition modifies the attribute depending on the layer.
             * Stops on last layer.
             */
            template<class... Params>
            void
            applyToAllLayers(Params &&...p) {
                try {
                    this->getNeighbour(DOWN)->setAttribute(std::forward<Params>(p)...);
                } catch (...) {} //We are on last layer nothing to do
            };

            template<class ApplyFunction>
            void
            setAttribute(ApplyFunction fun) {
                fun(this);
                applyToAllLayers(std::forward<ApplyFunction>(fun));
            };

        protected:
            const std::shared_ptr<std::vector<std::unique_ptr<NodeInterface>>> nodes;
            unordered_map<NeighbourPosition, large_num> neighbours;
            unordered_map<FlowType, ExternalFlow, FlowTypeHash> externalFlows;
            vector<t_meter> Zetas; // zeta surfaces (dimensionless density and elevation) of current node at t
            vector<t_meter> ZetasZero; // zeta surfaces (dimensionless density and elevation) of current node at t-1
            vector<t_dim> zones; // dimensionless density zones in current node | todo: do we need this or is it enough to access generally with densityProps.nusZones?
            int numOfExternalFlows{0};
            int numOfZetas{0};
            int numOfZones{0};
            bool nwt{false};
            bool initial_head{true};
            bool simpleDistance{false};
            bool simpleK{false};
            bool steadyState{false};

            friend class FluidMechanics;
            friend class VariableDensityFlow;

            FluidMechanics mechanics;
            VariableDensityFlow vdf;
            PhysicalProperties fields;

            /*
            * Short helper functions to make the code more readable
            */
            template<typename T, typename F>
            T get() { return fields.get<T, F>(); }

            template<typename T, typename F>
            void set(const T &value) { return fields.set<T, F>(value); }

            template<typename T, typename F>
            T getFrom(NodeInterface *nodeInterface, NeighbourPosition pos) {
                return nodeInterface->getNeighbour(pos)->getProperties().get<T, F>();
            }

            using map_itter = std::unordered_map<NeighbourPosition, large_num>::const_iterator;

            p_node &at(map_itter pos) { return nodes->at(pos->second); }

            template<typename T, typename F>
            T getAt(map_itter pos) {
                return at(pos)->get<T, F>();
            }

            /**
             * @brief Uses specific storage to calculate storativity
             * @return Flow budget for cell depending on head change
             */
            t_s_meter getStorageCapacity__Primary() noexcept {
                t_s_meter out = get<quantity<perUnit>, SpecificStorage>().value()
                                                                       * get<t_c_meter, VolumeOfCell>().value() * si::square_meter;
                return out;
            }

            /**
             * @brief Uses specific yield to calculate storativity
             * @return Flow budget for cell depending on head change
             */
            t_s_meter getStorageCapacity__Secondary() noexcept {
                return get<t_dim, SpecificYield>() * get<t_s_meter, Area>();
            }

            /**
             * @brief Uses specific yield with NWT smoother
             * @return Flow budget for cell depending on head change
             * Assumes that smooth function is linear
             */
            t_s_meter getStorageCapacity__SecondaryNWT() noexcept {
                return getStorageCapacity__Secondary()
                       * mechanics.smoothFunction__NWT(get<t_meter, Elevation>(), get<t_meter, VerticalSize>(),
                                                       get<t_meter, Head>())
                       + getStorageCapacity__Secondary();
            }

            /**
             * @brief Flow volume of cell
             * @return Flow volume
             */
            template<typename CompareFunction>
            t_vol_t getFlow(CompareFunction compare) noexcept {
                t_vol_t flow_v = 0.0 * si::cubic_meter / day;
                t_vol_t stFlow = 0.0 * si::cubic_meter / day;
                if (not steadyState) {
                    stFlow = getTotalStorageFlow();
                    if (compare(stFlow.value()))
                        flow_v += boost::units::abs(stFlow);
                }
                flow_v += getNonStorageFlow(std::forward<CompareFunction>(compare));
                return flow_v.value() * si::cubic_meter / day;
            }

        public:

            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version) {
                LOG(debug) << "Serializing abstract node";
                ar & nodes;
                ar & neighbours;
                ar & externalFlows;
                ar & numOfExternalFlows;
                ar & Zetas;
                ar & numOfZetas;
                ar & zones;
                ar & numOfZones;
                ar & nwt;
                ar & initial_head;
                ar & simpleDistance;
                ar & simpleK;
                ar & steadyState;
                ar & fields;
                ar & cached;
                ar & eq_flow;
            }


            /**
             * @brief Constructor of abstract class NodeInterface
             * @param nodes Vector of all other existing nodes
             * @param lat The latitude
             * @param lon The Longitude
             * @param area Area in m²
             * @param ArcID Unique ARC-ID specified by Kassel
             * @param ID Internal ID = Position in vector
             * @param K Hydraulic conductivity in meter/day (default)
             * @param stepModifier Modifies default step size of day (default=1)
             * @param aquiferDepth Vertical size of the cell
             * @param anisotropy Modifier for vertical conductivity based on horizontal
             * @param specificYield Yield of storage for dewatered conditions
             * @param specificStorage Specific storage - currently for confined and unconfined
             * @param confined Is node in a confined layer
             *
             */
            NodeInterface(NodeVector nodes,
                          double lat,
                          double lon,
                          t_s_meter area,
                          t_meter edgeLengthLeftRight,
                          t_meter edgeLengthFrontBack,
                          large_num ArcID,
                          large_num ID,
                          t_vel K,
                          int stepModifier,
                          double aquiferDepth,
                          double anisotropy,
                          double specificYield,
                          double specificStorage,
                          bool confined,
                          DensityProperties densityProperties);

            virtual ~NodeInterface() = default;

            large_num getID() { return get<large_num, ArcID>(); }

/*****************************************************************
Modify Properties
******************************************************************/

            /**
             * @brief Set elevation on top layer and propagate to lower layers
             * @param elevation The top elevation (e.g. from DEM)
             */
            void setElevation(t_meter elevation) {
                set < t_meter, Elevation > (elevation);
                set < t_meter, TopElevation > (elevation);
                applyToAllLayers([this](NodeInterface *node) {
                    try {
                        node->getProperties().set<t_meter, TopElevation>(getFrom<t_meter, TopElevation>(node, TOP));
                        node->getProperties().set<t_meter, Elevation>(
                                getFrom<t_meter, Elevation>(node, TOP) - getFrom<t_meter, VerticalSize>(node, TOP));
                    }
                    catch (...) {}
                });
            };

            /**
             * @brief Set slope from data on all layers
             * Slope input is in % but is required as absolut
             * thus: slope = slope_percent / 100
             * @param slope
             */
            void
            setSlope(double slope_percent) {
                set < t_dim, Slope > ((slope_percent / 100) * si::si_dimensionless);
                applyToAllLayers([slope_percent](NodeInterface *nodeInterface) {
                    try { // todo: should slope be added to the nodeInterace? currently only in PhyscalProperties
                        //nodeInterface->
                        //Slope(slope_percent);
                    }
                    catch (...) {}
                });
            };

            /**
             * @brief Set e-folding factor from data on all layers
             * @param e-fold
             */
            void setEfold(double efold) {
                set < t_meter, EFolding > (efold * si::meter);
                applyToAllLayers([efold](NodeInterface *nodeInterface) {
                    try { nodeInterface->setEfold(efold); }
                    catch (...) {}
                });
            };

            /**
             * @brief Calculated equilibrium groundwater-head from eq_wtd
             * Assumes that if initialhead = false that the eq_head is also used as initial head
             * @param head
             */
            void setEqHead(t_meter wtd) {
                t_meter eqhead = get<t_meter, Elevation>() - wtd;
                //if (eqhead.value() > 500) {
                //	eqhead = 500 * si::meter;
                //}
                set < t_meter, EQHead > (eqhead);
                setHead_direct(eqhead.value());
                //setHead_direct(get<t_meter, Elevation>().value());
                applyToAllLayers([eqhead](NodeInterface *nodeInterface) {
                    try {
                        nodeInterface->setHead_direct(eqhead.value());
                        nodeInterface->set<t_meter, EQHead>(eqhead);
                    }
                    catch (...) {}
                });
            }

            /**
             * Calculated equilibrium flow to neighbouring cells
             * Static thus calculated only once.
             *
             * Depends on: K in cell and eq_head in all 6 neighbours
             */
            bool cached{false};
            t_vol_t eq_flow{0 * si::cubic_meter / day};

            template<class HeadType>
            FlowInputHor createDataTuple(map_itter got) {
                t_meter edgeLength_self = 0 * si::meter;
                t_meter edgeLength_neig = 0 * si::meter;
                t_meter edgeWidth_self = 0 * si::meter;
                if (got->first == LEFT or got->first == RIGHT){
                    edgeLength_self = get<t_meter, EdgeLengthFrontBack>();
                    edgeLength_neig = getAt<t_meter, EdgeLengthFrontBack>(got);
                    edgeWidth_self = get<t_meter, EdgeLengthLeftRight>();
                } else { // if (got->first == FRONT or got->first == BACK)
                    edgeLength_self = get<t_meter, EdgeLengthLeftRight>();
                    edgeLength_neig = getAt<t_meter, EdgeLengthLeftRight>(got);
                    edgeWidth_self = get<t_meter, EdgeLengthFrontBack>();
                }
                return std::make_tuple(at(got)->getK(),
                                       getK(),
                                       edgeLength_neig, // length in front/back direction of neighbour node
                                       edgeLength_self, // length in front/back direction of this node
                        //getAt<t_meter, EdgeLengthLeftRight>(got), // width in left/right direction of neighbour node
                                       edgeWidth_self, // width in left/right direction of this node
                                       getAt<t_meter, HeadType>(got),
                                       get<t_meter, HeadType>(),
                                       getAt<t_meter, Elevation>(got),
                                       get<t_meter, Elevation>(),
                                       getAt<t_meter, VerticalSize>(got),
                                       get<t_meter, VerticalSize>(),
                                       get<bool, Confinement>());
            }

            FlowInputVert createDataTuple(map_itter got) {
                return std::make_tuple(at(got)->getK_vertical(),
                                       getK_vertical(),
                                       get<t_meter, VerticalSize>(),
                                       get<t_meter, Head>(),
                                       get<t_meter, Elevation>(),
                                       get<t_s_meter, Area>(),
                                       getAt<t_meter, Elevation>(got),
                                       getAt<t_meter, VerticalSize>(got),
                                       getAt<t_meter, Head>(got),
                                       get<bool, Confinement>());
            }

            /**
             * Calculate the lateral groundwater flow to the neighbouring nodes
             * Generic function used for calculating equilibrium and current step flow
             * @return
             */
            template<class HeadType>
            t_vol_t calcLateralFlows(bool onlyOut) {
                t_vol_t lateral_flow{0 * si::cubic_meter / day};
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT, NeighbourPosition::LEFT,
                         NeighbourPosition::RIGHT};

                for (const auto &position: possible_neighbours) {
                    std::unordered_map<NeighbourPosition, large_num>::const_iterator got = neighbours.find(position);
                    if (got == neighbours.end()) {//No neighbouring node
                    } else {
                        //There is a neighbour node
                        t_s_meter_t conductance;
                        //TODO check for option!
                        //if (get<int, Layer>() > 0) {
                        //    conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(got), get<t_meter, EFolding>(), getAt<t_meter, EFolding>(got));
                        //}else{
                        conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got));
                        //}

                        /*if (Zetas.size() > 0) { // Question: add fluxCorrection here?
                            t_meter fluxCorrection;
                            fluxCorrection = verticalFluxCorrectionTerm();

                            t_vol_t flow = conductance * (get<t_meter, HeadType>() - getAt<t_meter, HeadType>(got) -
                                                          fluxCorrection);
                        } else {*/
                        t_vol_t flow = conductance * (get<t_meter, HeadType>() - getAt<t_meter, HeadType>(got));
                        //}
                        if (onlyOut) {
                            if (flow.value() > 0) { lateral_flow = lateral_flow - flow; }
                        } else { lateral_flow = lateral_flow - flow; }
                    }
                }
                return lateral_flow;
            }

            /**
             * Calculate the equilibrium lateral flows
             * @return eq lateral flow
             */
            t_vol_t getEqFlow() noexcept {
                if (not cached) {
                    t_vol_t lateral_flow = calcLateralFlows<EQHead>(false);
                    NANChecker(lateral_flow.value(), "Eq Flow");
                    eq_flow = lateral_flow;
                    cached = true;
                }
                return eq_flow;
            }

            /**
             * Get the current lateral flow
             * @return
             */
            t_vol_t getLateralFlows() {
                return calcLateralFlows<Head>(false) * get<t_dim, StepModifier>();
            }

            /**
             * Get the current lateral out flows
             * @return
             */
            t_vol_t getLateralOutFlows() {
                return calcLateralFlows<Head>(true) * get<t_dim, StepModifier>();
            }


            /**
             * @brief Cuts off all heads above surface elevation
             * @warning Should only be used in spin up phase!
             * @return Bool if node was reset
             */
            bool resetFloodingHead() noexcept {
                auto elevation = get<t_meter, Elevation>();
                if (get<t_meter, Head>() > elevation) {
                    set < t_meter, Head > (elevation);
                    return true;
                }
                return false;
            }

            /**
             * @brief Scales river conduct by 50%
             * @warning Should only be used in spin up phase
             */
            void scaleRiverConduct() {
                eq_flow = eq_flow * 1.5;
            }

            /**
             * @brief Update the current head change (in comparison to last time step)
             * @note Should only be called at end of time step
             */
            void updateHeadChange() noexcept {
                set < t_meter, HeadChange_TZero > (
                        get<t_meter, Head>() - get<t_meter, Head_TZero>());
                set < t_meter, Head_TZero > (get<t_meter, Head>());
            }

            void initHead_t0() noexcept { set < t_meter, Head_TZero > (get<t_meter, Head>()); }

            void setHead_direct(double head) noexcept { set < t_meter, Head > (head * si::meter); }

            t_vel getK__pure() noexcept { return get<t_vel, K>(); }

            /**
             * @brief Get hydraulic conductivity
             * @return hydraulic conductivity (scaled by e-folding)
             */
            t_vel getK() noexcept {
                if (simpleK) { return get<t_vel, K>() * get<t_dim, StepModifier>(); }
                t_dim e_fold = 1 * si::si_dimensionless;
                if (get<int, Layer>() > 0) {
                    e_fold = efoldingFromData(get<t_meter, VerticalSize>());
                }
                t_vel out = get<t_vel, K>() * e_fold * get<t_dim, StepModifier>();
                if (out < 1e-20 * si::meter / day) {
                    out = 1e-20 * si::meter / day;
                }
                return out;
            }

            /**
             * @brief Get hydraulic vertical conductivity
             * @return hydraulic conductivity scaled by anisotropy (scaled by e-folding)
             */
            t_vel getK_vertical() noexcept { return (getK() / get<t_dim, Anisotropy>()) * get<t_dim, StepModifier>(); }

            /**
             * @brief Modify hydraulic conductivity (applied to all layers below)
             * @param New conductivity (if e-folding enabled scaled on layers)
             */
            void setK(t_vel conduct) {
                setK_direct(conduct);
                applyToAllLayers([&conduct](NodeInterface *nodeInterface) {
                    nodeInterface->getProperties().set<t_vel, K>(conduct);
                });
            }

            /**
             * @brief Modify hydraulic conductivity (no e-folding, no layers)
             * @param New conductivity
             */
            void setK_direct(t_vel conduct) { set < t_vel, K > (conduct); }

            void setSimpleK(){simpleK = true;}

            /**
             * @brief Get all outflow since simulation start
             */
            t_c_meter getOUT() noexcept { return get<t_c_meter, OUT>(); }

            /**
             * @brief Get all inflow since simulation start
             */
            t_c_meter getIN() noexcept { return get<t_c_meter, IN>(); }

            /**
             * @brief Toogle steady state simulation
             * @param onOFF true=on
             * Turns all storage equations to zero with no time steps
             */
            void toggleSteadyState(bool onOFF) { this->steadyState = onOFF; }

            void updateStepSize(double mod) { set < t_dim, StepModifier > (mod * si::si_dimensionless); }

            /**
             * @brief Storage capacity based on yield or specific storage
             * @return Potential flow budget when multiplied by head change
             * Uses an 0.001m epsilon to determine if a water-table condition is present.
             * If the layer is confined or not in water-table condition returns primary capacity.
             */
            t_s_meter getStorageCapacity() noexcept {
                t_meter epsilon = 0.001 * si::meter;
                //TODO make this an option!
                if (get<int, Layer>() == 0) {
                    //If this is the first layer always use equation for unconfined storage
                    //if (get<bool, Confinement>()) { return getStorageCapacity__Primary(); }
                    //else { return getStorageCapacity__Secondary();}
                    return getStorageCapacity__Secondary();
                } else {
                    //when we are in the second layer check if we are defined as confined
                    if (get<bool, Confinement>()) { return getStorageCapacity__Primary(); }
                }
                //we are not in the first layer and we are in unconfined conditions as specified by the user

                if (get<t_meter, Head>() + epsilon < get<t_meter, Elevation>()) {
                    //water-table condition
                    if (nwt) { return getStorageCapacity__SecondaryNWT(); }
                    return getStorageCapacity__Secondary();
                } else {
                    return getStorageCapacity__Primary();
                }
            }

            /**
             * @brief Get an external flow by its FlowType
             * @param type The flow type
             * @return Ref to external flow
             * @throw OutOfRangeException
             */
            ExternalFlow &getExternalFlowByName(FlowType type) {
                if (externalFlows.find(type) == externalFlows.end())
                    throw out_of_range("No such flow");
                return externalFlows.at(type);
            }

            /**
             * @brief Get an external flow volume by its FlowType
             * @param type The flow type
             * @return Flow volume
             */
            t_vol_t getExternalFlowVolumeByName(FlowType type) {
                for (const auto &flow : externalFlows) {
                    if (type == flow.second.getType()) {
                        t_vol_t out = calculateExternalFlowVolume(flow.second);
                        return out;
                    }
                }
                return 0 * si::cubic_meter / day;
            }

            /**
             * @brief Get flow budget based on head change
             * @return Flow volume
             * Note: Water entering storage is treated as an outflow (-), that is a loss of water from the flow system
             * while water released from storage is treated as inflow (+), that is a source of water to the flow system
             */
            t_vol_t getTotalStorageFlow() noexcept {
                return -getStorageCapacity() * get<t_meter, HeadChange_TZero>() / (day * get<t_dim, StepModifier>());
            }

            /**
             * @brief Get flow budget of a specific external flows
             * @param &flow A external flow
             * @return Flow volume
             * Note: Water entering storage is treated as an outflow (-), that is a loss of water from the flow system
             * while water released from storage is treated as inflow (+), that is a source of water to the flow system
             */
            t_vol_t calculateExternalFlowVolume(const ExternalFlow &flow) {
                if (is(flow.getType()).in(RECHARGE, NET_ABSTRACTION)) {
                    return flow.getRecharge() * get<t_dim, StepModifier>();
                }
                t_vol_t ex;
                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {
                    //ignore me there is no special_flow in this cell
                }
                t_dim slope = get<t_dim, Slope>();
                t_vol_t eqFlow = getEqFlow(); // get the equilibrium lateral flows
                if (is(flow.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, WETLAND, GLOBAL_WETLAND)) {
                    if (flow.flowIsHeadDependant(head)) {
                        ex = (flow.getP(eq_head, head, recharge, slope, eqFlow) * head +
                              flow.getQ(eq_head, head, recharge, slope, eqFlow)) * get<t_dim, StepModifier>();
                    } else { // flow is not head dependent when the head is below the bottom of the simulated cell
                        ex = (flow.getP(eq_head, head, recharge, slope, eqFlow) * flow.getBottom() +
                              flow.getQ(eq_head, head, recharge, slope, eqFlow)) * get<t_dim, StepModifier>();
                    }
                } else {
                    ex = (flow.getP(eq_head, head, recharge, slope, eqFlow) * head +
                          flow.getQ(eq_head, head, recharge, slope, eqFlow)) * get<t_dim, StepModifier>();
                }
                return ex;
            }

            /**
             * @brief Calculate dewatered flow
             * @return Flow volume per time
             * If a cell is dewatered but below a saturated or partly saturated cell:
             * this calculates the needed additional exchange volume
             */
            t_vol_t calculateDewateredFlow() noexcept {
                map_itter hasDown = neighbours.find(DOWN);
                map_itter hasUp = neighbours.find(TOP);
                t_vol_t out = 0 * si::cubic_meter / day;

                if (hasDown != neighbours.end()) {
                    t_meter elev = getAt<t_meter, Elevation>(hasDown);
                    t_meter head_n = getAt<t_meter, Head>(hasDown);
                    //Check if a dewatered condition is present
                    if (head_n < elev and get<t_meter, Head>() > elev) {
                        t_s_meter_t conductance_below = mechanics.calculateVerticalConductance(createDataTuple(hasDown));
                        out += conductance_below * (head_n - elev) * get<t_dim, StepModifier>();
                    }
                }

                if (hasUp != neighbours.end()) {
                    t_meter elev = getAt<t_meter, Elevation>(hasUp);
                    t_meter head_n = getAt<t_meter, Head>(hasUp);
                    //Check if a dewatered condition is present
                    if (get<t_meter, Head>() < get<t_meter, Elevation>() and head_n > elev) {
                        t_s_meter_t conductance_above =
                                mechanics.calculateVerticalConductance(createDataTuple(hasUp));
                        out += conductance_above * (get<t_meter, Elevation>() - get<t_meter, Head>()) *
                               get<t_dim, StepModifier>();
                    }
                }
                NANChecker(out.value(), "Dewatered flow");
                return out;
            }

            /**
             * @brief Get all current IN flow
             * @return Flow volume
             */
            t_vol_t getCurrentIN() noexcept { return getFlow([](double a) -> bool { return a > 0; }); }

            /**
             * @brief Get all current OUT flow
             * @return Flow volume
             */
            t_vol_t getCurrentOUT() noexcept { return getFlow([](double a) -> bool { return a < 0; }); }

            /**
             * @brief Tell cell to save its flow budget
             */
            void saveMassBalance() noexcept {
                fields.addTo<t_c_meter, OUT>(getCurrentOUT().value() * si::cubic_meter);
                fields.addTo<t_c_meter, IN>(getCurrentIN().value() * si::cubic_meter);
            }

            /**
             * @brief Add a neighbour
             * @param ID The internal ID and position in vector
             * @param neighbour The position relative to the cell
             */
            void setNeighbour(large_num ID, NeighbourPosition neighbour) { neighbours[neighbour] = ID; }

            int getNumofNeighbours() { return (int) neighbours.size(); }

            class NodeNotFoundException : public std::exception {
                virtual const char *what() const throw() { return "Node does not exist"; }
            };

            unordered_map<NeighbourPosition, large_num> getListOfNeighbours(){
                return neighbours;
            }
            /**
             * @brief Get a neighbour by position
             * @param neighbour The position relative to the cell
             * @return Pointer to cell object
             */
            NodeInterface *getNeighbour(NeighbourPosition neighbour) noexcept(false) {
                try {
                    large_num pos = neighbours.at(neighbour);
                    if (nodes->at(pos)->get<large_num, ID>() != pos) { throw NodeNotFoundException(); }
                    return nodes->at(pos).get();
                } catch (...) {
                    throw NodeNotFoundException();
                }
            }

            /**
             * @brief Add an external flow to the cell
             * @param type The flow type
             * @param flowHead The flow head
             * @param cond The conductance
             * @param bottom The bottom of the flow (e.g river bottom)
             * @return Number assigned by cell to flow
             */
            int addExternalFlow(FlowType type, t_meter flowHead, double cond, t_meter bottom) {
                if (hasTypeOfExternalFlow(type)) {
                    //LOG(debug) << "! adding a flow that is already existing for cell"
                    //currently it is assumed that only one external flow of one type is what we want
                    // FIXME if not we have to replace the enum with something different
                    removeExternalFlow(type);
                }

                if (type == RECHARGE or type == FAST_SURFACE_RUNOFF or type == NET_ABSTRACTION) {
                    externalFlows.insert(std::make_pair(type,
                                                        ExternalFlow(numOfExternalFlows, cond * (si::cubic_meter / day),
                                                                     type)));
                } else if (type == EVAPOTRANSPIRATION) {
                    externalFlows.insert(std::make_pair(type,
                                                        ExternalFlow(numOfExternalFlows, flowHead, bottom,
                                                                     cond * (si::cubic_meter / day))));
                }
                    // TODO Implementation of FLOODPLAIN_DRAIN
                    /* else if (type == FLOODPLAIN_DRAIN) {
                        externalFlows.insert(std::make_pair(type,
                                                            ExternalFlow(numOfExternalFlows, type,
                                                                            get<t_meter, Elevation>(),
                                                                            get<t_vel, K>() * get<t_meter,
                                                                            VerticalSize>(),
                                                                            bottom));

                    } */ else { // RIVER, RIVER_MM, DRAIN, WETLAND, GLOBAL_WETLAND, LAKE, GENERAL_HEAD_BOUNDARY
                    externalFlows.insert(std::make_pair(type,
                                                        ExternalFlow(numOfExternalFlows,
                                                                     type,
                                                                     flowHead,
                                                                     cond * (si::square_meter / day),
                                                                     bottom)));
                }
                numOfExternalFlows++;
                if(numOfExternalFlows != externalFlows.size()){
                    LOG(debug) << "Printing flows ";
                    for(auto const& imap: externalFlows)
                        LOG(debug) << " " << imap.first;
                    throw "Number of external flows don't match";
                }
                return numOfExternalFlows;
            }

            /**
             * @brief Remove an external flow to the cell by id
             * @param ID The flow id
             */
            void removeExternalFlow(FlowType type) {
                if (externalFlows.erase(type)) {
                    numOfExternalFlows = numOfExternalFlows - 1;
                }
                if(numOfExternalFlows != externalFlows.size()){
                    LOG(debug) << "Printing flows ";
                    for(auto const& imap: externalFlows)
                        LOG(debug) << " " << imap.first;
                    throw "Number of external flows don't match";
                }
            }

            /**
             * @brief The number of external flows
             */
            int getNumOfExternalFlows(){return numOfExternalFlows;}


            /**
             * @brief Check for an external flow by type
             * @param type The flow type
             * @return bool
             */
            bool hasTypeOfExternalFlow(FlowType type) {
                if(externalFlows.find(type) == externalFlows.end()){
                    return false;
                }
                return true;
            }

            /**
            * @brief Add a density zone to the cell
            * @param nus dimensionless density of the zone
            * @return number of density zones in the cell
            */
            int addZone(t_dim nus){
                zones.push_back(nus);
                sort(zones.begin(), zones.end());
                numOfZones++;
                return numOfZones;
            }

            // todo: removeZone like removeZeta below

            /**
            * @brief Add a zeta surface to the cell
            * @param height the zeta surface height in meters
            * @return number of zeta surfaces in the cell
            */
            int addZetaSurface(t_meter height){
                Zetas.push_back(height);
                sort(Zetas.begin(), Zetas.end(), greater<t_meter>());
                numOfZetas++;
                return numOfZetas;
            }

            /**
            * @brief Set the height of zeta surface n
            * @param n the zeta ID
            * @param height the zeta surface height in meters
            */
            void setZeta(int n, t_meter height){
                if (n > Zetas.size() - 1){
                    Zetas.push_back(height);
                } else {
                    Zetas[n] = height;
                }
            }

            /**
            * @brief Remove a zeta surface of the cell by zeta id
            * @param zeta_ID The zeta id
            */
            void removeZeta(int id) {
                // Todo: caution! this may delete a surface between two surface | add exception
                Zetas.erase(Zetas.begin() + id);
                numOfZetas = numOfZetas - 1;
            }

            /**
            * @brief The number of density surfaces
            */
            int getNumOfZetas() { return (int) numOfZetas;}

            /**
            * @brief get zeta surfaces
            */
            std::vector<t_meter> getZetas() noexcept { return Zetas;}

            /**
             * @brief set zeta surfaces at t-1
             */
            void setZetasZero(std::vector<t_meter> zetas){
                ZetasZero.clear();
                for (int n = 0; n < zetas.size() - 1; n++){
                    ZetasZero.push_back(zetas[n]);
                }
            }

            /**
            * @brief get density zones
            */
            std::vector<t_dim> getZones() noexcept { return zones;}

            /**
             * @brief Modify effective porosity (applied to all layers below)
             * @param New effective porosity
             */
            void setEffectivePorosity(t_dim effectivePorosity) {
                setEffectivePorosity_direct(effectivePorosity);
                applyToAllLayers([&effectivePorosity](NodeInterface *nodeInterface) {
                    nodeInterface->getProperties().set<t_dim, EffectivePorosity>(effectivePorosity);
                });
            }

            /**
             * @brief Modify hydraulic conductivity (no e-folding, no layers)
             * @param New conductivity
             */
            void setEffectivePorosity_direct(t_dim effectivePorosity) { set<t_dim, EffectivePorosity>(effectivePorosity); }

            /**
             * @brief Updates GW recharge
             * Currently assumes only one recharge as external flow!
             * @param amount The new flow amount
             * @param Should the recharge in the dynamic rivers be locked or updated by this change?
             */
            void updateUniqueFlow(double amount, FlowType flow = RECHARGE, bool lock = true) {
                if (lock and flow == RECHARGE) {
                    if (hasTypeOfExternalFlow(RIVER_MM)) {
                        //get current recharge and lock it bevor setting new recharge
                        //in arid regions recharge might be 0
                        t_vol_t recharge{0 * si::cubic_meter /day};
                        if(hasTypeOfExternalFlow(RECHARGE)){recharge = getExternalFlowByName(RECHARGE).getRecharge();}
                        getExternalFlowByName(RIVER_MM).setLock();
                        getExternalFlowByName(RIVER_MM).setLockRecharge(recharge);
                        //also lock conductance value
                        getExternalFlowByName(RIVER_MM).getERC(recharge,get<t_meter, EQHead>(),get<t_meter, Head>(),getEqFlow());
                    }
                }
                if (hasTypeOfExternalFlow(flow)) {
                    removeExternalFlow(flow);
                }

                addExternalFlow(flow, 0 * si::meter, amount, 0 * si::meter);
                if(numOfExternalFlows != externalFlows.size()){
                    throw "Number of external flows don't match";
                }
            }


            /**
             * Scale dynamic rivers for sensitivity
             * @param mult
             */
            void scaleDynamicRivers(double mult) {
                if (hasTypeOfExternalFlow(RIVER_MM)) {
                    getExternalFlowByName(RIVER_MM).setMult(mult);
                }
                return;
            }

            /**
             * @brief Update wetlands, lakes
             * @param amount
             * @param type
             */
            void updateExternalFlowConduct(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = getExternalFlowByName(type).getFlowHead();
                    double conduct = getExternalFlowByName(type).getConductance().value() * amount;
                    t_meter bottom = getExternalFlowByName(type).getBottom();
                    removeExternalFlow(type);
                    addExternalFlow(type, flowHead, conduct, bottom);
                }
            }

            /**
            * @brief Multiplies flow head for Sensitivity An. wetlands, lakes, rivers
            * @param amount
            * @param type
            */
            void updateExternalFlowFlowHead(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = getExternalFlowByName(type).getFlowHead() * amount;
                    double conduct = getExternalFlowByName(type).getConductance().value();
                    t_meter bottom = getExternalFlowByName(type).getBottom();
                    removeExternalFlow(type);
                    addExternalFlow(type, flowHead, conduct, bottom);
                }
            }

            /**
            * @brief Sets flowHead An. wetlands, lakes, rivers
            * @param amount
            * @param type
            */
            void setExternalFlowFlowHead(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = amount * si::meter;
                    double conduct = getExternalFlowByName(type).getConductance().value();
                    t_meter bottom = getExternalFlowByName(type).getBottom();
                    removeExternalFlow(type);
                    addExternalFlow(type, flowHead, conduct, bottom);
                }
            }

            /**
           * @brief adds delta to flowHead An. wetlands, lakes, rivers
           * @note Also checks for locked recharge
           * @param amount
           * @param type
           */
            void addExternalFlowFlowHead(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    ExternalFlow& externalFlow = getExternalFlowByName(type);
                    t_meter delta{amount * si::meter};
                    t_meter bottom{externalFlow.getBottom()};
                    t_meter flowHead{externalFlow.getFlowHead() + delta};
                    //The river is dry
                    if(std::isnan(amount)){ flowHead = bottom; }
                    double conduct{externalFlow.getConductance().value()};
                    bool lock{externalFlow.getLock()};
                    t_vol_t recharge{externalFlow.getLockRecharge()};
                    t_s_meter_t l_cond{externalFlow.getLockConduct()};
                    removeExternalFlow(type);
                    NANChecker(flowHead.value(), "Stage value");
                    NANChecker(l_cond.value(), "Conduct value");
                    NANChecker(bottom.value(), "Bottom value");

                    addExternalFlow(type, flowHead, conduct, bottom);
                    if (lock) {
                        getExternalFlowByName(type).setLock();
                        getExternalFlowByName(type).setLockRecharge(recharge);
                        getExternalFlowByName(type).setLockConduct(l_cond);
                    }
                }
            }

            /**
             * @brief Update lake bottoms
             * Used for sensitivity
             * @param amount
             */
            void updateLakeBottoms(double amount) {
                if (hasTypeOfExternalFlow(LAKE)) {
                    t_meter flowHead = getExternalFlowByName(LAKE).getFlowHead();
                    double conduct = getExternalFlowByName(LAKE).getConductance().value();
                    t_meter bottom = getExternalFlowByName(LAKE).getBottom() * amount;
                    removeExternalFlow(LAKE);
                    addExternalFlow(LAKE, flowHead, conduct, bottom);
                }
            }

            /**
             * @brief Check for type river
             * @return bool
             */
            bool hasRiver() { return hasTypeOfExternalFlow(RIVER); }

            /**
             * @brief Check for type GHB
             * @return bool
             */
            bool hasGHB() { return hasTypeOfExternalFlow(GENERAL_HEAD_BOUNDARY); }


            /**
             * @brief Get Q part (external sources) of flow equations
             * @return volume over time
             */
            t_vol_t getQ() noexcept { // Question: Rename to getExternalSources() ?
                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me cell has no special_flow
                }
                t_dim slope = get<t_dim, Slope>();
                t_vol_t eqFlow = getEqFlow();
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                for (const auto &flow : externalFlows) {
                    out += flow.second.getQ(eq_head, head, recharge, slope, eqFlow) * get<t_dim, StepModifier>();
                }
                return out;
            }

            t_s_meter_t getEffectivePorosityTerm(){ // todo: debugging
                t_s_meter_t out = (get<t_dim, EffectivePorosity>() *
                                   (get<t_meter, EdgeLengthLeftRight>() * get<t_meter, EdgeLengthFrontBack>())) /
                                  (day * get<t_dim, StepModifier>());
                return out;
            }

            t_vol_t getZetaMovementRHS(int zetaID){ // todo: debugging, call it from a solver (in MF: each layer and in each layer each zeta surface are solved individually!)
                //
                t_vol_t porosityTerm = getEffectivePorosityTerm() * Zetas[zetaID];
                t_vol_t sourceTermBelowZeta = getSourceTermBelowZeta(zetaID); // in MF: with SWIHCOF and BRHS of current layer and zeta
                t_vol_t pseudoSource_Zetas = getPseudoSource_Zetas(zetaID);
                t_vol_t out = porosityTerm - sourceTermBelowZeta + pseudoSource_Zetas;
                return out;
            }

            /**
             * @brief The matrix entry for the zeta surface
             * @return map <CellID,Conductance>
             * The left hand side of the equation
             */
            std::unordered_map<large_num, t_s_meter_t>  getConductance_ZetaMovement(int zetaID){ // Question: Rename to getLeftHandSide_Zeta?
                // todo: potential to enhance computational speed: change zoneConductanceCum to compute a single parameter instead of a vector
                // todo: debugging
                size_t numC = 5;
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(numC);

                DensityProperties densityProps = get<DensityProperties, densityProperties>();
                vector<t_dim> eps = densityProps.getEps();
                vector<t_dim> delnus = densityProps.getDelnus();
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};

                // pseudo source term calculation
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) {
                        //No neighbouring node
                        continue;
                    } else {
                        //There is a neighbour node

                        t_meter edgeLength_neig = 0 * si::meter;
                        t_meter edgeLength_self = 0 * si::meter;

                        edgeLength_self = getEdgeLengthSelf(got);
                        edgeLength_neig = getEdgeLengthNeig(got);

                        vector<t_meter> zetasNeig = at(got)->Zetas;
                        std::vector<t_meter> zoneThicknesses = vdf.calculateZoneThicknesses(
                                Zetas, zetasNeig, edgeLength_neig, edgeLength_self);
                        std::vector<t_s_meter_t> zoneConductances =
                                vdf.calculateDensityZoneConductances(
                                        zoneThicknesses,
                                        mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got)));
                        std::vector<t_s_meter_t> zoneConductancesCum =
                                vdf.calculateCumulativeDensityZoneConductances(zoneConductances);

                        t_s_meter_t zetaMovementConductance =
                                vdf.calculateZetaMovementConductances( zetaID,zoneConductances,
                                                                       zoneConductancesCum,delnus,eps);

                        NANChecker(zetaMovementConductance.value(), "zetaMovementConductance");
                        out[nodes->at(got->second)->get<large_num, ID>()] = move(zetaMovementConductance);
                    }
                }

                t_s_meter_t tmp_c = 0 * (si::square_meter / day);
                for (const auto &storedZetaMovementCond : out) { tmp_c = tmp_c - storedZetaMovementCond.second; }

                t_s_meter_t effectivePorosityTerm = getEffectivePorosityTerm();
                tmp_c = tmp_c + effectivePorosityTerm;
                NANChecker(effectivePorosityTerm.value(), "effectivePorosityTerm");
                out[get<large_num, ID>()] = tmp_c;
                return out;
            }

            t_vol_t getVerticalLeakage(){ // todo debugging
                t_vol_t out = 0 * (si::cubic_meter / day);
                t_s_meter_t verticalConductance;
                t_vol_t verticalFlux;
                t_vol_t verticalFluxCorrected;

                t_meter head_dif = 0 * si::meter;
                for (int n = 0; n <= zones.size() - 1; n++){
                    head_dif -= zones[n] * (Zetas[n+1] - Zetas[n]);
                }

                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::TOP, NeighbourPosition::DOWN};
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) {//No top node
                    } else {//Current node has a top node
                        verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                        //vector<t_meter> zetas_neig = at(got)->Zetas;
                        vector<t_dim> zonesTopNode = at(got)->zones;
                        t_dim nusBottomOfTopNode = zonesTopNode.back(); // in SWI2: NUBOT_i,j,k-1
                        t_dim nusTopOfThisNode = zones.front();  // in SWI2: NUTOP_i,j,k
                        verticalFlux = (verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(got)));
                        verticalFluxCorrected = verticalFlux - verticalFluxCorrection({position}); // (head_dif + (0.5 * (zetas_neig.back() + Zetas.front())) * (nusBottomOfTopNode + nusTopOfThisNode));
                        if (got->first == NeighbourPosition::TOP){
                            out += verticalFluxCorrected;
                        } else { // NeighbourPosition::DOWN
                            out -= verticalFluxCorrected;
                        }
                        // todo implement:
                        //  nubelbot = NUTOP(j,i,k+1)
                        //  IF ((qzbot.LT.0).AND.(nubelbot.LT.NUS(iz)).AND.(NUTOP(j,i,k).LE.nubelbot)) THEN
                        //      BRHS(j,i,iz) = BRHS(j,i,iz) + 0 // already implemented: ELSE BRHS(j,i,iz) = BRHS(j,i,iz) + qzbot
                        // and
                        //  nuontop = NUBOT(j,i,k-1)
                        //  IF ((qztop.LT.0).AND.(nuontop.GE.NUS(iz)).AND.NUBOT(j,i,k).GE.nuontop)) THEN
                        //      BRHS(j,i,iz) = BRHS(j,i,iz) + qztop

                    }
                }
                return out;
            }


            /**
             * The source term below a zeta surface (in SWI2: G)
             * return volume per time
             */
            t_vol_t getSourceTermBelowZeta(int zetaID){ // in SWI2: G todo debugging
                //
                // get RHS of PREVIOUS time step (without VDF terms pseudo source term and flux correction)
                // todo make sure this is from previous time step!!
                t_vol_t rhs_old = get<t_vol_t, RHSConstantDensity_TZero>(); // in SWI2 code: RHSPRESWI // Equation.getRHS() returns a Matrix<double, Dynamic, 1>

                // get HCOF of NEW time step
                t_s_meter_t hcof = mechanics.getHCOF(steadyState,
                                                     get<t_dim, StepModifier>(),
                                                     getStorageCapacity(),
                                                     getP());

                // get vertical leakage of NEW time step (the only term of G that depends on zones of density
                t_vol_t verticalLeakageTerm = getVerticalLeakage();

                // G = RHS(before VDF) - HCOF_(i,j,k,n)*h^(m)_(i,j,k) + (verticalLeakage_(i,j,k-1,n) - verticalLeakage_(i,j,k,n))
                t_vol_t out = rhs_old - hcof * get<t_meter, Head>() + verticalLeakageTerm; // in SWI2 code, (rhs_old - hcof * get<t_meter, Head>())) is multiplied by "fact = thickb / thick" that depends on IZONENR
                return out;
            }

            t_vol_t getPseudoSource_Zetas(int zetaID) { // for one zone
                // todo debug and add notes/description
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                DensityProperties densityProps = get<DensityProperties, densityProperties>();
                vector<t_dim> eps = densityProps.getEps();
                vector<t_dim> delnus = densityProps.getDelnus();
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
                t_meter edgeLength_neig;
                t_meter edgeLength_self;

                // pseudo source term calculation
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) {
                        continue;
                    }
                    edgeLength_neig = getEdgeLengthNeig(got);
                    edgeLength_self = getEdgeLengthSelf(got);

                    vector<t_meter> zetasNeig = at(got)->Zetas;
                    std::vector<t_meter> zoneThicknesses = vdf.calculateZoneThicknesses(
                            Zetas, zetasNeig, edgeLength_neig, edgeLength_self);
                    std::vector<t_s_meter_t> zoneConductances =
                            vdf.calculateDensityZoneConductances(
                                    zoneThicknesses,
                                    mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got)));
                    std::vector<t_s_meter_t> zoneConductancesCum =
                            vdf.calculateCumulativeDensityZoneConductances(zoneConductances);

                    for (int i = 0; i <= Zetas.size() - 2; i++) {
                        // 1st part: cumulative zone conductances and the difference between the new heads
                        out += zoneConductancesCum[i] * (getAt<t_meter, Head>(got) - get<t_meter, Head>());
                        // 2nd part: density variation in this zone (eps), zone conductance & "old" zeta surface heights
                        if (eps[i] != 0 and i == zetaID) {
                            out += eps[i] *
                                   (zoneConductances[i] * ((Zetas[i] - zetasNeig[i + 1]) - (Zetas[i] - Zetas[i + 1])));
                        }
                        // 3rd part: density variation (delnus) in all zones except this one, cumulative zone conductance and "old" zeta surface heights
                        if (delnus[i] != 0 and i != zetaID) {
                            out -= delnus[i] * (zoneConductancesCum[i] * ((zetasNeig[i] - Zetas[i])));
                        }
                    }
                }
                return out;
            }


            /**
                 * The pseudo source term for the flow equation, only used if variable density flow is active:
                 * This accounts for the effects of variable density flow
                 * @return volume per time
                 */
            t_vol_t getPseudoSource_Flow() {
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                DensityProperties densityProps = get<DensityProperties, densityProperties>();
                vector<t_dim> eps = densityProps.getEps();
                vector<t_dim> delnus = densityProps.getDelnus();
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};

                // pseudo source term calculation
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) {
                        continue;
                    }
                    t_meter edgeLength_neig = 0 * si::meter;
                    t_meter edgeLength_self = 0 * si::meter;

                    edgeLength_neig = getEdgeLengthNeig(got);
                    edgeLength_self = getEdgeLengthSelf(got);

                    vector<t_meter> zetasNeig = at(got)->Zetas;
                    std::vector<t_meter> zoneThicknesses = vdf.calculateZoneThicknesses(
                            Zetas, zetasNeig, edgeLength_neig, edgeLength_self);
                    std::vector<t_s_meter_t> zoneConductances =
                            vdf.calculateDensityZoneConductances(
                                    zoneThicknesses,
                                    mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got)));
                    std::vector<t_s_meter_t> zoneConductancesCum =
                            vdf.calculateCumulativeDensityZoneConductances(zoneConductances);

                    for (int i = 0; i <= Zetas.size() - 2; i++){
                        if (eps[i] != 0) {
                            out += eps[i] * (zoneConductances[i] * ((Zetas[i] - zetasNeig[i+1]) - (Zetas[i] - Zetas[i+1])));
                        }
                        if (delnus[i] != 0) {
                            out -= delnus[i] * (zoneConductancesCum[i] * ((zetasNeig[i] - Zetas[i])));
                        }
                    }
                }
                return out;
            }

            /**
             * @brief calculates the flux correction term in vertical direction
             * (in SWI2 documentation: BOUY; in code: QLEXTRA)
             * @out meter
             */
            t_meter verticalFluxCorrectionTerm(){
                DensityProperties densityProps = get<DensityProperties, densityProperties>();
                vector<t_dim> delnus = densityProps.getDelnus();
                t_meter out = 0 * si::meter;
                // dimensionless density at the top of current node (in SWI2: NUTOP)
                t_dim nusTopOfThisNode = zones.front();  // in SWI2: NUTOP_i,j,k
                for(int i = 1; i <= zones.size(); i++) {
                    if (Zetas[i] == get<t_meter, Elevation>()) { // if there is a ZETA surface at the top of current node (in SWI2: IPLPOS_(i,j,k,n) = 1)
                        nusTopOfThisNode += delnus[i];
                    }
                }
                // find the top neighbor
                std::unordered_map<NeighbourPosition, large_num>::const_iterator got =
                        neighbours.find(NeighbourPosition::TOP);
                if (got == neighbours.end()) {//No top node
                } else {//Current node has a top node
                    vector<t_dim> zonesTopNode = at(got)->zones;
                    vector<t_meter> zetasTopNode = at(got)->Zetas;
                    // dimensionless density at the bottom of the top node (in SWI2: NUBOT)
                    t_dim nusBottomOfTopNode = zonesTopNode.back(); // in SWI2: NUBOT_i,j,k-1
                    for(int i = 1; i <= delnus.size(); i++) {
                        t_meter bottomOfTopNode = getAt<t_meter, Elevation>(got) - getAt<t_meter, VerticalSize>(got);
                        t_meter headTopNode = getAt<t_meter, Head>(got);
                        if (zetasTopNode[i] == bottomOfTopNode or // if there is a ZETA surface at the bottom of the top node (in SWI2: IPLPOS_(i,j,k-1,n) = 2)
                            headTopNode < bottomOfTopNode){ // or the hydraulic head is below the bottom of the node
                            nusBottomOfTopNode -= delnus[i];
                        }
                    }
                    // first part of the flux correction term
                    for (int i = 0; i <= zonesTopNode.size(); i++){
                        out -= zonesTopNode[i] * (zetasTopNode[i] - zetasTopNode[i+1]); // Question: how to deal with this: in documentation is, BOUY is calculated with the simple sum (would be out +=), MODFLOW code for headdiff is as implemented (like out -=)
                    }
                    // second part of the flux correction term
                    out += 0.5 * (zetasTopNode.back() - Zetas.front()) * (nusBottomOfTopNode + nusTopOfThisNode); // Question: how to deal with this: in documentation, BOUY is calculated with a - between NUBOT and NUTOP, but in the code there is a - in the calculation of QLEXTRA
                }
                return out;
            }

            /**
             * @brief Calculates the vertical flux correction for variable density flow (in SWI2 documentation: BOUY; in SWI2 code: QLEXTRA)
             * @param possible_neighbours neighbours (TOP and/or DOWN) to/from which flow should be corrected
             * @return volume per time
             */
            t_vol_t verticalFluxCorrection(std::forward_list<NeighbourPosition> possible_neighbours) noexcept { // todo test for case where top and down cells exist
                t_vol_t out = 0 * (si::cubic_meter / day);
                t_meter fluxCorrectionTerm;
                t_s_meter_t verticalConductance;

                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) { // no neighbour at position
                    } else if (got->first == NeighbourPosition::DOWN) {
                        fluxCorrectionTerm = verticalFluxCorrectionTerm(); // flux correction term for this node
                        verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                        out -= fluxCorrectionTerm * verticalConductance; // -=
                    } else { // got->first == NeighbourPosition::TOP
                        fluxCorrectionTerm = at(got)->verticalFluxCorrectionTerm(); // flux correction term for top node
                        verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                        out += fluxCorrectionTerm * verticalConductance; // +=
                    }
                }
                return out;
            }

            /**
             * @brief Updating zeta surface heights at the top after a solution is found for the groundwater heads.
             * Zeta surface heights at the top are set to the new groundwater heads.
             */
            void updateZetaSurfaceHeights(){
                t_meter bottomOfNode = get<t_meter, Elevation>() - get<t_meter, VerticalSize>();
                t_meter topOfNode = get<t_meter, Elevation>();
                t_meter updatedZeta;
                // if groundwater head is BELOW the top of the node
                if (get<t_meter, Head>() < topOfNode){
                    // if groundwater head is ABOVE the bottom of the node
                    if (get<t_meter, Head>() > bottomOfNode){
                        updatedZeta = get<t_meter, Head>();
                        // if groundwater head is BELOW OR EQUAL to the bottom of the node
                    } else { // if bottomOfNode <= get<t_meter, Head>()
                        updatedZeta = bottomOfNode;
                    }
                    // update the first zeta surface
                    setZeta(0, updatedZeta);
                    // update all other zeta surfaces that are ABOVE the updated first zeta surface
                    for (int n = 1; n < Zetas.size() - 1; n++){
                        if (Zetas[n] > updatedZeta) { setZeta(n, updatedZeta); }
                    }
                    // if groundwater head is ABOVE the top of the node
                } else { // get<t_meter, Head>() > topOfNode
                    // clip zeta to the top of the node
                    setZeta(0, topOfNode);
                }
            }

            /*
             * @brief adjust zeta surface heights after the flow and zeta equations have found solutions
             */
            void adjustZetaHeights(){
                verticalZetaMovement();
                horizontalZetaMovement();
                clipZetaHeights();
                correctCrossingZetas();
                preventZetaLocking();
            }


            void verticalZetaMovement() { // todo: call this, debug

                // skip nodes where head is below the bottom of the node
                t_meter bottom_self = get<t_meter, Elevation>() - get<t_meter, VerticalSize>();
                if (get<t_meter, Head>() < bottom_self) {
                    t_vol_t verticalFluxTop; // in SWI2: qztop
                    t_s_meter_t verticalConductanceTop;
                    t_meter deltaZeta;
                    map_itter top = neighbours.find(NeighbourPosition::TOP);
                    // skip nodes that do not have a top neighbour
                    if (top == neighbours.end()) { // no neighbour at position
                    } else {
                        // calculate flux through the top
                        verticalConductanceTop = mechanics.calculateVerticalConductance(createDataTuple(top));
                        verticalFluxTop = verticalConductanceTop * (get<t_meter, Head>() - getAt<t_meter, Head>(top)) -
                                          verticalFluxCorrection({NeighbourPosition::TOP});

                        for (int n = 0; n < Zetas.size() - 1; n++) {
                            // zeta only moves through a node surface if if there is a ZETA surface
                            // - at the top of current node (in SWI2: IPLPOS_(i,j,k,n) = 1)
                            // - and at the bottom of the top node (in SWI2: IPLPOS_(i,j,k-1,n) = 2)
                            t_meter bottomOfTopNode = getAt<t_meter, Elevation>(top) - getAt<t_meter, VerticalSize>(top);
                            if (Zetas[n] >= get<t_meter, Elevation>() and
                                at(top)->Zetas[n] >= bottomOfTopNode) {
                                // if vertical flux through the top of the node is positive...
                                if (verticalFluxTop > (0 * si::cubic_meter / day)) {
                                    deltaZeta = (verticalFluxTop * (day * get<t_dim, StepModifier>())) /
                                                (get<t_s_meter, Area>() * getAt<t_dim, EffectivePorosity>(top));
                                    // ...lift zeta height of the lowest zeta surface in this node
                                    at(top)->setZeta(n, at(top)->Zetas.back() + deltaZeta);
                                // if vertical flux through the top of the node is negative...
                                } else if (verticalFluxTop < (0 * si::cubic_meter / day)) {
                                    deltaZeta = (verticalFluxTop * (day * get<t_dim, StepModifier>())) /
                                                (get<t_s_meter, Area>() * get<t_dim, EffectivePorosity>());
                                    // ...lower zeta height of the highest zeta surface in this node
                                    setZeta(n, Zetas.back() + deltaZeta);
                                }
                            }
                        }
                    }
                }
            }


            /**
             * @brief tracking the tips and toes of the zeta surfaces (at top and bottom of the current node)
             */
            void horizontalZetaMovement(){ // in SWI2 code: SSWI2_HORZMOVE, in documentation: "Tip and Toe Tracking"
                // todo do not change tips and toes at boundary nodes
                // todo call this, debug
                vector<t_meter> zetas_neig; // zetas of the neighbor in current direction
                vector<t_meter> zetas_neig2; // zetas of the neighbor in the opposite direction
                std::unordered_map<NeighbourPosition, NeighbourPosition> oppositeDirections;
                oppositeDirections[NeighbourPosition::BACK] = NeighbourPosition::FRONT;
                oppositeDirections[NeighbourPosition::FRONT] = NeighbourPosition::BACK;
                oppositeDirections[NeighbourPosition::LEFT] = NeighbourPosition::RIGHT;
                oppositeDirections[NeighbourPosition::RIGHT] = NeighbourPosition::LEFT;
                DensityProperties densityProps = get<DensityProperties, densityProperties>();
                t_dim maxToeSlope = densityProps.getMaxToeSlope();
                t_dim maxTipSlope = densityProps.getMaxTipSlope();

                t_meter edgeLength_self;
                t_meter edgeLength_neig;
                t_meter edgeLength_neig2; // edge length of opposite neighbour node
                t_meter maxDeltaZeta; // maximum height difference of a zeta surface n between adjacent nodes
                t_meter referenceElevation; // top/bottom of the node for tip/toe tracking
                t_dim slopeAdjustmentFraction = 0.1 * si::si_dimensionless; // slope adjustment fraction (in SWI2: ALPHA)
                t_dim minDepthThreshold = 0.1 * si::si_dimensionless; // minimum depth threshold for zeta surface toes (in SWI2: BETA)
                t_dim effPor_self; // effective porosity of this node
                t_dim effPor_neig; // effective porosity of neighbouring node
                t_dim effPor_neig2; // effective porosity of opposite neighbouring node
                t_meter zetaChange_self; // zeta height
                t_meter zetaChange_neig;
                t_meter top_self = get<t_meter, Elevation>();
                t_meter bottom_self = get<t_meter, Elevation>() - get<t_meter, VerticalSize>();
                t_meter top_neig;
                t_meter bottom_neig;

                std::unordered_map<string, t_dim> trackMap;
                trackMap["Toe"] = maxToeSlope;
                trackMap["Tip"] = maxTipSlope;
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
                for (const auto &tracker : trackMap) {
                    for (const auto &position : possible_neighbours) {
                        map_itter got = neighbours.find(position);
                        if (got == neighbours.end()) { // no neighbour at position
                        } else {
                            zetas_neig = at(got)->Zetas;
                            bottom_neig = getAt<t_meter, Elevation>(got) - getAt<t_meter, VerticalSize>(got);
                            top_neig = getAt<t_meter, Elevation>(got);

                            for (int n = 0; n < Zetas.size() -1; n++){

                                if (Zetas[n] > bottom_self and Zetas[n] < top_self) {
                                    // get reference elevation
                                    if (tracker.first == "Toe" and zetas_neig[n] == bottom_neig) {
                                        referenceElevation = bottom_self; // node bottom
                                    } else if (tracker.first == "Tip" and zetas_neig[n] == top_neig) {
                                        referenceElevation = top_self; // node top
                                    } else {
                                        continue;
                                    }

                                    // get length of the edges
                                    edgeLength_self = getEdgeLengthSelf(got);
                                    edgeLength_neig = getEdgeLengthNeig(got);

                                    // get max delta of zeta between nodes
                                    maxDeltaZeta = 0.5 * (edgeLength_self + edgeLength_neig) * tracker.second;

                                    // raise/lower zeta surfaces
                                    if (Zetas[n] - referenceElevation > maxDeltaZeta) {
                                        effPor_self = get<t_dim, EffectivePorosity>();
                                        effPor_neig = getAt<t_dim, EffectivePorosity>(got);
                                        // if tracking tip/toe: raise/lower this zeta surface in this node (ZETA_(i,j,k,n))
                                        zetaChange_self = slopeAdjustmentFraction * maxDeltaZeta *
                                                          ((effPor_neig * edgeLength_neig) /
                                                           ((effPor_self * edgeLength_self) + (effPor_neig * edgeLength_neig)));
                                        // for tracking tip/toe: lower/raise this zeta surface in neighbouring node (e.g. ZETA_(i,j+1,k,n))
                                        zetaChange_neig = slopeAdjustmentFraction * maxDeltaZeta *
                                                          ((effPor_self * edgeLength_self) /
                                                           ((effPor_self * edgeLength_self) + (effPor_neig * edgeLength_neig)));
                                    }


                                    if (tracker.first == "Toe" and zetas_neig[n] == bottom_neig) {
                                        if((Zetas[n] - Zetas.back()) > maxDeltaZeta) {
                                            setZeta(n, Zetas[n] - zetaChange_self);
                                            setZeta(n, at(got)->Zetas[n] + zetaChange_neig);
                                        }
                                    } else if (tracker.first == "Tip" and zetas_neig[n] == top_neig){
                                        if((Zetas[n] - Zetas.back()) > maxDeltaZeta) {
                                            setZeta(n, Zetas[n] + zetaChange_self);
                                            setZeta(n, at(got)->Zetas[n] - zetaChange_neig);
                                        }
                                    } else {
                                        continue;
                                    }

                                    if ((Zetas[n] - Zetas.back()) < (minDepthThreshold * zetaChange_neig)){
                                        zetas_neig = at(got)->Zetas;
                                        if (zetas_neig[n] > bottom_self and zetas_neig[n] < top_self) {
                                            // change zeta in other direction neighbour
                                            map_itter got2 = neighbours.find(oppositeDirections[position]);
                                            zetas_neig2 = at(got2)->Zetas;
                                            edgeLength_neig2 = getEdgeLengthNeig(got2);
                                            effPor_neig2 = getAt<t_dim, EffectivePorosity>(got2);
                                            setZeta(n, zetas_neig2[n] + ((Zetas[n] - Zetas.back()) *
                                                                         (edgeLength_self * effPor_self) /
                                                                         (edgeLength_neig2 * effPor_neig2)));
                                            setZeta(n, zetas_neig[zetas_neig.size() - 1]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /**
             * @brief clip zeta surfaces that are outside of the bounds (bounds: first and last zeta surface)
             * in SWI2: SSWI2_ZETACLIP
             */
            void clipZetaHeights(){
                for (int n = 0; n < Zetas.size() - 1; n++){
                    if (Zetas[n] < Zetas.back()) {
                        setZeta(n, Zetas.back());
                    }
                    if (Zetas[n] > Zetas.front()) {
                        setZeta(n, Zetas.front());
                    }
                }
            }

            /**
             * @brief adjust zeta surfaces if they have crossed or their height difference is below a minimum value
             * in SWI2: SSWI2_ZETACROSS
             */
            void correctCrossingZetas(){
                // todo: debug with this info:
                //  n | n2      | n3
                //  1 | -       | -
                //  2 | 1       | 1,2,3
                //  3 | 2,1     | 2,3,4; 1,2,3,4
                //  4 | 3,2,1   | 3,4,5; 2,3,4,5; 1,2,3,4,5
                t_meter zetaDifferenceCap = 0.001 * si::meter;
                t_meter zetaAverage;
                t_meter zetaSum;
                int n2_min;
                int n2_max;
                for (int n = 1; n < Zetas.size() - 2; n++){
                    // if zeta surface n is very close to or lower than the zeta surface that SHOULD be below (n+1)
                    if (Zetas[n] - Zetas[n+1] < zetaDifferenceCap) {
                        // make the zeta height of both surfaces their average
                        zetaAverage = 0.5 * (Zetas[n] + Zetas[n+1]);
                        setZeta(n, zetaAverage);
                        setZeta(n+1, zetaAverage);
                        // if there are zeta surfaces above that could possibly be crossing
                        if (n >= 2){
                            for (int x = 1; x < n-1; x++){
                                // create a range from n_min (n-x) to n_max (n+1)
                                n2_min = n-x;
                                n2_max = n+1;
                                // if a zeta surface above (n-x) crosses or is very close to zeta surface below (n+1)
                                //  (which is now at the same, averaged, height of zeta surface n)
                                if(Zetas[n2_min] - Zetas[n2_max] < zetaDifferenceCap){
                                    // calculate the average height from that zeta surface above (n-x) until
                                    //  the zeta surface below (n+1)
                                    for(int n2 = n2_min; n2 <= n2_max; n2++){ zetaSum += Zetas[n2]; }
                                    zetaAverage = zetaSum / ((n2_max-n2_min) * si::si_dimensionless); // todo check whether that is true
                                    // set every zeta surface between n-1 and n+1 to the averaged value
                                    for(int n2 = n2_min; n2 <= n2_max; n2++){ setZeta(n2, zetaAverage); }
                                }
                            }
                        }
                    }
                }
            }

            /**
             * @brief
             * in SWI2: SSWI2_ANTILOCKMIN
             */
            void preventZetaLocking(){
                t_meter swiLock = 0.001 * si::meter;
                t_dim slopeAdjustmentFraction = 0.1 * si::si_dimensionless; // slope adjustment fraction (in SWI2: ALPHA)
                t_meter maxDeltaZeta;
                t_meter topOfNode_self;
                t_meter bottomOfNode_self;
                t_meter topOfNode_neig;
                t_meter bottomOfNode_neig;
                t_meter maxAdjustment_self;
                t_meter maxAdjustment_neig;
                t_meter deltaZeta_self;
                t_meter deltaZeta_neig;
                t_meter edgeLength_self;
                t_meter edgeLength_neig;
                t_dim effectivePorosity_self;
                t_dim effectivePorosity_neig;

                for (int n = 1; n < Zetas.size() - 1; n++){
                    topOfNode_self = get<t_meter, Elevation>();
                    bottomOfNode_self = get<t_meter, Elevation>() - get<t_meter, VerticalSize>();
                    if (Zetas[n] == topOfNode_self or Zetas[n] == bottomOfNode_self) {
                        // iterate through horizontal neighbours
                        std::forward_list<NeighbourPosition> possible_neighbours =
                                {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                                 NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
                        for (const auto &position: possible_neighbours) {
                            map_itter got = neighbours.find(position);
                            if (got == neighbours.end()) { //No horizontal neighbouring node
                            } else {
                                // determine max delta zeta
                                maxAdjustment_self = get<t_meter, VerticalSize>() * slopeAdjustmentFraction;
                                maxAdjustment_neig = getAt<t_meter, VerticalSize>(got) * slopeAdjustmentFraction;
                                if (swiLock > maxAdjustment_self or swiLock > maxAdjustment_neig) {
                                    maxDeltaZeta = std::min(maxAdjustment_self, maxAdjustment_neig);
                                } else {
                                    maxDeltaZeta = swiLock;
                                }

                                // calculate delta zeta of node and neighbour
                                edgeLength_self = getEdgeLengthSelf(got);
                                edgeLength_neig = getEdgeLengthNeig(got);
                                deltaZeta_self = maxDeltaZeta * (effectivePorosity_neig * edgeLength_neig) /
                                                 (effectivePorosity_self * edgeLength_self +
                                                  effectivePorosity_neig * edgeLength_neig);
                                deltaZeta_neig = maxDeltaZeta * (effectivePorosity_self * edgeLength_self) /
                                                 (effectivePorosity_self * edgeLength_self +
                                                  effectivePorosity_neig * edgeLength_neig);

                                topOfNode_neig = getAt<t_meter, Elevation>(got);
                                bottomOfNode_neig = getAt<t_meter, Elevation>(got) - getAt<t_meter, VerticalSize>(got);
                                // if a zeta surface is at the BOTTOM of this node and at the TOP of the neighbour
                                // else if a zeta surface is at the TOP of this node and at the BOTTOM of the neighbour
                                if (Zetas[n] == bottomOfNode_self and at(got)->Zetas[n] == topOfNode_neig) { // IPLPOS_neig = 1 and IPLPOS_self = 2
                                    // adjust zeta surface heights
                                    setZeta(n, Zetas[n] + deltaZeta_self);
                                    at(got)->setZeta(n, Zetas[n] - deltaZeta_neig);
                                } else if (Zetas[n] == topOfNode_self and at(got)->Zetas[n] == bottomOfNode_neig) { // Zetas[n] == topOfNode_self // IPLPOS_neig = 2 and IPLPOS_self = 1
                                    // adjust zeta surface heights
                                    setZeta(n, Zetas[n] - deltaZeta_self);
                                    at(got)->setZeta(n, Zetas[n] + deltaZeta_neig);
                                }
                            }
                        }
                    }
                }
            }
            // todo list:
            //  flow equation
            //  - UPDATE THE UPPER ZETA SURFACE TO THE TOP OF THE WATER TABLE FOR UNCONFINED CONDITIONS // SSWI2_UPZ1()
            //  - CHECK TIP AND TOE SLOPE FOR SWI2 PACKAGE // SSWI2_CHKSLOPE
            //  Someday:
            //  - INSTANTANEOUS MIXING TERMS FOR ZONE BUDGETS // SSWI2_IMIX(A)
            //  - CALCULATE PRE TIP TOE TRACKING CHANGE IN ZONE THICKNESS (ZONECHG1) // SSWI2_ZCHG(A)
            //  - CALCULATE POST TIP TOE TRACKING CHANGE IN ZONE THICKNESS (ZONECHG2) // SSWI2_ZCHG(A)

            t_meter getEdgeLengthNeig(map_itter got){
                if (got->first == NeighbourPosition::LEFT or got->first == NeighbourPosition::RIGHT){
                    return getAt<t_meter, EdgeLengthLeftRight>(got);
                } else { // NeighbourPosition::FRONT or NeighbourPosition::BACK
                    return getAt<t_meter, EdgeLengthFrontBack>(got);
                }
            }

            t_meter getEdgeLengthSelf(map_itter got){
                if (got->first == NeighbourPosition::LEFT or got->first == NeighbourPosition::RIGHT){
                    return get<t_meter, EdgeLengthLeftRight>();
                } else { // NeighbourPosition::FRONT or NeighbourPosition::BACK
                    return get<t_meter, EdgeLengthFrontBack>();
                }
            }

            /**
             * @brief Get P part of flow equations
             * @return volume over time
             */
            t_s_meter_t getP() noexcept {
                t_s_meter_t out = 0.0 * (si::square_meter / day);
                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me
                }
                t_dim slope = get<t_dim, Slope>();
                t_vol_t eqFlow = getEqFlow();
                for (const auto &flow : externalFlows) {
                    if (is(flow.second.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, WETLAND, GLOBAL_WETLAND)) {
                        if (flow.second.flowIsHeadDependant(get<t_meter, Head>())) {
                            out += flow.second.getP(eq_head, head, recharge, slope, eqFlow) * get<t_dim, StepModifier>();
                        }
                    } else {
                        out += flow.second.getP(eq_head, head, recharge, slope, eqFlow) * get<t_dim, StepModifier>();
                    }
                }
                return out;
            }

            /**
             * @brief Get flow which is not groundwater head dependent
             * @return volume over time
             * Flow can be added to constant flows on right side of the equations
             * If head is above river bottom for example
             */
            t_vol_t calculateNotHeadDependandFlows() noexcept {
                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me
                }
                t_dim slope = get<t_dim, Slope>();
                t_vol_t eqFlow = getEqFlow();
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                //Q part is already subtracted in RHS
                for (const auto &flow : externalFlows) {
                    if (is(flow.second.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, WETLAND, GLOBAL_WETLAND)) {
                        if (not flow.second.flowIsHeadDependant(get<t_meter, Head>())) {
                            out += flow.second.getP(eq_head, head, recharge, slope, eqFlow) * get<t_dim, StepModifier>() *
                                   flow.second.getBottom();
                        }
                    }
                }
                return out;
            }

            /**
             * @brief The jacobian entry for the cell (NWT approach)
             * @return map <CellID,Conductance>
             */
            std::unordered_map<large_num, t_s_meter_t> getJacobian() noexcept {
                std::unordered_map<large_num, t_s_meter_t> out;
                size_t numC = 7;
                out.reserve(numC);
                unordered_map<large_num, t_s_meter_t> map = getConductance();
                t_s_meter_t tmp_hcofCRVC = map[get<large_num, ID>()];  // Question: is this any map or from row above?
                double head_diff{0.0};
                for (auto &ele : map) {
                    map[ele.first] = ele.second * mechanics.getDerivate__NWT(
                            nodes->at(ele.first)->get<t_meter, Elevation>(),
                            nodes->at(ele.first)->get<t_meter, VerticalSize>(),
                            nodes->at(ele.first)->get<t_meter, Head>()
                    );
                    if (ele.first != get<large_num, ID>()) {
                        head_diff = nodes->at(ele.first)->get<t_meter, Head>().value()
                                    - get<t_meter, Head>().value();
                        out[ele.first] =
                                (ele.second.value() + ele.second.value() * mechanics.getDerivate__NWT(
                                        nodes->at(ele.first)->get<t_meter, Elevation>(),
                                        nodes->at(ele.first)->get<t_meter, VerticalSize>(),
                                        nodes->at(ele.first)->get<t_meter, Head>()
                                ) * head_diff) * si::square_meter / day;
                        out[get<large_num, ID>()] += map[ele.first].value() * head_diff * si::square_meter / day;;
                    }
                }
                out[get<large_num, ID>()] +=
                        ((mechanics.getHCOF(steadyState,
                                            get<t_dim, StepModifier>(),
                                            getStorageCapacity(),
                                            getP()).value() * mechanics.getDerivate__NWT(
                                get<t_meter, Elevation>(),
                                get<t_meter, VerticalSize>(),
                                get<t_meter, Head>()
                        ) * get<t_meter, Head>().value())
                         + tmp_hcofCRVC.value() - (getRHS().value() * mechanics.getDerivate__NWT(
                                get<t_meter, Elevation>(),
                                get<t_meter, VerticalSize>(),
                                get<t_meter, Head>()
                        ))) * si::square_meter / day;;
                return out;
            }

            /**
             * @brief The matrix entry for the cell
             * @return map <CellID,Conductance>
             * The left hand side of the equation
             */
            std::unordered_map<large_num, t_s_meter_t> getConductance() { // Question: rename to getLeftHandSide?
                size_t numC = 7;
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(numC);

                //Get all conductances from neighbouring cells
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::TOP, NeighbourPosition::BACK, NeighbourPosition::DOWN, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};

                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    t_s_meter_t conduct = 0; // Question: move this to after the else? And multiply 0 with the unit?
                    if (got == neighbours.end()) {
                        //No neighbouring node
                        continue;
                    } else {
                        //There is a neighbour node
                        if (got->first == TOP or got->first == DOWN) {
                            conduct = mechanics.calculateVerticalConductance(createDataTuple(got));
                        } else {
                            //TODO check for option!
                            //if (get<int, Layer>() > 0) {
                            //    conduct = mechanics.calculateEFoldingConductance(createDataTuple<Head>(got), get<t_meter, EFolding>(), getAt<t_meter, EFolding>(got));
                            //}else{
                            conduct = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got));
                            //}
                        }
                        NANChecker(conduct.value(), "Conductances");
                        /*if(conduct.value() == 0){
                            LOG(numerics) << "conductance to neighbour is 0";
                        }*/
                        out[nodes->at(got->second)->get<large_num, ID>()] = move(conduct);
                    }
                }

                t_s_meter_t tmp_c = 0; // Question: multiply with units?

                for (const auto &c : out) { tmp_c = tmp_c - c.second; }
                t_s_meter_t hcof = mechanics.getHCOF(steadyState,
                                                     get<t_dim, StepModifier>(),
                                                     getStorageCapacity(),
                                                     getP());
                tmp_c = tmp_c + hcof;
                NANChecker(tmp_c.value(), "HCOF");
                //if(tmp_c.value() == 0){
                //	LOG(numerics) << "HCOF term is 0";
                //}
                out[get<large_num, ID>()] = tmp_c;
                return out;
            };

            /**
             * @brief The right hand side of the equation
             * @return volume per time
             */
            t_vol_t getRHS() {
                t_vol_t externalFlows = -getQ();
                t_vol_t dewateredFlow = calculateDewateredFlow();
                t_vol_t rivers = calculateNotHeadDependandFlows();
                t_vol_t storageFlow =
                        getStorageCapacity() * (get<t_meter, Head_TZero>() / (day* get<t_dim, StepModifier>()));
                DensityProperties densityProps = get<DensityProperties, densityProperties>();
                if (steadyState) {
                    storageFlow = 0 * (si::cubic_meter / day);
                }

                t_vol_t out = externalFlows + dewateredFlow - rivers - storageFlow;

                if (densityProps.isDensityVariable()){ // todo maybe make condition: zones.size() > 1
                    // save RHS without variable density terms for calculation of zeta movement at next time step
                    set < t_vol_t, RHSConstantDensity_TZero > (out);
                    // calculate variable density terms
                    t_vol_t pseudoSourceFlow = getPseudoSource_Flow();
                    t_vol_t fluxCorrection = verticalFluxCorrection({NeighbourPosition::TOP, NeighbourPosition::DOWN});
                    // add variable density terms to RHS
                    out += pseudoSourceFlow + fluxCorrection;
                }

                NANChecker(out.value(), "RHS");
                return out;
            }

            /**
             * @brief The right hand side of the equation (NWT)
             * @return volume per time
             */
            double getRHS__NWT() noexcept {
                double out{0.0};
                const unordered_map<large_num, t_s_meter_t> &map = getConductance();
                const double sum = std::accumulate(std::begin(map), std::end(map), 0,
                                                   [this](const size_t previous,
                                                          const std::pair<int, t_s_meter_t>
                                                          &p) {
                                                       return previous + (p.second.value() *
                                                                          mechanics.smoothFunction__NWT(
                                                                                  nodes->at(
                                                                                          p.first)->get<t_meter, Elevation>(),
                                                                                  nodes->at(
                                                                                          p.first)->get<t_meter, VerticalSize>(),
                                                                                  nodes->at(
                                                                                          p.first)->get<t_meter, Head>()));
                                                   });
                out = sum - getRHS().value();
                return out;
            }

            void setHead(t_meter head) noexcept {
                __setHead(head);
            }

            t_meter calcInitialHead(t_meter initialParam) noexcept { return __calcInitialHead(initialParam); }

            bool isStaticNode() noexcept { return __isStaticNode(); }

            PhysicalProperties &getProperties() { return fields; } // Question: rename to getNodeProperties?

            void enableNWT() { nwt = true; }

            /**
             * @brief Calculate non storage related in and out flow
             * @return Flow volume
             */
            template<typename CompareFunction>
            t_vol_t getNonStorageFlow(CompareFunction compare) noexcept {
                t_vol_t out = 0.0 * si::cubic_meter / day;
                t_vol_t ex;
                for (const auto &flow : externalFlows) {
                    ex = calculateExternalFlowVolume(flow.second);
                    if (compare(ex.value())) {
                        out += boost::units::abs(ex);
                    }
                }
                t_vol_t dwateredFlow = -calculateDewateredFlow();
                if (compare(dwateredFlow.value())) {
                    out += boost::units::abs(dwateredFlow);
                }
                return out;
            }

            /**
             * Calculate the lateral flow velocity
             * @param pos
             * @return
             */
            quantity<Velocity> getVelocity(map_itter pos) {
                t_vol_t lateral_flow{0 * si::cubic_meter / day};
                t_meter vert_size = get<t_meter, VerticalSize>();
                t_s_meter_t conductance;
                //TODO check for option!
                //if (get<int, Layer>() > 0) {
                //    conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(pos), get<t_meter, EFolding>(), getAt<t_meter, EFolding>(pos));
                //}else{
                conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(pos));
                //}

                lateral_flow = conductance * (get<t_meter, Head>() - getAt<t_meter, Head>(pos));
                return lateral_flow / (vert_size * vert_size);
            }

            /**
             * @brief Calculate flow velocity for flow tracking
             * Vx and Vy represent the flow velocity in x and y direction.
             * A negative value represents a flow in the oposite direction.
             * @return Velocity vector (x,y)
             */
            std::pair<double, double> getVelocityVector() {
                quantity<Velocity> Vx{0};
                quantity<Velocity> Vy{0};

                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) {
                        continue;
                    }
                    if (got->first == NeighbourPosition::LEFT) {
                        Vx += -getVelocity(got);
                    }
                    if (got->first == NeighbourPosition::BACK) {
                        Vy += -getVelocity(got);
                    }

                    if (got->first == NeighbourPosition::RIGHT) {
                        Vx += getVelocity(got);
                    }
                    if (got->first == NeighbourPosition::FRONT) {
                        Vy += getVelocity(got);
                    }
                }

                return std::make_pair(Vx.value(), Vy.value());
            };

        };

/**
 * @class StandardNode
 * A standard groundwater node
 */
        class StandardNode : public NodeInterface {
        public:
            StandardNode(std::shared_ptr<std::vector<std::unique_ptr<NodeInterface>>> nodes,
                         double lat,
                         double lon,
                         t_s_meter area,
                         t_meter edgeLengthLeftRight,
                         t_meter edgeLengthFrontBack,
                         large_num ArcID,
                         large_num ID,
                         t_vel K,
                         int stepmodifier,
                         double aquiferDepth,
                         double anisotropy,
                         double specificYield,
                         double specificStorage,
                         bool confined,
                         DensityProperties densityProps)
                    : NodeInterface(nodes, lat, lon, area, edgeLengthLeftRight, edgeLengthFrontBack, ArcID, ID, K,
                                    stepmodifier, aquiferDepth, anisotropy, specificYield, specificStorage, confined, densityProps) {}

        private:
            // implementation
            friend class NodeInterface;

            //Learning weight
            t_dim weight = 0.1;

            /**
             * @brief Update heads after one or multiple inner iterations
             * @param head
             */
            virtual void __setHead(t_meter delta) {
                NANChecker(delta.value(), "Set Head");
                //t_meter deltaH__old = get<t_meter, HeadChange>();
                t_meter current_head = get<t_meter, Head>();
                //t_meter delta = head - current_head;
                set<t_meter, HeadChange>(delta);
                set<t_meter, Head>(current_head + delta);

                /*	if (nwt){
                        //TODO move to Numerics
                        //Underrelaxation with delta-bar-delta from Smith 1993:

                        t_dim gamma = 0 * si::si_dimensionless;
                        t_dim theta = 0.9;
                        t_dim kappa = 0.00001;
                        t_dim momentum = 0.1;
                        t_meter deltaH = (1 - gamma) * get<t_meter, HeadChange
                        >() + gamma * deltaH__old;
                        if (weight < 1)
                        {
                            weight = weight + kappa;
                        }
                        else
                        {
                            weight = weight - theta * weight;
                        }
                        t_meter h = head + weight * get<t_meter, HeadChange>() + momentum * deltaH;
                        set<t_meter, HeadChange>(h - get<t_meter, Head>());
                        set<t_meter, Head>(h);
                    }*/
            };

            virtual t_meter
            __calcInitialHead(t_meter initialParam) {
                t_meter elevation = get<t_meter, TopElevation>();
                if (elevation >= initialParam) {
                    return elevation - initialParam;
                }
                return elevation;
            }

            virtual bool
            __isStaticNode() { return false; }

            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version) {
                boost::serialization::void_cast_register<NodeInterface, StandardNode>();
                boost::serialization::void_cast_register<StandardNode, NodeInterface>();
                boost::serialization::base_object<NodeInterface>(*this);
                LOG(debug) << "Serializing Standard Node";
                ar & nodes;
                ar & neighbours;
                ar & externalFlows;
                ar & numOfExternalFlows;
                ar & nwt;
                ar & initial_head;
                ar & simpleDistance;
                ar & simpleK;
                ar & steadyState;
                ar & fields;
                ar & cached;
                ar & eq_flow;
            }

        };

/**
 * @class StaticHeadNode
 * A node without changing head
 * Can be used as boundary condition
 */
        class StaticHeadNode : public NodeInterface {
        public:
            StaticHeadNode(std::shared_ptr<std::vector<std::unique_ptr<NodeInterface>>> nodes,
                           large_num ID,
                           t_s_meter area,
                           t_meter edgeLengthLeftRight,
                           t_meter edgeLengthFrontBack,
                           DensityProperties densityProps) // Question: densityProps required here? If yes: what values?
                    : NodeInterface(
                    nodes,
                    0,
                    0,
                    area,
                    edgeLengthLeftRight,
                    edgeLengthFrontBack,
                    ID,
                    ID,
                    0.3 * (si::meter / day), 1, 100, 10, 0.15, 0.000015, true, densityProps) {}

        private:
            friend class NodeInterface;

            virtual void __setHead(t_meter head) {
                set<t_meter, HeadChange>(0 * si::meter);
            };

            virtual t_meter __calcInitialHead(t_meter initialParam) { return 0; }

            virtual bool __isStaticNode() { return true; }

            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version) {
                boost::serialization::void_cast_register<NodeInterface, StaticHeadNode>();
                boost::serialization::void_cast_register<StaticHeadNode, NodeInterface>();
                boost::serialization::base_object<NodeInterface>(*this);
                LOG(debug) << "Serializing Static Node";
                ar & nodes;
                ar & neighbours;
                ar & externalFlows;
                ar & numOfExternalFlows;
                ar & nwt;
                ar & initial_head;
                ar & simpleDistance;
                ar & simpleK;
                ar & steadyState;
                ar & fields;
                ar & cached;
                ar & eq_flow;
            }

        };

    }
}

#endif //NODE_HPP