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
 * Neighbouring positions for cells (side-view)
 *     TOP
 * LEFT * RIGHT
 *     DOWN
 *
 * Top view:
 *     FRONT
 * LEFT * RIGHT
 *     BACK
 *
 *  Top view with refinement into 4:
 *          FRONTLEFT FRONTRIGHT
 * LEFTFRONT      ******        RIGHTFRONT
 *                *node*
 * LEFTBACK       ******        RIGHTBACK
 *          BACKLEFT  BACKRIGHT
 *
 *  Top view with refinement into 9:
 *          FRONTLEFT FRONTFRONT FRONTRIGHT
 * LEFTFRONT           ******              RIGHTFRONT
 * LEFTLEFT            *node*              RIGHTRIGHT
 * LEFTBACK            ******              RIGHTBACK
 *          BACKLEFT  BACKBACK  BACKRIGHT
 *
 *  Top view with refinement into 16:
 *              FRONTLEFT FRONTFRONTLEFT FRONTFRONTRIGHT FRONTRIGHT
 * LEFTFRONT                        ******                          RIGHTFRONT
 * LEFTLEFTFRONT                    *node*                          RIGHTRIGHTFRONT
 * LEFTLEFTBACK                     ******                          RIGHTRIGHTBACK
 * LEFTBACK                         ******                          RIGHTBACK
 *              BACKLEFT  BACKBACKLEFT   BACKBACKRIGHT  BACKRIGHT
 *
 */
        enum NeighbourPosition {
            TOP = 1,
            DOWN,               // 2
            RIGHT,              // 3
            LEFT,               // 4
            FRONT,              // 5
            BACK,               // 6
            // for refined nodes:
            RIGHTFRONT,         // 7
            RIGHTBACK,          // 8
            LEFTFRONT,          // 9
            LEFTBACK,           // 10
            FRONTLEFT,          // 11
            FRONTRIGHT,         // 12
            BACKLEFT,           // 13
            BACKRIGHT,          // 14
            // for nodes refined into 9:
            RIGHTRIGHT,         // 15
            LEFTLEFT,           // 16
            FRONTFRONT,         // 17
            BACKBACK,           // 18
            // for nodes refined into 16:
            RIGHTRIGHTFRONT,    // 19
            RIGHTRIGHTBACK,     // 20
            LEFTLEFTFRONT,      // 21
            LEFTLEFTBACK,       // 22
            FRONTFRONTLEFT,     // 23
            FRONTFRONTRIGHT,    // 24
            BACKBACKLEFT,       // 25
            BACKBACKRIGHT       // 26
        };
    }
}

namespace std {
    template<>
    struct hash<GlobalFlow::Model::NeighbourPosition> {
        typedef GlobalFlow::Model::NeighbourPosition argument_type;
        typedef std::underlying_type<argument_type>::type underlying_type;

        std::size_t operator()(const argument_type &arg) const {
            std::hash<underlying_type> hasher;
            return hasher(static_cast< underlying_type >( arg ));
        }
    };
}

namespace GlobalFlow {
    namespace Model {

        using namespace boost::units;

        class NodeInterface;

        using p_node = std::unique_ptr<GlobalFlow::Model::NodeInterface>;
        using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

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
             * E-Folding function as defined by Ying Fan et al. e^(-Depth / E-Folding Depth).
             */
            t_dim efoldingFromData(t_meter depth) {
                t_meter folding = get<t_meter, EFolding>();
                if (folding == 0.0 * si::meter)
                    return 1 * si::si_dimensionless;
                //Alter if a different size should be used and not full vertical size
                //depth = (depth / (2 * si::si_dimensionless));
                t_dim out = exp(-depth / folding);
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
            std::unordered_map<NeighbourPosition, large_num> neighbours;
            std::unordered_map<NeighbourPosition, large_num> horizontal_neighbours;
            std::unordered_map<NeighbourPosition, large_num> vertical_neighbours;
            std::vector<NeighbourPosition> finer_neighbours;

            std::unordered_map<FlowType, ExternalFlow, FlowTypeHash> externalFlows;
            int numOfExternalFlows{0};
            bool initial_head{true};
            bool simpleDistance{false};

            friend class FluidMechanics;

            FluidMechanics mechanics;
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

            p_node & at(large_num nodeID) { return nodes->at(nodeID); }

            template<typename T, typename F>
            T getAt(large_num nodeID) {
                return at(nodeID)->get<T, F>();
            }

            /**
             * @brief Uses specific storage to calculate storativity
             * @return Squaremeter
             */
            t_s_meter getStorageCapacity__Primary() noexcept {
                t_s_meter out = get<quantity<perUnit>, SpecificStorage>().value()
                                    * get<t_c_meter, VolumeOfCell>().value() * si::square_meter;
                return out;
            }

            /**
             * @brief Uses specific yield to calculate storativity
             * @return Squaremeter
             */
            t_s_meter getStorageCapacity__Secondary() noexcept {
                return get<t_dim, SpecificYield>() * get<t_s_meter, Area>();
            }

            /**
             * @brief Flow volume IN or OUT of cell
             * @return Flow volume
             */
            template<typename CompareFunction>
            t_vol_t getFlow(CompareFunction compare) noexcept {
                t_vol_t out = 0.0 * si::cubic_meter / day;
                t_vol_t storageFlow = getStorageFlow();
                if (compare(storageFlow.value())) {
                    out += boost::units::abs(storageFlow);
                }
                out += getNonStorageFlow(std::forward<CompareFunction>(compare));
                return out;
            }

        public:
            virtual ~NodeInterface() = default;

            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version) {
                LOG(debug) << "Serializing abstract node";
                ar & nodes;
                ar & neighbours;
                ar & externalFlows;
                ar & numOfExternalFlows;
                ar & initial_head;
                ar & simpleDistance;
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
             * @param SpatID Unique Spat-ID
             * @param ID Internal ID = Position in vector
             * @param K Hydraulic conductivity in meter/day (default)
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
                          large_num SpatID,
                          large_num ID,
                          t_vel K,
                          t_meter head,
                          double aquiferDepth,
                          double anisotropy,
                          double specificYield,
                          double specificStorage,
                          bool useEfolding,
                          bool confined,
                          large_num refID,
                          large_num refinementInto,
                          bool isSteadyState,
                          bool isDensityVariable,
                          std::vector<t_dim> delnus,
                          std::vector<t_dim> nusInZones,
                          double effPorosity,
                          double maxTipSlope,
                          double maxToeSlope,
                          double minDepthFactor,
                          double slopeAdjFactor,
                          t_meter vdfLock,
                          int sourceZoneGHB,
                          int sourceZoneRecharge);

/*****************************************************************
Get Properties
******************************************************************/

            large_num getID() { return get<large_num, ID>(); }

            t_vel getK__pure() noexcept { return get<t_vel, K>(); }

            /**
             * @brief Get hydraulic conductivity
             * @return hydraulic conductivity (scaled by e-folding)
             */
            t_vel getK() noexcept {
                if (get<bool, UseEfolding>()){
                    t_dim e_fold{1 * si::si_dimensionless};
                    if (get<int, Layer>() > 0) {
                        e_fold = efoldingFromData(get<t_meter, VerticalSize>());
                    }
                    t_vel out = get<t_vel, K>() * e_fold * get<t_dim, StepSize>();
                    if (out < 1e-20 * si::meter / day) {
                        out = 1e-20 * si::meter / day;
                    }
                    return out;
                } else {
                    return get<t_vel, K>() * get<t_dim, StepSize>();
                }
            }

/*****************************************************************
Set Properties
******************************************************************/

            /**
             * @brief Set elevation on top layer and propagate to lower layers
             * @param elevation The top elevation (e.g. from DEM)
             */
            void setElevation_allLayers(t_meter elevation) {
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
            void setEqHead_allLayers(t_meter wtd) {
                t_meter eqhead = get<t_meter, Elevation>() - wtd;
                set < t_meter, EQHead > (eqhead);
                setHead(eqhead);
                applyToAllLayers([eqhead](NodeInterface *nodeInterface) {
                    try {
                        nodeInterface->setHead(eqhead);
                        nodeInterface->set<t_meter, EQHead>(eqhead);
                    }
                    catch (...) {}
                });
            }

            /**
             * @brief Update the current head change (in comparison to last time step)
             * @note Should only be called at end of time step
             */
            void updateHeadChange_TZero() noexcept {
                set < t_meter, HeadChange_TZero > (getHead() - getHead_TZero());
            }

            void updateHead_TZero() noexcept { setHead_TZero(getHead()); }

            void setHead_TZero(t_meter head) noexcept { set<t_meter, Head_TZero>(head); }

            void setHead_TZero_allLayers(t_meter head) noexcept {
                setHead_TZero(head);
                applyToAllLayers([&head](NodeInterface *nodeInterface) {
                    nodeInterface->setHead_TZero(head);
                });
            }



            void setHead(t_meter head) noexcept {
                NANChecker(head.value(), "Set Head");
                set<t_meter, Head>(head);
            }

            void setHead_allLayers(t_meter head) noexcept {
                setHead(head);
                applyToAllLayers([&head](NodeInterface *nodeInterface) {
                    nodeInterface->setHead(head);
                });
            }

/*****************************************************************
Helpers
******************************************************************/

            /**
             * Calculated equilibrium flow to neighbouring cells
             * Static thus calculated only once.
             *
             * Depends on: K in cell and eq_head in all 6 neighbours
             */
            bool cached{false};
            t_vol_t eq_flow{0 * si::cubic_meter / day};

            template<class HeadType>
            FlowInputHor createDataTuple(NeighbourPosition neigPos, large_num neigNodeID) {
                return std::make_tuple(at(neigNodeID)->getK(),
                                       getK(),
                                       at(neigNodeID)->getNodeLength(neigPos), // length of neighbour node (parallel to direction)
                                       getNodeLength(neigPos), // length of this node (parallel to direction)
                                       std::min(getNodeWidth(neigPos), at(neigNodeID)->getNodeWidth(neigPos)), // width of smaller node (perpendicular to direction)
                                       getAt<t_meter, HeadType>(neigNodeID),
                                       get<t_meter, HeadType>(),
                                       getAt<t_meter, Elevation>(neigNodeID),
                                       get<t_meter, Elevation>(),
                                       getAt<t_meter, VerticalSize>(neigNodeID),
                                       get<t_meter, VerticalSize>(),
                                       get<bool, Confinement>());
            }

            FlowInputVert createDataTuple(large_num neigNodeID) {
                return std::make_tuple(at(neigNodeID)->getK_vertical(),
                                       getK_vertical(),
                                       get<t_meter, VerticalSize>(),
                                       getAt<t_meter, VerticalSize>(neigNodeID),
                                       get<t_meter, Head>(),
                                       getAt<t_meter, Head>(neigNodeID),
                                       get<t_meter, Elevation>(),
                                       getAt<t_meter, Elevation>(neigNodeID),
                                       get<t_s_meter, Area>(),
                                       get<bool, Confinement>());
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

/*****************************************************************
Calculate
******************************************************************/

            /**
             * Calculate the lateral groundwater flow to the neighbouring nodes
             * Generic function used for calculating equilibrium and current step flow
             * @return
             */
            template<class HeadType>
            t_vol_t calcLateralFlows(bool onlyOut) {
                t_vol_t lateral_flow{0 * si::cubic_meter / day};
                t_s_meter_t conductance;

                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                        conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(neigPos, neigNodeID), get<t_meter, EFolding>(), getAt<t_meter, EFolding>(neigNodeID));
                    } else {
                        conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                    }

                    t_vol_t flow = conductance * (get<t_meter, HeadType>() - getAt<t_meter, HeadType>(neigNodeID));

                    if (onlyOut) {
                        if (flow.value() > 0) { lateral_flow -= flow; }
                    } else { lateral_flow -= flow; }
                }
                return lateral_flow;
            }

            /**
             * Calculate the lateral groundwater flow to (-)/from (+) a neighbouring node
             * @return
             */
            std::unordered_map< NeighbourPosition, double> getFlowToOrFromNeighbours() {
                std::unordered_map< NeighbourPosition, double> mapNeighboursToFlows;
                t_vol_t lateralFlow{0 * si::cubic_meter / day};
                t_s_meter_t conductance{0 * si::square_meter / day};
                std::vector<NeighbourPosition> neigPos_LRFB = getNeigPos_LRFB();
                for (const auto &neigPos: neigPos_LRFB) { // todo only LRFB or maybe also other nodes?
                    if (horizontal_neighbours.find(neigPos) == horizontal_neighbours.end()) {
                        mapNeighboursToFlows[neigPos] = std::nan("1") ;
                    } else {
                        large_num neigNodeID = horizontal_neighbours[neigPos];
                        if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                            conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(neigPos, neigNodeID),
                                                                                 get<t_meter, EFolding>(),
                                                                                 getAt<t_meter, EFolding>(neigNodeID));
                        } else {
                            conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                        }
                        // if this head is larger than neihgbour head, flow is out, which is negative
                        lateralFlow = -conductance * (getHead() - getAt<t_meter, Head>(neigNodeID)) *
                                      get<t_dim, StepSize>();
                        mapNeighboursToFlows[neigPos] = lateralFlow.value() ;
                    }


                }
                return mapNeighboursToFlows;
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
             * @return lateral flows
             */
            t_vol_t getLateralFlows() {
                return calcLateralFlows<Head>(false) * get<t_dim, StepSize>();
            }

            /**
             * Get the current lateral out flows
             * @return lateral outflows
             */
            t_vol_t getLateralOutFlows() {
                return calcLateralFlows<Head>(true) * get<t_dim, StepSize>();
            }

            /**
             * @brief Get hydraulic vertical conductivity
             * @return hydraulic conductivity scaled by anisotropy (scaled by e-folding)
             */
            t_vel getK_vertical() noexcept { return (getK() / get<t_dim, Anisotropy>()); }

            /**
             * @brief Modify hydraulic conductivity (applied to all layers below)
             * @param conduct conductivity (if e-folding enabled scaled on layers)
             */
            void setK_allLayers(t_vel conduct) {
                setK(conduct);
                applyToAllLayers([&conduct](NodeInterface *nodeInterface) {
                    nodeInterface->setK(conduct);
                });
            }

            /**
             * @brief Modify hydraulic conductivity (no e-folding, no layers)
             * @param conduct conductivity
             */
            void setK(t_vel conduct) { set < t_vel, K > (conduct); }

            int getSpatID() {return (int) get<large_num, SpatID>();}

            double getLat() {return get<double, Lat>();}

            double getLon() {return get<double, Lon>();}

            large_num getRefID() {return get<large_num, RefID>(); }

            t_s_meter getArea(){return get<t_s_meter, Area>();}

            t_meter getVerticalSize(){return get< t_meter, VerticalSize >(); }

            t_dim getEffectivePorosity(){return get<t_dim, EffectivePorosity>();}

            t_meter getEdgeLengthLeftRight(){return get<t_meter, EdgeLengthLeftRight>();}

            t_meter getEdgeLengthFrontBack(){return get<t_meter, EdgeLengthFrontBack>();}

            t_meter getElevation(){return get<t_meter, Elevation>();}

            t_meter getBottom(){return get<t_meter, Elevation>() - get<t_meter, VerticalSize>();}

            t_meter getHead(){ return get<t_meter, Head>(); }

            t_meter getHead_TZero() noexcept { return get<t_meter, Head_TZero>(); }

            int getLayer(){ return get<int, Layer>(); }

            large_num getRefinedInto() { return get<large_num, RefinedInto>(); }

            /**
             * @brief Get all outflow since simulation start
             */
            t_c_meter getOUT() noexcept { return get<t_c_meter, OUT>(); }

            /**
             * @brief Get all inflow since simulation start
             */
            t_c_meter getIN() noexcept { return get<t_c_meter, IN>(); }

            void updateStepSize(double stepSize) { set < t_dim, StepSize > (stepSize * si::si_dimensionless); }

            void updateIsSteadyState(bool isSteadyState) { set<bool, IsSteadyState>(isSteadyState); }

            void updateIsDensityVariable(bool isDensityVariable) { set<bool, IsDensityVariable>(isDensityVariable); }

            /**
             * @brief Storage capacity based on yield or specific storage
             * @return Potential flow budget when multiplied by head change
             * Uses an 0.001m epsilon to determine if a water-table condition is present.
             * If the layer is confined or not in water-table condition returns primary capacity.
             */
            t_s_meter getStorageCapacity() noexcept {
                if (get<bool, IsSteadyState>()) {
                    return 0 * si::square_meter;
                }

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
                //we are not in the first layer, and we are in unconfined conditions as specified by the user

                if (get<t_meter, Head>() + epsilon < get<t_meter, Elevation>()) {
                    //water-table condition
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
            ExternalFlow & getExternalFlowByName(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type);
                } else {
                    throw std::out_of_range("No such flow: " + std::to_string(type) + " at nodeID " + std::to_string(getID()));
                }
            }

            /**
             * @brief Get an external flow volume by its FlowType
             * @param type The flow type
             * @return Flow volume
             */
            t_vol_t getExternalFlowVolumeByName(FlowType type) {
                for (const auto & flow : externalFlows) {
                    if (type == flow.second.getType()) {
                        return calculateExternalFlowVolume(flow.second);
                    }
                }
                return 0 * si::cubic_meter / day;
            }

            /**
             * @brief Get conductance of an external flow by its FlowType
             * @param type The flow type
             * @return Conductance
             */
            double getExternalFlowConductance(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type).getConductance().value();
                } else {
                    return std::nan("1");
                }
            }

            /**
             * @brief Get head of an external flow by its FlowType
             * @param type The flow type
             * @return Head
             */
            double getExternalFlowElevation(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type).getFlowHead().value();
                } else {
                    return std::nan("1");
                }
            }

            /**
            * @brief Get bottom of an external flow by its FlowType
            * @param type The flow type
            * @return Bottom
            */
            t_meter getExternalFlowBottom(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type).getBottom();
                } else {
                    return 0 * si::meter;
                }
            }

            /**
             * @brief Get flow budget based on head change
             * @return Flow volume
             * Note: Water entering storage is treated as an outflow (-), that is a loss of water from the flow system
             * while water released from storage is treated as inflow (+), that is a source of water to the flow system
             */
            t_vol_t getStorageFlow() noexcept {
                return -getStorageCapacity() * get<t_meter, HeadChange_TZero>() / (day * get<t_dim, StepSize>());
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
                    return flow.getRecharge() * get<t_dim, StepSize>();
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
                t_vol_t eqFlow = getEqFlow();
                if (is(flow.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                    if (flow.isFlowHeadDependent(head)) {
                        ex = (flow.getP(eq_head, head, recharge, eqFlow) * head +
                              flow.getQ(eq_head, head, recharge, eqFlow)) * get<t_dim, StepSize>();
                    } else { // flow is not head dependent when the head is below the bottom of the simulated cell
                        ex = (flow.getP(eq_head, head, recharge, eqFlow) * flow.getBottom() +
                              flow.getQ(eq_head, head, recharge, eqFlow)) * get<t_dim, StepSize>();
                    }
                } else {  // GENERAL_HEAD_BOUNDARY (Question: what about FLOODPLAIN_DRAIN, EVAPOTRANSPIRATION, FAST_SURFACE_RUNOFF)
                    //LOG(debug) << "getting GHB flow";
                    ex = (flow.getP(eq_head, head, recharge, eqFlow) * head + // = - ghb_conductance * gw_head
                          flow.getQ(eq_head, head, recharge, eqFlow)) // = ghb_conductance * ghb_elevation (often = 0)
                                  * get<t_dim, StepSize>();
                    //LOG(debug) << "GHB flow:" << ex.value();
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
                t_vol_t out = 0 * si::cubic_meter / day;

                if (neighbours.find(DOWN) != neighbours.end()) {
                    large_num nodeIDDown = neighbours.find(DOWN)->first;
                    auto elev = getAt<t_meter, Elevation>(nodeIDDown);
                    auto head_n = getAt<t_meter, Head>(nodeIDDown);
                    //Check if a dewatered condition is present
                    if (head_n < elev and get<t_meter, Head>() > elev) {
                        t_s_meter_t conductance_below = mechanics.calculateVerticalConductance(createDataTuple(nodeIDDown));
                        out += conductance_below * (head_n - elev) * get<t_dim, StepSize>();
                    }
                }

                if (neighbours.find(TOP) != neighbours.end()) {
                    large_num nodeIDTop = neighbours.find(TOP)->first;
                    t_meter elev = getAt<t_meter, Elevation>(nodeIDTop);
                    t_meter head_n = getAt<t_meter, Head>(nodeIDTop);
                    //Check if a dewatered condition is present
                    if (get<t_meter, Head>() < get<t_meter, Elevation>() and head_n > elev) {
                        t_s_meter_t conductance_above =
                                mechanics.calculateVerticalConductance(createDataTuple(nodeIDTop));
                        out += conductance_above * (get<t_meter, Elevation>() - get<t_meter, Head>()) *
                               get<t_dim, StepSize>();
                    }
                }
                NANChecker(out.value(), "Dewatered flow");
                return out;
            }

            /**
             * @brief Get all current IN flow
             * @return Flow volume
             */
            t_c_meter getCurrentIN() noexcept {
                return getFlow([](double a) -> bool { return a > 0; }) * day;
            }

            /**
             * @brief Get all current OUT flow
             * @return Flow volume
             */
            t_c_meter getCurrentOUT() noexcept {
                return -getFlow([](double a) -> bool { return a < 0; }) * day;
            }

            /**
             * @brief Tell cell to save its flow budget
             */
            void saveMassBalance() noexcept {
                fields.addTo<t_c_meter, OUT>(getCurrentOUT());
                fields.addTo<t_c_meter, IN>(getCurrentIN());
            }


            /**
             * @brief Calculates zone change budget of a node at localZetaID (for variable density flow budget)
             * @return Zone change budget
             * @note in SWI2 code: SSWI2_ZCHG
             */
            t_c_meter calculateZoneChange(int localZetaID) {
                t_c_meter out = 0.0 * si::cubic_meter;

                t_meter zeta = getZeta(localZetaID);
                t_meter zetaBelow = getZeta(localZetaID + 1);
                // if zeta below next zeta surface
                if (zeta < zetaBelow){ // Question: potential problem here? (by clangevin)
                    // set zeta to height of next zeta
                    zeta = zetaBelow;
                }

                // if head is below bottom of node
                if (getHead() < getBottom()) {
                    // set zetas to bottom of node (zone thickness = zero)
                    zeta = getBottom();
                    zetaBelow = getBottom();
                }

                t_meter zetaOld = getZetaTZero(localZetaID);
                t_meter zetaBelowOld = getZetaTZero(localZetaID + 1);
                // if old zeta is below next old zeta surface
                if (zetaOld < zetaBelowOld) {
                    // set zeta old to height of next zeta old
                    zetaOld = zetaBelowOld;
                }

                // if old head is below bottom of node
                if (getHead_TZero() < getBottom()) {
                    // set zetas to bottom of node (zone thickness = zero)
                    zetaOld = getBottom();
                    zetaBelowOld = getBottom();
                }

                // calculate zone change ((zeta-zetaBelow) * Area is the current zone volume, with old zetas old vol)
                // not dividing by (day * step size) since we want to get the volumetric budget over full time step
                out += (getEffectivePorosity() * get<t_s_meter, Area>()) *
                        ((zeta - zetaBelow) - (zetaOld - zetaBelowOld));
                return out;
            }


            /**
             * @brief Instantaneous mixing of water, as described in SWI2 documentation under "Vertical Leakage Between
             * Aquifers": (2) when freshwater leaks up into an aquifer containing only saline water, that freshwater is
             *                added as saline water
             *            (4) when saline water leaks down into an aquifer containing only freshwater, that saline water
             *                is added as freshwater
             * @param localZetaID
             * @note in SWI2 code: SSWI2_IMIX, comments referring to lines (at each "if"), refer to lines in gwf2swi27.f
             */
            t_c_meter calculateInstantaneousMixing(int localZetaID) {
                t_c_meter out = 0.0 * si::cubic_meter;
                // skip nodes that do not have a down neighbour
                if (neighbours.find(DOWN) == neighbours.end()) { return out; } // line 4142

                // skip if dimensionless density at bottom of this node is below or equal to
                // dimension less density at the top of down node
                large_num nodeIDDown = neighbours[DOWN];
                if (getNusBot() <= at(nodeIDDown)->getNusTop()){ return out; } // line 4150

                // skip if head in this or neighbour node is below node bottom
                if (getHead() < getBottom() or at(nodeIDDown)->getHead() < at(nodeIDDown)->getBottom()){ return out; } // lines 4140 and 4145

                // skip if localZetaID is not involved
                auto nusInZones = getNusInZones();
                if (nusInZones[localZetaID] != getNusBot()) { return out; } // line 4156
                if (nusInZones[localZetaID] != at(nodeIDDown)->getNusTop()) { return out; } // line 4162

                // calculate the flux down
                t_c_meter fluxDown = getFluxDown().value() * si::cubic_meter;

                // if flux down is positive
                if (fluxDown > (0 * si::cubic_meter)) { // line 4168
                    // skip if dim-less density at bottom of down neighbour is greater or equal to
                    // dim-less density at the bottom of this node
                    if (at(nodeIDDown)->getNusBot() >= getNusBot()) { return out; } // line 4169

                    // skip if dim-less density at top of this node is smaller or equal to
                    // dim-less density at the top of neighbour node
                    if (getNusTop() <= at(nodeIDDown)->getNusTop()){ return out; } // line 4172
                }

                return fluxDown;
            }

            t_c_meter getInstantaneousMixing(bool in) noexcept {
                t_c_meter iMix_in;
                t_c_meter iMix_out;
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; ++localZetaID) {
                    t_c_meter iMix = calculateInstantaneousMixing(localZetaID);
                    if (iMix.value() > 0) { iMix_in += iMix; } else { iMix_out += iMix; }
                }
                if (in) { return iMix_in; } else { return iMix_out; }
            }

            /**
             * @brief save volumetric density zone change between last and new time step
             */
            void saveZoneChange() noexcept {
                set<t_c_meter, ZCHG_IN>(getZoneChange(true));
                set<t_c_meter, ZCHG_OUT>(getZoneChange(false));
            }

            t_c_meter getZoneChange(bool in) noexcept {
                t_c_meter zoneChange_in = 0 * si::cubic_meters;
                t_c_meter zoneChange_out = 0 * si::cubic_meters;
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; ++localZetaID) {
                    t_c_meter zoneChange = calculateZoneChange(localZetaID);
                    if (zoneChange.value() > 0) { zoneChange_in += zoneChange; } else { zoneChange_out += zoneChange; }
                }
                if (in) { return zoneChange_in; } else { return zoneChange_out; }
            }

            t_c_meter getTipToeTrackingZoneChange(bool in) {
                t_c_meter result = 0 * si::cubic_meters;
                t_c_meter tttOut = getZoneChange(false) - getZCHG_OUT();
                t_c_meter tttIn = getZoneChange(true) - getZCHG_IN();
                if (in) {
                    if (tttOut.value() > 0) { result += tttOut; } // higher by 3 orders of magnitude
                    if (tttIn.value() > 0) { result += tttIn; }
                } else {
                    if (tttOut.value() < 0) { result += tttOut; } // lower
                    if (tttIn.value() < 0) { result += tttIn; }
                }
                return result;
            }

            /**
             * @brief save variable density flow mass balance in node property (called before and after adjustment)
             * @param
             */
            t_c_meter getCurrentIN_VDF() {
                t_c_meter vdfIn = 0 * si::cubic_meter;
                bool in = true;
                vdfIn += getZCHG_IN(); // add current zone change before tip toe tracking
                vdfIn += getInstantaneousMixing(in); // add instantaneous mixing
                vdfIn += getTipToeTrackingZoneChange(in); // add zone change from tiptoetracking
                return vdfIn;
            }

            t_c_meter getCurrentOUT_VDF() {
                t_c_meter vdfOut = 0 * si::cubic_meter;
                bool out = false;
                vdfOut += getZCHG_OUT(); // add current zone change before tip toe tracking
                vdfOut += getInstantaneousMixing(out); // add instantaneous mixing
                vdfOut += getTipToeTrackingZoneChange(out); // add zone change from tiptoetracking
                return vdfOut;
            }

            t_c_meter getZCHG_OUT() { return get<t_c_meter, ZCHG_OUT>(); }

            t_c_meter getZCHG_IN() { return get<t_c_meter, ZCHG_IN>(); }

            void saveGNCMassBalance(){
                t_c_meter gncOut = getGNC_OUT();
                t_c_meter gncIn = getGNC_IN();

                t_c_meter gncFromUnrefinedNodes = getGNCFromNodes().value() * si::cubic_meter;
                t_c_meter gncToRefinedNode = getGNCToRefinedNode().value() * si::cubic_meter;

                if (gncFromUnrefinedNodes.value() < 0) {
                    gncOut += gncFromUnrefinedNodes;
                } else {
                    gncIn += gncFromUnrefinedNodes;
                }

                if (gncToRefinedNode.value() < 0) {
                    gncOut += gncToRefinedNode;
                } else {
                    gncIn += gncToRefinedNode;
                }

                set<t_c_meter, GNC_OUT>(gncOut);
                set<t_c_meter, GNC_IN>(gncIn);
            }

            t_c_meter getGNC_OUT() { return get<t_c_meter, GNC_OUT>(); }

            t_c_meter getGNC_IN() { return get<t_c_meter, GNC_IN>(); }

            /**
             * @brief Add a neighbour
             * @param ID The internal ID and position in vector
             * @param neighbour The position relative to the cell
             */
            void setNeighbour(large_num ID, NeighbourPosition neighbourPosition) {
                neighbours[neighbourPosition] = ID;
                if (neighbourPosition == NeighbourPosition::TOP or neighbourPosition == NeighbourPosition::DOWN) {
                    vertical_neighbours[neighbourPosition] = ID;
                } else {
                    horizontal_neighbours[neighbourPosition] = ID;
                }
                if (neighbourPosition != NeighbourPosition::TOP and neighbourPosition != NeighbourPosition::DOWN and
                    neighbourPosition != NeighbourPosition::FRONT and neighbourPosition != NeighbourPosition::BACK and
                    neighbourPosition != NeighbourPosition::LEFT and neighbourPosition != NeighbourPosition::RIGHT) {
                    finer_neighbours.push_back(neighbourPosition);
                }
            }

            /**
             * @brief Add multiple neighbours in one direction
             * @param nodeIDs The neighbouring nodeIDs
             * @param neighbour The neighbour position relative to the cell
             */
            void setNeighbours(std::unordered_map<large_num, large_num> nodeIDs_neig,
                               NeighbourPosition neighbourPosition) {
                if (nodeIDs_neig.empty()) { // do nothing
                    return;
                } else if (nodeIDs_neig.size() == 1) { // if only one node at this spatID
                    setNeighbour(nodeIDs_neig.at(0), neighbourPosition);
                } else { // if more than one node at this spatID
                    std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
                            mapToNeig;

                    if (nodeIDs_neig.size() == 4){
                        mapToNeig = getMapNeighboursRefinedToFour(neighbourPosition);
                    } else if (nodeIDs_neig.size() == 9){
                        mapToNeig = getMapNeighboursRefinedToNine(neighbourPosition);
                    } else if (nodeIDs_neig.size() == 16){
                        mapToNeig = getMapNeighboursRefinedToSixteen(neighbourPosition);
                    }else {
                        return;
                    }

                    // loop over neigbour nodes at one spatID
                    for (auto nodeID_neig: nodeIDs_neig) {
                        //LOG(debug) << "setting neighbours for nodeID " << nodeID.second;
                        auto refID = nodes->at(nodeID_neig.second)->get<large_num, RefID>();
                        try{
                            auto refNeigPos = mapToNeig.at(neighbourPosition).at(refID);
                            setNeighbour(nodeID_neig.second, refNeigPos);
                        } catch (const std::out_of_range &ex) {
                            continue;
                        }
                    }
                }
            }

            std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
            getMapNeighboursRefinedToFour(Model::NeighbourPosition neighbourPosition){
                /*
                 * RefIDs:
                 *                  (this)
                 *                  BACK
                 *                  ||  ||
                 *                  \/  \/
                 *              =>  1   2  <=
                 * (this) RIGHT     (neig)    LEFT (this)
                 *              =>  3   4  <=
                 *                  /\  /\
                 *                  ||  ||
                 *                  FRONT
                 *                  (this)
                 */

                // define a map that maps neighbour position to refID to refined neighbour position
                std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
                        mapNeigPosToRefIdToRefNeigPos = { { Model::FRONT, { {3, Model::FRONTLEFT},  {4, Model::FRONTRIGHT} } },
                                                          { Model::BACK,  { {1, Model::BACKLEFT},   {2, Model::BACKRIGHT}  } },
                                                          { Model::RIGHT, { {1, Model::RIGHTFRONT}, {3, Model::RIGHTBACK}  } },
                                                          { Model::LEFT,  { {2, Model::LEFTFRONT},  {4, Model::LEFTBACK}   } } };
                return mapNeigPosToRefIdToRefNeigPos;
            }

            std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
            getMapNeighboursRefinedToNine(Model::NeighbourPosition neighbourPosition){
                /*
                 * RefIDs:
                 *                  (this)
                 *                  BACK
                 *                  ||  || ||
                 *                  \/  \/ \/
                 *              =>  1   2  3 <=
                 * (this) RIGHT =>  4  [5] 6 <= LEFT (this)
                 *                    (neig)
                 *              =>  7   8  9 <=
                 *                  /\  /\ /\
                 *                  ||  || ||
                 *                  FRONT
                 *                  (this)
                 */

                // define a map that maps neighbour position to refID to refined neighbour position
                std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
                        mapNeigPosToRefIdToRefNeigPos = {
                        { Model::FRONT, { {7, Model::FRONTLEFT},  {8, Model::FRONTFRONT}, {9, Model::FRONTRIGHT} } },
                        { Model::BACK,  { {1, Model::BACKLEFT},   {2, Model::BACKBACK},  {3, Model::BACKRIGHT}  } },
                        { Model::RIGHT, { {1, Model::RIGHTFRONT}, {4, Model::RIGHTRIGHT}, {7, Model::RIGHTBACK}  } },
                        { Model::LEFT,  { {3, Model::LEFTFRONT},  {6, Model::LEFTLEFT},  {9, Model::LEFTBACK}   } } };
                return mapNeigPosToRefIdToRefNeigPos;
            }

            std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
            getMapNeighboursRefinedToSixteen(Model::NeighbourPosition neighbourPosition){
                /*
                 * RefIDs:
                 *                  (this)
                 *                  BACK
                 *                  ||  ||  ||  ||
                 *                  \/  \/  \/  \/
                 *              =>   1   2   3   4   <=
                 * (this) RIGHT =>   5   6   7   8   <= LEFT (this)
                 *                      (neig)
                 *              =>   9  10  11  12   <=
                 *              =>  13  14  15  16   <=
                 *                  /\  /\  /\  /\
                 *                  ||  ||  ||  ||
                 *                  FRONT
                 *                  (this)
                 */

                // define a map that maps neighbour position to refID to refined neighbour position
                std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, Model::NeighbourPosition>>
                        mapNeigPosToRefIdToRefNeigPos = {
                        { Model::FRONT, { {13, Model::FRONTLEFT},       {14, Model::FRONTFRONTLEFT},
                                                {15, Model::FRONTFRONTRIGHT}, {16, Model::FRONTRIGHT} } },
                        { Model::BACK,  { {1, Model::BACKLEFT},         {2, Model::BACKBACKLEFT},
                                                {3, Model::BACKBACKRIGHT},    {4, Model::BACKRIGHT}  } },
                        { Model::RIGHT, { {1, Model::RIGHTFRONT},       {5, Model::RIGHTRIGHTFRONT},
                                                {9, Model::RIGHTRIGHTBACK},   {13, Model::RIGHTBACK}  } },
                        { Model::LEFT,  { {4, Model::LEFTFRONT},        {8, Model::LEFTLEFTFRONT},
                                                {12, Model::LEFTLEFTBACK},    {16, Model::LEFTBACK}   } } };
                return mapNeigPosToRefIdToRefNeigPos;
            }

            int getNumOfNeighbours() { return (int) neighbours.size(); }

            int getNumOfHorizontalNeighbours() { return (int) horizontal_neighbours.size(); }

            class NodeNotFoundException : public std::exception {
                virtual const char *what() const throw() { return "Node does not exist"; }
            };

            std::unordered_map<NeighbourPosition, large_num> getListOfNeighbours(){
                return neighbours;
            }

            std::unordered_map<NeighbourPosition, large_num> getListOfHorizontalNeighbours(){
                return horizontal_neighbours;
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
                    //LOG(debug) << "Adding flow " << std::to_string(type) << " that already existed at nodeID: " << getID();
                    //currently it is assumed that only one external flow of one type is what we want
                    // FIXME if not we have to replace the enum with something different
                    removeExternalFlow(type);
                }

                if (type == RECHARGE or type == FAST_SURFACE_RUNOFF or type == NET_ABSTRACTION) {
                    externalFlows.insert(
                            std::make_pair(type, ExternalFlow(numOfExternalFlows, cond * (si::cubic_meter / day), type)));
                } else if (type == EVAPOTRANSPIRATION) {
                    externalFlows.insert(std::make_pair(type,
                                                        ExternalFlow(numOfExternalFlows, flowHead, bottom,
                                                                     cond * (si::cubic_meter / day))));
                    /** TODO Implementation of FLOODPLAIN_DRAIN
                    *} else if (type == FLOODPLAIN_DRAIN) {
                    *    externalFlows.insert(std::make_pair(type,
                    *                                        ExternalFlow(numOfExternalFlows, type,
                    *                                                        get<t_meter, Elevation>(),
                    *                                                        get<t_vel, K>() * get<t_meter,
                    *                                                        VerticalSize>(),
                    *                                                        bottom));
                    */
                } else { // RIVER, RIVER_MM, DRAIN, WETLAND, GLOBAL_WETLAND, LAKE, GLOBAL_LAKE, GENERAL_HEAD_BOUNDARY
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
             * @param type The flow id
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
             * @return int
             */
            int getNumOfExternalFlows(){ return numOfExternalFlows; }


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
             * @brief Ghost Node Correction following MODFLOW-USG documentation by Panday et al. (2013) and Fortran file
             * "disu2gncn.f" of the MODFLOW-USG package. Using symmetric implementation, adding to RHS (lines 208-211)
             * @return
             */
            t_vol_t getGNCFromNodes(){
                //LOG(debug) << "getGNCFromNodes";
                return -calculateGhostNodeCorrection(finer_neighbours);
            }

            /**
             * @brief Ghost Node Correction following MODFLOW-USG documentation by Panday et al. (2013) and Fortran file
             * "disu2gncn.f" of the MODFLOW-USG package. Using symmetric implementation, adding to RHS (lines 208-211)
             * @return
             */
            t_vol_t getGNCToRefinedNode(){
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                if (getRefinedInto() == 1) { return out;} // if this node is not refined -> return

                // get neighbours at left, right, front or back
                std::vector<NeighbourPosition> neigPos_LRFB = getNeigPos_LRFB();

                for (const auto &neigPos: neigPos_LRFB) {
                    if (horizontal_neighbours.find(neigPos) != horizontal_neighbours.end()) {
                        // if neighbour is unrefined
                        large_num neigNodeID = horizontal_neighbours[neigPos];
                        if (at(neigNodeID)->getRefinedInto() == 1) { // if neig is not refined
                            // get the refined neighbour position this node has relative to that unrefined neighbour
                            NeighbourPosition refNeigPos = getRefNeigPosToUnrefNeig(get<large_num, RefID>(), neigPos,
                                                                                    getRefinedInto());
                            //LOG(debug) << "getGNCToRefinedNode, nodeID: " << getID() << ", neig nodeID: " << at(neig)->getID();
                            out += at(neigNodeID)->calculateGhostNodeCorrection({refNeigPos});
                        }
                    }
                }
                return out;
            }

            /**
             * @brief Get the refined neighbour position this node has relative to that unrefined neighbour
             * @param refID refinement Identifier
             * @param neigPos neighbour position
             * @param refinedInto number of nodes the original cell is split into
             * @return
             */
            NeighbourPosition getRefNeigPosToUnrefNeig(large_num refID, NeighbourPosition neigPos,
                                                       large_num refinedInto){
                std::unordered_map<large_num, std::unordered_map<NeighbourPosition, NeighbourPosition>>
                        mapRefIdToNeigToRefNeig;
                if (refinedInto == 4) {
                    mapRefIdToNeigToRefNeig = {
                            {1, { {NeighbourPosition::FRONT, NeighbourPosition::BACKLEFT},
                                        {NeighbourPosition::LEFT,  NeighbourPosition::RIGHTFRONT} } },
                            {2, { {NeighbourPosition::FRONT, NeighbourPosition::BACKRIGHT},
                                        {NeighbourPosition::RIGHT, NeighbourPosition::LEFTFRONT} } },
                            {3, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTLEFT},
                                        {NeighbourPosition::LEFT,  NeighbourPosition::RIGHTBACK} } },
                            {4, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTRIGHT},
                                        {NeighbourPosition::RIGHT, NeighbourPosition::LEFTBACK} } } };
                } else if (refinedInto == 9) {
                    mapRefIdToNeigToRefNeig = {
                            {1, { {NeighbourPosition::FRONT, NeighbourPosition::BACKLEFT},
                                        {NeighbourPosition::LEFT,  NeighbourPosition::RIGHTFRONT} } },
                            {2, { {NeighbourPosition::FRONT, NeighbourPosition::BACKBACK} } },
                            {3, { {NeighbourPosition::FRONT,  NeighbourPosition::BACKRIGHT},
                                        {NeighbourPosition::RIGHT,  NeighbourPosition::LEFTFRONT} } },
                            {4, { {NeighbourPosition::LEFT,  NeighbourPosition::RIGHTRIGHT} } },
                            {5, { {} } }, // has only refined neighbours
                            {6, { {NeighbourPosition::RIGHT, NeighbourPosition::LEFTLEFT} } },
                            {7, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTLEFT},
                                        {NeighbourPosition::LEFT,  NeighbourPosition::RIGHTBACK} } },
                            {8, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTFRONT} } },
                            {9, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTRIGHT},
                                        {NeighbourPosition::RIGHT, NeighbourPosition::LEFTBACK} } } };
                } else if (refinedInto == 16) {
                    mapRefIdToNeigToRefNeig = {
                            {1, { {NeighbourPosition::FRONT,  NeighbourPosition::BACKLEFT},
                                        {NeighbourPosition::LEFT,   NeighbourPosition::RIGHTFRONT} } },
                            {2, { {NeighbourPosition::FRONT,  NeighbourPosition::BACKBACKLEFT} } },
                            {3, { {NeighbourPosition::FRONT,  NeighbourPosition::BACKBACKRIGHT} } },
                            {4, { {NeighbourPosition::FRONT,  NeighbourPosition::BACKRIGHT},
                                        {NeighbourPosition::RIGHT,  NeighbourPosition::LEFTFRONT} } },
                            {5, { {NeighbourPosition::LEFT,   NeighbourPosition::RIGHTRIGHTFRONT} } },
                            {6, { {} } }, // has only refined neighbours
                            {7, { {} } }, // has only refined neighbours
                            {8, { {NeighbourPosition::RIGHT,  NeighbourPosition::LEFTLEFTFRONT} } },
                            {9, { {NeighbourPosition::LEFT,   NeighbourPosition::RIGHTRIGHTBACK} } },
                            {10, { {} } }, // has only refined neighbours
                            {11, { {} } }, // has only refined neighbours
                            {12, { {NeighbourPosition::RIGHT, NeighbourPosition::LEFTLEFTBACK} } },
                            {13, { {NeighbourPosition::BACK, NeighbourPosition::FRONTLEFT},
                                        {NeighbourPosition::LEFT, NeighbourPosition::RIGHTBACK} } },
                            {14, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTFRONTLEFT} } },
                            {15, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTFRONTRIGHT} } },
                            {16, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTRIGHT},
                                        {NeighbourPosition::RIGHT, NeighbourPosition::LEFTBACK} } } };
                }
                return mapRefIdToNeigToRefNeig.at(refID).at(neigPos);
            }

            t_vol_t calculateGhostNodeCorrection(const std::vector<NeighbourPosition>& possibleRefNeigPos){
                t_vol_t gnc;
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                t_meter head = getHead();
                t_dim multiplierContributor{};
                t_dim multiplierNodeInner{};
                t_dim multiplierNodeOuter{};
                t_s_meter_t transmissivitySelf = get<t_meter, VerticalSize>() * getK();
                t_s_meter_t contributorConductance = 0.0 * (si::square_meter / day);
                t_s_meter_t nodeConductance = 0.0 * (si::square_meter / day);
                std::vector<NeighbourPosition> potentialContributors;

                int neigRefInto{};
                int contrRefInto{};
                for (const auto &refNeigPos: possibleRefNeigPos) {
                    if (horizontal_neighbours.find(refNeigPos) != horizontal_neighbours.end()) {
                        large_num refNeigNodeID = horizontal_neighbours[refNeigPos];
                        potentialContributors = getPotentialContributors(refNeigPos);
                        for (const auto &contrPos: potentialContributors) {
                            if (horizontal_neighbours.find(contrPos) == horizontal_neighbours.end()) { // at model boundary: no contributor
                                continue; // todo implement impact of GHB (not required if all nodes at coast/GHB are refined)
                            }
                            large_num contrNodeID = horizontal_neighbours[contrPos]; // contributor is named "j" in USG doc

                            contrRefInto = (int) at(contrNodeID)->getRefinedInto();
                            multiplierContributor = ( 1.0 / (2.0 * sqrt(contrRefInto) )) * si::si_dimensionless;

                            neigRefInto = (int) at(contrNodeID)->getRefinedInto();
                            if (neigRefInto == 4) {
                                multiplierNodeInner = (1.0 / 4.0) * si::si_dimensionless;
                                multiplierNodeOuter = (1.0 / 4.0) * si::si_dimensionless;
                            } else if (neigRefInto == 9) {
                                multiplierNodeInner = (2.0 / 6.0) * si::si_dimensionless;
                                multiplierNodeOuter = (1.0 / 6.0) * si::si_dimensionless;
                            } else if (neigRefInto == 16) {
                                if (refNeigPos > 18) { // todo debug
                                    multiplierNodeInner = (1.0 / 8.0) * si::si_dimensionless;
                                    multiplierNodeOuter = (3.0 / 8.0) * si::si_dimensionless;
                                } else {
                                    multiplierNodeInner = (3.0 / 8.0) * si::si_dimensionless;
                                    multiplierNodeOuter = (1.0 / 8.0) * si::si_dimensionless;
                                }
                            }

                            t_s_meter_t transmissivityNeig =
                                    getAt<t_meter, VerticalSize>(contrNodeID) * at(contrNodeID)->getK();

                            if (transmissivityNeig != 0 * si::square_meter / day and
                                transmissivitySelf != 0 * si::square_meter / day) {
                                t_meter nodeWidth = std::min(getNodeWidth(contrPos),
                                                             at(contrNodeID)->getNodeWidth(contrPos));
                                // conductance from contributor node to ghost node
                                // if refined into four, multiplierContributor is 0.5, multiplierNodeOuter is 0.25
                                contributorConductance = nodeWidth *
                                                         ((transmissivitySelf * transmissivityNeig)
                                                          / (transmissivitySelf *
                                                             at(contrNodeID)->getNodeLength(contrPos) *
                                                             multiplierContributor +
                                                             transmissivityNeig * getNodeLength(contrPos) *
                                                             multiplierNodeOuter));
                            }
                            // conductance from this node's center to ghost node inside this node
                            // if refined into four, multiplierNodeInner is 0.25
                            nodeConductance = getNodeWidth(contrPos) *
                                              (transmissivitySelf / (getNodeLength(contrPos) * multiplierNodeInner));

                            // the alpha coefficient is used to weigh influence on ghost node height difference
                            t_dim alpha = contributorConductance / (contributorConductance + nodeConductance);
                            //LOG(debug) << "alpha = " << alpha << ", contributorConductance = " << contributorConductance.value();

                            // conductance to the refined neighbour node
                            t_s_meter_t condRefNeig =
                                    mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(refNeigPos, refNeigNodeID));

                            // calculate ghost node correction
                            gnc = condRefNeig * (alpha * (getHead() - at(contrNodeID)->getHead()));
                            out += gnc;
                            //LOG(debug) << "GNC for nodeID "     << get<large_num, ID>() <<
                            //              " to refined nodeID " << getAt<large_num, ID>(refinedNeig) <<
                            //              " with contributor "  << getAt<large_num, ID>(contributor) <<
                            //              " = "                 << gnc.value();
                        }
                    }
                }
                return out;
            }

            static std::vector<NeighbourPosition>
            getPotentialContributors(NeighbourPosition refinedNeigPos){
                if (refinedNeigPos == NeighbourPosition::FRONTLEFT or
                    refinedNeigPos == NeighbourPosition::BACKLEFT or
                    refinedNeigPos == NeighbourPosition::FRONTFRONTLEFT or
                    refinedNeigPos == NeighbourPosition::BACKBACKLEFT) {
                    return {NeighbourPosition::LEFT, NeighbourPosition::LEFTFRONT, NeighbourPosition::LEFTBACK,
                            NeighbourPosition::LEFTLEFTFRONT, NeighbourPosition::LEFTLEFTBACK};
                } else if (refinedNeigPos == NeighbourPosition::FRONTRIGHT or
                           refinedNeigPos == NeighbourPosition::BACKRIGHT or
                           refinedNeigPos == NeighbourPosition::FRONTFRONTRIGHT or
                           refinedNeigPos == NeighbourPosition::BACKBACKRIGHT) {
                    return {NeighbourPosition::RIGHT, NeighbourPosition::RIGHTFRONT, NeighbourPosition::RIGHTBACK,
                            NeighbourPosition::RIGHTRIGHTFRONT, NeighbourPosition::RIGHTRIGHTBACK};
                } else if (refinedNeigPos == NeighbourPosition::LEFTFRONT or
                           refinedNeigPos == NeighbourPosition::RIGHTFRONT or
                           refinedNeigPos == NeighbourPosition::LEFTLEFTFRONT or
                           refinedNeigPos == NeighbourPosition::RIGHTRIGHTFRONT) {
                    return {NeighbourPosition::FRONT, NeighbourPosition::FRONTLEFT, NeighbourPosition::FRONTRIGHT,
                            NeighbourPosition::FRONTFRONTRIGHT, NeighbourPosition::FRONTFRONTLEFT};
                } else if (refinedNeigPos == NeighbourPosition::LEFTBACK or
                           refinedNeigPos == NeighbourPosition::RIGHTBACK or
                           refinedNeigPos == NeighbourPosition::LEFTLEFTBACK or
                           refinedNeigPos == NeighbourPosition::RIGHTRIGHTBACK) {
                    return {NeighbourPosition::BACK, NeighbourPosition::BACKLEFT, NeighbourPosition::BACKRIGHT,
                            NeighbourPosition::BACKBACKLEFT, NeighbourPosition::BACKBACKRIGHT};
                } else {
                    std::vector<NeighbourPosition> emptyVector{};
                    return emptyVector;
                }
            }

            static std::vector<NeighbourPosition>
            getNeigPos_LRFB(){
                return {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                        NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
            }


            void initializeZetas() {
                std::vector<t_meter> zetas;
                t_meter zetaTop;
                // add top zeta (at node elevation or gw head or node bottom)
                if (getHead() > getBottom()) {
                    if (getHead() > getElevation()) {
                        // top zeta is max at node elevation
                        zetaTop = getElevation();
                    } else {
                        // top zeta is at gw head
                        zetaTop = getHead();
                    }
                // if head is below node bottom
                } else {
                    zetaTop = getBottom();
                }
                zetas.push_back(zetaTop);
                // add bottom zeta at node bottom
                zetas.push_back(getBottom());
                setZetas(zetas);

                // add two zeta changes
                std::vector<t_meter> zetasChange;
                zetasChange.push_back(0 * si::meter);
                zetasChange.push_back(0 * si::meter);
                set<std::vector<t_meter>,ZetasChange>(zetasChange);
            }

            /**
             * @brief Add a zeta surface to the cell (bounded by elevation at top and by cell bottom at bottom).
             * @param localZetaID
             * @param zeta the zeta surface height in meters
             */
            void addZeta(int localZetaID, t_meter zeta){
                NANChecker(zeta.value(), "zeta (in addZeta)");
                if (localZetaID == 0) { throw "localZetaID should not be 0 when adding zetas"; }
                auto zetas = getZetas();
                auto zetasChange = getZetasChange();

                auto zetaLimited = applyZetaLimits(localZetaID, zeta);
                zetas.insert(zetas.begin() + localZetaID, zetaLimited);
                zetasChange.insert(zetasChange.begin() + localZetaID, 0 * si::meter);

                // check if zetas vector is sorted
                for (int zetaID = 0; zetaID < zetas.size() - 1; ++zetaID) {
                    if (zetas[zetaID] < zetas[zetaID+1]) {
                        LOG(userinfo) << "At nodeID " << getID() << ": zeta at ID= " << zetaID << ":" <<
                                      zetas[zetaID].value() << ", zeta at ID+1: " <<  zetas[zetaID+1].value();
                        throw "Vector of zetas needs to be sorted!";
                    }
                }

                setZetas(zetas);
                set<std::vector<t_meter>,Zetas_TZero>(zetas);
                set<std::vector<t_meter>,ZetasChange>(zetasChange);
            }

            /**
             * @brief Update zetas after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            void setZeta(int localZetaID, t_meter zeta) {
                NANChecker(zeta.value(), "height (in setZeta)");
                auto zetas = getZetas();
                if (localZetaID < zetas.size()) {
                    t_meter zetaLimited = applyZetaLimits(localZetaID, zeta);
                    zetas[localZetaID] = zetaLimited;
                    setZetas(zetas);
                } else {
                    throw "localZetaID too large in setZeta";
                }

                //LOG(debug) << "nodeID: " << getID() << ", setZeta[" << localZetaID << "] = " << getZeta(localZetaID).value();
            }


            t_meter applyZetaLimits(int localZetaID, t_meter zeta) {
                // the top zeta (localZetaID = 0) is limited by the top and bottom elevation of the node
                if (localZetaID == 0 and zeta <= getElevation() and zeta >= getBottom()) {
                    return zeta;
                }

                // all other zetas are limited by the zeta height of the top and bottom zetas
                // (in SWI2: lines 660-680)
                auto zetas = getZetas();
                if (zeta > zetas.front() - get<t_meter, VDFLock>()) {
                    zeta = zetas.front();
                    // limit zeta height to bottom zeta (in SWI2: lines 660-680)
                } else if (zeta < zetas.back() + get<t_meter, VDFLock>()) {
                    zeta = zetas.back();
                }
                return zeta;
            }

            /**
             * @brief Update zetas after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            void setZetas(std::vector<t_meter> zetas) { set<std::vector<t_meter>, Zetas>(zetas); }

            /**
             * @brief Update zeta change after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            void setZetaAndZetaChange(int localZetaID, t_meter potentialZetaChange){
                NANChecker(potentialZetaChange.value(), "potentialZetaChange (in setZetaAndZetaChange)");
                auto zetasChange = getZetasChange();
                if (localZetaID < zetasChange.size()) {
                    // get old zeta
                    t_meter oldZeta = getZeta(localZetaID);

                    // set new zeta using potential zeta change
                    t_meter potentialNewZeta = oldZeta + potentialZetaChange;
                    setZeta(localZetaID, potentialNewZeta);

                    // calculate and set actual zeta change
                    t_meter actualZetaChange = getZeta(localZetaID) - oldZeta;
                    zetasChange[localZetaID] = actualZetaChange;
                    set<std::vector<t_meter>, ZetasChange>(zetasChange);
                } else {
                    throw "localZetaID too large in setZetaChange";
                }
            }

            void setZetasTZero(){ set < std::vector<t_meter>, Zetas_TZero > (getZetas()); }

            /**
             * @brief get all zeta surface heights
             * @return vector<meter>
             */
            std::vector<t_meter> getZetas() { return get<std::vector<t_meter>, Zetas>();}

            /**
             * @brief get one zeta surface height
             * @param localZetaID
             * @return meter
             */
            t_meter getZeta(int localZetaID) {

                if (localZetaID < getZetas().size()){
                    auto zetas = getZetas();
                    auto zeta = zetas[localZetaID];
                    return zeta;
                } else {
                    throw "Not set at nodeID " + std::to_string(getID()) +
                          ": Zetas[localZetaID = " + std::to_string(localZetaID) + "]";
                }
            }

            /**
             * @brief get all zeta surface changes
             * @return vector<meter>
             */
            std::vector<t_meter> getZetasChange() { return get<std::vector<t_meter>, ZetasChange>();}

            /**
             * @brief get zeta surface change
             * @return meter
             */
            t_meter getZetaChange(int localZetaID){
                auto zetasChange = get<std::vector<t_meter>, ZetasChange>();
                if (localZetaID < zetasChange.size()){
                    return zetasChange[localZetaID];
                } else {
                    throw "Not set at nodeID " + std::to_string(getID()) +
                          ": ZetasChange[localZetaID = " + std::to_string(localZetaID) + "]";
                }
            }

            std::vector<t_meter> getZetasTZero() { return get<std::vector<t_meter>, Zetas_TZero>();}

            t_meter getZetaTZero(int localZetaID){
                if (localZetaID < get<std::vector<t_meter>, Zetas_TZero>().size()){
                    return get<std::vector<t_meter>, Zetas_TZero>()[localZetaID];
                } else {
                    throw "Not set at nodeID " + std::to_string(getID()) +
                          ": Zetas_TZero[localZetaID = " + std::to_string(localZetaID) + "]";
                }
            }

            std::vector<t_dim> getNusInZones() { return get<std::vector<t_dim>, NusInZones>();}

            int getSourceZoneGHB(){ return get<int, SourceZoneGHB>(); }

            int getSourceZoneRecharge(){ return get<int, SourceZoneRecharge>(); }

            /**
             * @brief Set effective porosity in node and layers below
             * @param effectivePorosity effective porosity in node
             */
            void setEffectivePorosity(t_dim effectivePorosity) {
                set<t_dim, EffectivePorosity>(effectivePorosity);
                applyToAllLayers([&effectivePorosity](NodeInterface *nodeInterface) {
                    nodeInterface->getProperties().set<t_dim, EffectivePorosity>(effectivePorosity);
                });
            }

            /**
             * @brief Updates GW recharge
             * Curently assumes only one recharge as external flow!
             * @param amount The new flow amount
             * @param Should the recharge in the dynamic rivers be locked or updated by this change?
             */
            void updateUniqueFlow(double amount, FlowType flow = RECHARGE, bool lock = true) {
                if (lock and flow == RECHARGE) {
                    if (hasTypeOfExternalFlow(RIVER_MM)) {
                        //get current recharge and lock it before setting new recharge
                        //in arid regions recharge might be 0
                        t_vol_t recharge{0 * si::cubic_meter /day};
                        if(hasTypeOfExternalFlow(RECHARGE)){recharge = getExternalFlowByName(RECHARGE).getRecharge();}
                        //also lock conductance value
                        getExternalFlowByName(RIVER_MM).getERC(recharge,get<t_meter, EQHead>(),get<t_meter, Head>(),getEqFlow());
                        getExternalFlowByName(RIVER_MM).setLockRecharge(recharge);
                        getExternalFlowByName(RIVER_MM).setLock(); //locks conductance to steady state conductance and inhibits updates later
                        //!comment! if this code is deactivated locked conductance and locked recharge is lost if flow is removed in addExternalFlowFlowHead but not important bc. never used (in calcERC) if conductance should be changed by calcERC
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
            }

            /**
             * @brief Update wetlands, lakes conduct if conduct is set to 0 due to reductionFactor -> reset to initial conduct from steady state to allow GW flow into SWB
             * @param reductionFactor 0 <= reductionFactor <= 1
             * @param type type of flow
             */
            void updateExternalFlowConduct(double reductionFactor, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = getExternalFlowByName(type).getFlowHead();
                    t_meter bottom = getExternalFlowByName(type).getBottom();
                    double RiverDepth = getExternalFlowByName(type).getRiverDepthSteadyState();
                    double initConduct = getExternalFlowByName(type).getInitConductance().value();
                    double conduct = initConduct;
                    if (reductionFactor != 0.)
                        conduct *= reductionFactor;
                    // else conduct does not change and is initial value
                    removeExternalFlow(type);
                    addExternalFlow(type, flowHead, conduct, bottom);
                    getExternalFlowByName(type).setInitConductance(initConduct);
                    getExternalFlowByName(type).setRiverDepthSteadyState(RiverDepth);
                }
            }

            /**
             *@brief saves conductance of steady state solution
             */
            void saveNodeSteadyStateConduct() {
                double conduct;
                if (hasTypeOfExternalFlow(Model::FlowType::LAKE)) {
                    conduct = getExternalFlowByName(Model::FlowType::LAKE).getConductance().value();
                    getExternalFlowByName(Model::FlowType::LAKE).setInitConductance(conduct);
                }
                if (hasTypeOfExternalFlow(Model::FlowType::WETLAND)) {
                    conduct = getExternalFlowByName(Model::FlowType::WETLAND).getConductance().value();
                    getExternalFlowByName(Model::FlowType::WETLAND).setInitConductance(conduct);
                }
                if (hasTypeOfExternalFlow(Model::FlowType::GLOBAL_LAKE)) {
                    conduct = getExternalFlowByName(Model::FlowType::GLOBAL_LAKE).getConductance().value();
                    getExternalFlowByName(Model::FlowType::GLOBAL_LAKE).setInitConductance(conduct);
                }
                if (hasTypeOfExternalFlow(Model::FlowType::GLOBAL_WETLAND)) {
                    conduct = getExternalFlowByName(Model::FlowType::GLOBAL_WETLAND).getConductance().value();
                    getExternalFlowByName(Model::FlowType::GLOBAL_WETLAND).setInitConductance(conduct);
                }
            }

            /**
            *@brief saves river depth (flow head - bottom) of steady state solution
            */
            void saveNodeSteadyStateRiverDepth() {
                double RiverDepth;
                if (hasTypeOfExternalFlow(Model::FlowType::RIVER_MM)){
                    RiverDepth = getExternalFlowByName(Model::FlowType::RIVER_MM).getFlowHead().value() - getExternalFlowByName(Model::FlowType::RIVER_MM).getBottom().value();
                    if (RiverDepth < 1.)    // else if 0 head could never change
                        RiverDepth = 1.;
                    getExternalFlowByName(Model::FlowType::RIVER_MM).setRiverDepthSteadyState(RiverDepth);
                }
            };

            /**
            * @brief Multiplies flow head for Sensitivity An. wetlands, lakes, rivers
            * @param amount
            * @param type
            */
            void updateExternalFlowFlowHead(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = getExternalFlowByName(type).getFlowHead() * amount; //TODO: if amount puts flowhead down head might be under bottom
                    double conduct = getExternalFlowByName(type).getConductance().value();
                    double RiverDepth = getExternalFlowByName(type).getRiverDepthSteadyState();
                    t_meter bottom = getExternalFlowByName(type).getBottom();
                    double initConduct = getExternalFlowByName(type).getInitConductance().value();
                    removeExternalFlow(type);
                    addExternalFlow(type, flowHead, conduct, bottom);
                    getExternalFlowByName(type).setInitConductance(initConduct);
                    getExternalFlowByName(type).setRiverDepthSteadyState(RiverDepth);
                }
            }

            /**
            * @brief Sets flowHead for wetlands, lakes, rivers; only used for global lakes at the moment
            * @param amount
            * @param type
            */
            void setExternalFlowFlowHead(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = amount * si::meter;
                    double conduct = getExternalFlowByName(type).getConductance().value();
                    double RiverDepth = getExternalFlowByName(type).getRiverDepthSteadyState();
                    t_meter bottom = getExternalFlowByName(type).getBottom();
                    double initConduct = getExternalFlowByName(type).getInitConductance().value();
                    removeExternalFlow(type);
                    if (flowHead.value() < bottom.value())
                        flowHead = bottom;
                    addExternalFlow(type, flowHead, conduct, bottom);
                    getExternalFlowByName(type).setInitConductance(initConduct);
                    getExternalFlowByName(type).setRiverDepthSteadyState(RiverDepth);
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
                    double initConduct = externalFlow.getInitConductance().value();
                    double RiverDepth = getExternalFlowByName(type).getRiverDepthSteadyState();
                    //The river is dry
                    if(std::isnan(amount)){ flowHead = bottom; }
                    double conduct{externalFlow.getConductance().value()};
                    bool lock{externalFlow.getLock()};
                    t_vol_t recharge{externalFlow.getLockRecharge()};
                    t_s_meter_t l_cond{externalFlow.getLockConduct()};
                    removeExternalFlow(type);
                    if (flowHead.value() < bottom.value())
                        flowHead = bottom;
                    NANChecker(flowHead.value(), "Stage value");
                    NANChecker(l_cond.value(), "Conduct value");
                    NANChecker(bottom.value(), "Bottom value");

                    addExternalFlow(type, flowHead, conduct, bottom);
                    getExternalFlowByName(type).setInitConductance(initConduct);
                    getExternalFlowByName(type).setRiverDepthSteadyState(RiverDepth);
                    if (lock) {
                        getExternalFlowByName(type).setLock();
                        getExternalFlowByName(type).setLockRecharge(recharge);
                        getExternalFlowByName(type).setLockConduct(l_cond);
                    }
                }
            }

            /**
             * @brief Multiplies flow bottom for Sensitivity An. wetlands, lakes, rivers
             * @param amount
             * @param type
             */
            void updateExternalFlowBottom(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = getExternalFlowByName(type).getFlowHead();
                    double conduct = getExternalFlowByName(type).getConductance().value();
                    double RiverDepth = getExternalFlowByName(type).getRiverDepthSteadyState();
                    t_meter bottom = getExternalFlowByName(type).getBottom() * amount; // TODO: if amount puts bottom upwards head might be under bottom
                    double initConduct = getExternalFlowByName(type).getInitConductance().value();
                    removeExternalFlow(type);
                    addExternalFlow(type, flowHead, conduct, bottom);
                    getExternalFlowByName(type).setInitConductance(initConduct);
                    getExternalFlowByName(type).setRiverDepthSteadyState(RiverDepth);
                }
            }

            /**
            * @brief Set bottom for wetlands, lakes, rivers
            * @param amount
            * @param type
            */
            void setExternalFlowBottom(double amount, FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    t_meter flowHead = getExternalFlowByName(type).getFlowHead();
                    double conduct = getExternalFlowByName(type).getConductance().value();
                    double RiverDepth = getExternalFlowByName(type).getRiverDepthSteadyState();
                    t_meter bottom = amount * si::meter;
                    double initConduct = getExternalFlowByName(type).getInitConductance().value();
                    removeExternalFlow(type);
                    if (flowHead.value() < bottom.value())
                        flowHead = bottom;
                    addExternalFlow(type, flowHead, conduct, bottom);
                    getExternalFlowByName(type).setInitConductance(initConduct);
                    getExternalFlowByName(type).setRiverDepthSteadyState(RiverDepth);
                }
            }

            /**
             * @brief Check for type river
             * @return bool
             */
            bool hasRiver() { return hasTypeOfExternalFlow(RIVER); }

            /**
            * @brief Check for type river
            * @return bool
            */
            bool hasRiver_MM() { return hasTypeOfExternalFlow(RIVER_MM); }

            /**
             * @brief Check for type GHB
             * @return bool
             */
            bool hasGHB() { return hasTypeOfExternalFlow(GENERAL_HEAD_BOUNDARY); }

            /**
             * @brief Get Q part (external sources) of flow equations
             * @return volume over time
             */
            t_vol_t getQ() noexcept {
                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me
                }
                t_vol_t eqFlow = getEqFlow();
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                for (const auto &flow : externalFlows) {
                    out += flow.second.getQ(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>();
                }
                return out;
            }

            /**
             * @brief Get Q part (external sources) source part of RHS of current density zone
             * @return volume over time
             */
            t_vol_t getQ(int localZetaID) noexcept {
                int densityZone = localZetaID; // the denity zone has the same ID as the zeta interface
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                t_vol_t Q = 0.0 * (si::cubic_meter / day);

                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = getHead(); // need the previous head
                t_meter flowHeight;
                t_meter bottomOfFlowInZone;
                t_meter flowHeightInZone;
                t_dim flowFractionOfZone;
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me cell has no special_flow
                }
                t_vol_t eqFlow = getEqFlow();
                for (const auto &[flowType, extFlow] : externalFlows) {
                    Q = extFlow.getQ(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>();
                    if (is(flowType).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                        // if this is zone 0
                        //if (densityZone == 0) { out += Q; }
                        out += Q * getFlowFractionOfZone(localZetaID, extFlow);
                    } else if (flowType == RECHARGE) {
                        // if recharge is into node (+) and this is zone 0
                        if ((Q.value() > 0 and densityZone == getSourceZoneRecharge())) {
                            out += Q;
                        }
                    } else if (flowType == GENERAL_HEAD_BOUNDARY) {
                        // if GHB flow is into node (gw_head < ghb_elevation) and ...
                        if (getHead() < extFlow.getFlowHead()) {
                            // ... this is the zone of GHB inflow
                            if (densityZone == getSourceZoneGHB()) {
                                out += Q;
                            }
                        // if GHB flow is out of node (gw_head > ghb_elevation) and ...
                        } else if (getHead() > extFlow.getFlowHead()) {
                            // ... interface above is not active (-> no zeta above active)
                            if (!isZetaActive(localZetaID - 1)) {
                                out += Q;
                            }
                        }
                    }
                }
                return out;
            }

            /**
             * @brief The effective porosity term for this node
             * @return square meter over time
             * @note In SSWI2: SWIHCOF (line 3744)
             */
            t_s_meter_t getEffectivePorosityTerm(){ // computed independent of steady or transient flow (see SWI2 doc "Using the SWI2 Package")
                t_s_meter_t out = (getEffectivePorosity() * get<t_s_meter, Area>()) /
                                  day;
                                  //(day * get<t_dim, StepSize>());
                NANChecker(out.value(), "getEffectivePorosityTerm");
                return out;
            }

            /**
             * @brief Dimensionless density (nus) at the top of this node
             * @return dimensionless
             * @note in SWI2: NUTOP, lines 1230-1249
             */
            t_dim getNusTop(){
                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                t_dim out = nusInZones.front();
                for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++){
                    if (isZetaAtTop(localZetaID)){
                        auto delnus = get<std::vector<t_dim>, Delnus>();
                        out += delnus[localZetaID];
                    }
                }
                return out;
            }

            /**
             * @brief Dimensionless density (nus) at the bottom of this node
             * @return dimensionless
             * @note in SWI2: NUBOT, lines 1230-1249
             */
            t_dim getNusBot(){
                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                t_dim out = nusInZones.back();
                for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++){
                    if (isZetaAtBottom(localZetaID)){
                        auto delnus = get<std::vector<t_dim>, Delnus>();
                        out -= delnus[localZetaID];
                    }
                }
                return out;
            }

            /**
             * @brief Checking if zeta surface is at node top or at head
             * @param localZetaID
             * @return bool
             */
            bool isZetaAtTop(int localZetaID){
                return (getZeta(localZetaID) >= getZetas().front());
            }

            /**
             * @brief Checking if zeta surface is at node bottom
             * @param localZetaID
             * @return bool
             * @note using ZetasTZero since they are the zetas from the beginning of the current iterative equation solving process
             */
            bool isZetaAtBottom(int localZetaID){
                return (getZeta(localZetaID) <= getZetas().back() or
                        getHead() < getBottom());
            }

            /**
             * @brief Checking if zeta surface between top and bottom
             * @param localZetaID
             * @return bool
             */
            bool isZetaBetween(int localZetaID){
                return (!isZetaAtTop(localZetaID) and !isZetaAtBottom(localZetaID));
            }

            /**
             * @brief Checking if any zeta surface is at node top or bottom
             * @return bool
             */
            bool isAnyZetaBetween(){
                for (int localZetaID = 0; localZetaID < getZetas().size(); ++localZetaID){
                    if (isZetaBetween(localZetaID)){ return true;}
                }
                return false;
            }

            /**
             * @brief Checking if zeta surface is active (between top and bottom and porosity above zero)
             * @param localZetaID
             * @return bool
             */
            bool isZetaActive(int localZetaID){
                return (isZetaBetween(localZetaID) and getEffectivePorosity().value() > 0);
            }

            /**
             * @brief Checking if zeta surface is the highest active zeta
             * @param localZetaID
             * @return bool
             */
            bool isThisTheHighestActiveZeta(int localZetaID){
                if (isZetaActive(localZetaID)) {
                    // if this zeta is active -> localZetaID >= 1
                    for (int i = localZetaID-1; i > 0; --i) {
                        if (isZetaActive(i)) { // if a higher zeta interface is active
                            return false;
                        }
                    }
                    // no higher active zeta interface exists
                    return true;
                } else {
                    // this zeta interface is not even active
                    return false;
                }
            }

            /**
             * @brief The source flow below a zeta surface (for the right hand side in the zeta equation)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note like G in SWI2 doc, but without vertical leakage; in SWI2 code: lines 3523-3569
             * G = RHS (of flow, for constant density) - HCOF_(i,j,k,n)*h^(m)_(i,j,k) + (verticalLeakage_(i,j,k-1,n) - verticalLeakage_(i,j,k,n))
             */
            t_vol_t getSources(int localZetaID){
                // We use the following sink/source concept:
                //  - the source zones of GHB are specified in config and/or in file, sink of GHB is the top zone
                //  - the source/sink zone of SWBs depend on interface heights
                //  - the source zone recharge is the freshwater zone
                // todo: compute BUFF with SSWI2_BDCH for constant head cells
                t_vol_t sources = 0.0 * (si::cubic_meter / day);

                // if the new groundwater head is above or equal to the node bottom
                if (getHead() >= getBottom()) { // lines 3532-3536
                    // todo add GhostNodeSources (apply getZoneFraction to it)
                    t_vol_t externalSources = - getQ(localZetaID) - getP(localZetaID) * getHead();

                    t_vol_t storageChange = getStorageCapacity() * (getHead() - getHead_TZero()) /
                                                    (day * get<t_dim, StepSize>());

                    sources = externalSources + storageChange * getZoneFraction(localZetaID); // see line 3535-3536
                }
                //if (sources.value() != 0){
                //    LOG(debug) << "getSources at node " << getID() << ": " << sources.value();
                //}
                return sources;
            }

            t_dim getZoneFraction(int localZetaID) {
                if (localZetaID >= getZetasTZero().size() -1) { throw "getZoneFraction(): localZetaID too high";}
                // return: height of this zone / height of all zones
                return (getZetaTZero(localZetaID) - getZetaTZero(localZetaID + 1)) /
                       (getZetasTZero().front() - getZetasTZero().back());
            }

            /**
             * @brief The pseudo source for a zeta surface (for the right hand side in the zeta equation)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note in SWI2 code lines 3574-3635, using SSWI2_SD and SSWI2_SR)
             */
            t_vol_t getPseudoSourceBelowZeta(int localZetaID) {
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<std::vector<t_dim>, Delnus>();

                // pseudo source term calculation (in 2 parts)
                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    // on right side of equation -> need to use zeta interface heights of previous time step
                    std::vector<t_meter> zoneThicknessesTZero = calculateZoneThicknessesTZero(neigPos, neigNodeID);
                    // calculating zone conductances for pseudo source term calculation
                    std::vector<t_s_meter_t> zoneConductancesTZero = getZoneConductances(neigPos, neigNodeID,
                                                                                        zoneThicknessesTZero);
                    t_s_meter_t zoneCondCumTZero = getZoneConductanceCum(localZetaID,zoneConductancesTZero);
                    if (isZetaActive(localZetaID) and // if IPLPOS == 0 (line 3571)
                        at(neigNodeID)->isZetaActive(localZetaID)) { // if neighbouring IPLPOS == 0 (e.g. line 2094)
                        //%% head part %%
                        t_vol_t head_part = -zoneCondCumTZero * (getAt<t_meter, Head>(neigNodeID) - getHead());
                        out += head_part;
                        //LOG(debug) << "head_part: " << head_part.value() << std::endl;

                        t_s_meter_t zoneCondCumTZeroDelnus;
                        for (int zetaID = 0; zetaID < getZetasTZero().size() - 1; zetaID++) {
                            //%% delnus part %%
                            if (zetaID < localZetaID) {
                                zoneCondCumTZeroDelnus = zoneCondCumTZero;
                            } else if (zetaID == localZetaID) {
                                continue;
                            } else if (zetaID > localZetaID) {
                                zoneCondCumTZeroDelnus = getZoneConductanceCum(zetaID, zoneConductancesTZero);
                            }
                            t_vol_t delnus_part = -delnus[zetaID] * zoneCondCumTZeroDelnus *
                                                  (at(neigNodeID)->getZetaTZero(zetaID) - getZetaTZero(zetaID));
                            out += delnus_part;
                            //LOG(debug) << "delnus_part (zetaID = " << zetaID << "): " << delnus_part.value() << std::endl;
                        }
                    }
                }
                NANChecker(out.value(), "getPseudoSourceBelowZeta");
                return out;
            }

            /**
             * @brief Flow at tips and toes
             * @param got neighbouring node
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             */
            t_vol_t getFluxHorizontal(NeighbourPosition neigPos, large_num neigNodeID, int localZetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<std::vector<t_dim>, Delnus>();

                // on right hand side of equation -> we need the zone conductances of previous time step
                std::vector<t_meter> zoneThicknessesTZero = calculateZoneThicknessesTZero(neigPos, neigNodeID);
                //LOG(debug) << "zoneThicknessesTZero: " << zoneThicknessesTZero[0].value() << ", " <<
                //                                          zoneThicknessesTZero[1].value();
                std::vector<t_s_meter_t> zoneConductancesTZero = getZoneConductances(neigPos, neigNodeID,
                                                                                zoneThicknessesTZero);

                t_s_meter_t zoneCondCumTZero = getZoneConductanceCum(localZetaID, zoneConductancesTZero);
                // %%head part %%
                t_vol_t head_part = zoneCondCumTZero * (getAt<t_meter, Head>(neigNodeID) - get<t_meter, Head>());
                out += head_part;
                //LOG(debug) << "head_part (tip/toe): " << head_part.value() << std::endl;

                // %%delnus part %%
                t_s_meter_t zoneCondCumTZeroDelnus;
                for (int zetaID = 0; zetaID < getZetas().size() - 1; zetaID++) {
                    if (zetaID <= localZetaID){
                        zoneCondCumTZeroDelnus = zoneCondCumTZero;
                    } else {
                        zoneCondCumTZeroDelnus = getZoneConductanceCum(zetaID,zoneConductancesTZero);
                    }
                    t_vol_t delnus_part = delnus[zetaID] * zoneCondCumTZeroDelnus *
                                          (at(neigNodeID)->getZetaTZero(zetaID) - getZetaTZero(zetaID));
                    out += delnus_part;
                    //LOG(debug) << "delnus_part (tip/toe): " << delnus_part.value() << std::endl;
                }
                return out;
            }

            /**
             * @brief Specification of boundary condition at tips and toes
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note adapted from SWI2 code lines 3637-3703 (includes usage of SSWI2_QR and SSWI2_QC)
             */
            t_vol_t getTipToeFlow(int localZetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);

                // tip and toe flow calculation
                if (isZetaActive(localZetaID)) { // if IPLPOS == 0 (line 3571)
                    for (auto const &[neigPos, neigNodeID]: neighbours) {
                        if (!at(neigNodeID)->isZetaActive(localZetaID)) { // if IPLPOS == 0 (line 3571)
                            if (neigPos == NeighbourPosition::TOP) {
                                // vertical leakage to TOP neighbour
                                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                                if (getFluxTop() < 0 * (si::cubic_meter / day) and
                                    at(neigNodeID)->getNusBot() >= nusInZones[localZetaID] and
                                    getNusBot() >= at(neigNodeID)->getNusBot()) { // IF ((qztop.LT.0).AND.(NUBOT(i,j,k-1).GE.NUS(iz)).AND.(NUBOT(j,i,k).GE.NUBOT(i,j,k-1))) THEN
                                    out += getFluxTop(); // in SWI2: qztop
                                }
                            } else if (neigPos == NeighbourPosition::DOWN) {
                                // vertical leakage to DOWN neighbour
                                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                                if(getFluxDown() < 0 * (si::cubic_meter / day) and
                                   at(neigNodeID)->getNusTop() < nusInZones[localZetaID] and
                                   getNusTop() <= at(neigNodeID)->getNusTop()) { // IF ((qzbot.LT.0).AND.(NUTOP(i,j,k+1).LT.NUS(iz)).AND.(NUTOP(j,i,k).LE.NUTOP(i,j,k+1))) THEN
                                    continue;
                                } else{
                                    out += getFluxDown(); // in SWI2: qzbot
                                }
                            } else { // if neighbour is horizontal
                                out -= getFluxHorizontal(neigPos, neigNodeID, localZetaID); // SSWI2_QR (left/right), SSWI2_QC (front/back)
                            }
                        }
                    }
                }
                return out;
            }

            /**
             * @brief The density zone conductance for variable density flow calculations (both flow and zeta equations)
             * @param got neighbouring node
             * @return vector with entries in square meter per time
             */
            std::vector<t_s_meter_t> getZoneConductances(NeighbourPosition neigPos, large_num neigNodeID,
                                                         std::vector<t_meter> zoneThicknesses) {
                std::vector<t_s_meter_t> out;
                t_s_meter_t conductance;
                t_s_meter_t zoneConductance;

                // calculate the zone thickness sum
                t_meter sumOfZoneThicknesses = 0 * si::meter;
                for (const auto &zoneThickness : zoneThicknesses) {
                    sumOfZoneThicknesses += zoneThickness;
;                }

                // calculate the density zone conductances
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; localZetaID++) {
                    //LOG(debug) << "zoneThicknesses[" << localZetaID << "]:" << zoneThicknesses[localZetaID].value();
                    zoneConductance = 0 * si::square_meter / day;
                    if (sumOfZoneThicknesses != (0 * si::meter)) { // adapted from SWI2 code line 1159
                        conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                        //LOG(debug) << "conductance:" << conductance.value();
                        zoneConductance = conductance * (zoneThicknesses[localZetaID] / sumOfZoneThicknesses);
                        //LOG(debug) << "zoneConductance[" << localZetaID << "] :" << zoneConductance.value();
                    }
                    else {
                        if (1 < localZetaID < (getZetas().size() - 2)) {
                            if(!isZetaActive(localZetaID) or // adapted from SWI2 code lines 1212-1222
                               !at(neigNodeID)->isZetaActive(localZetaID) or
                               !isZetaActive(localZetaID + 1) or
                               !at(neigNodeID)->isZetaActive(localZetaID + 1)) {
                                // this section is reached if getZetas().size() >= 4 and numberOfZones >= 3
                                zoneConductance = 0 * si::square_meter / day; // zone is inactive or across layers
                            }
                        }
                    }
                    //LOG(debug) << "zoneConductance[id = " << localZetaID << "] : " << zoneConductance.value() << std::endl;
                    out.push_back(zoneConductance);
                    NANChecker(zoneConductance.value(), "zoneConductance");
                }
                return out;
            }

            std::vector<t_meter> calculateZoneThicknesses(NeighbourPosition neigPos, large_num neigNodeID) {
                std::vector<t_meter> zoneThicknesses;
                t_meter zoneThickness;
                t_meter deltaZeta;
                t_meter deltaZeta_neig;
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; localZetaID++) {
                    deltaZeta = getZeta(localZetaID) - getZeta(localZetaID + 1);
                    deltaZeta_neig = at(neigNodeID)->getZeta(localZetaID) - at(neigNodeID)->getZeta(localZetaID + 1);
                    if (deltaZeta <= (0 * si::meter) or deltaZeta_neig <= (0 * si::meter)){ // adapted from SWI2 code line 1149
                        zoneThickness = 0 * si::meter;
                    } else {
                        zoneThickness = ((getLengthNeig(neigPos, neigNodeID) * deltaZeta) +
                                         (getNodeLength(neigPos) * deltaZeta_neig)) /
                                        (getLengthNeig(neigPos, neigNodeID) + getNodeLength(neigPos));
                    }
                    NANChecker(zoneThickness.value(), "zoneThickness");
                    zoneThicknesses.push_back(zoneThickness);
                }
                return zoneThicknesses;
            }

            std::vector<t_meter> calculateZoneThicknessesTZero(NeighbourPosition neigPos, large_num neigNodeID) {
                std::vector<t_meter> zoneThicknessesTZero;
                t_meter zoneThicknessTZero;
                t_meter deltaZeta;
                t_meter deltaZeta_neig;
                for (int localZetaID = 0; localZetaID < getZetasTZero().size() - 1; localZetaID++) {
                    deltaZeta = getZetaTZero(localZetaID) - getZetaTZero(localZetaID + 1);
                    deltaZeta_neig = at(neigNodeID)->getZetaTZero(localZetaID) -
                                     at(neigNodeID)->getZetaTZero(localZetaID + 1);
                    if (deltaZeta <= (0 * si::meter) or deltaZeta_neig <= (0 * si::meter)){ // adapted from SWI2 code line 1149
                        zoneThicknessTZero = 0 * si::meter;
                    } else {
                        zoneThicknessTZero = ((getLengthNeig(neigPos, neigNodeID) * deltaZeta) +
                                              (getNodeLength(neigPos) * deltaZeta_neig)) /
                                                        (getLengthNeig(neigPos, neigNodeID) + getNodeLength(neigPos));
                    }
                    NANChecker(zoneThicknessTZero.value(), "zoneThicknessTZero");
                    zoneThicknessesTZero.push_back(zoneThicknessTZero);
                }
                return zoneThicknessesTZero;
            }

            /**
             * @brief Sum of density zone conductances below current local zeta id
             * @param localZetaID zeta surface id in this node
             * @param zoneConductances density zone conductances
             * @return square meter per time
             */
            t_s_meter_t getZoneConductanceCum(int localZetaID, std::vector<t_s_meter_t> zoneConductances) {
                // calculate the sum of density zone conductances below a zeta surface n and add to vector out
                t_s_meter_t out = 0 * si::square_meter / day;

                for (int zoneID = localZetaID; zoneID < getZetas().size() - 1; zoneID++) {
                    out += zoneConductances[zoneID];
                };
                //LOG(debug) << "zoneConductanceCum:" << out.value() << std::endl;
                NANChecker(out.value(), "getZoneConductanceCum");
                return out;
            }

            /**
             * @brief The pseudo source for the flow equation, only used if variable density flow is active
             * @return volume per time
             * @note This accounts for the effects of variable density flow (in SWI2 code: QREXTRA/QFEXTRA)
             */
            t_vol_t getPseudoSourceNode() {
                t_vol_t out = 0.0 * (si::cubic_meter / day);

                auto delnus = get<std::vector<t_dim>, Delnus>();

                // pseudo source term calculation
                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    if (isAnyZetaBetween() or at(neigNodeID)->isAnyZetaBetween()) { // check if there are any active zeta surfaces
                        // in flow equation -> zeta interface heights are not updated between iterations
                        // -> here, calculateZoneThicknesses and calculateZoneThicknessesTZero give the same results
                        std::vector<t_meter> zoneThicknesses = calculateZoneThicknesses(neigPos, neigNodeID);
                        std::vector<t_s_meter_t> zoneConductances = getZoneConductances(neigPos, neigNodeID,
                                                                                        zoneThicknesses);
                        for (int zetaID = 0; zetaID < getZetas().size() - 1; zetaID++) {
                            t_s_meter_t zoneConductanceCum = getZoneConductanceCum(zetaID, zoneConductances);
                            t_vol_t pseudoSource = delnus[zetaID] * zoneConductanceCum *
                                                   (at(neigNodeID)->getZeta(zetaID) - getZeta(zetaID));
                            out -= pseudoSource;
                        }
                    }
                }
                //LOG(debug) << "getPseudoSourceNode: " << out.value() << std::endl;
                return out;
            }

            /**
             * @brief The flux correction term in vertical direction
             * @return volume per time
             * @note in SWI2 documentation: CV*BOUY; in code: QLEXTRA
             */
            t_vol_t getVerticalFluxCorrection(){
                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                t_meter headdiff = 0 * si::meter;
                t_vol_t out = 0 * (si::cubic_meter / day);

                if (neighbours.find(TOP) != neighbours.end()){ //Current node has a top node
                    large_num neigNodeID = neighbours[TOP];
                    // first part of the flux correction term
                    for (int localZetaID = 0; localZetaID < getZetasTZero().size() - 1; localZetaID++){
                        headdiff -= nusInZones[localZetaID] *
                                    (at(neigNodeID)->getZetaTZero(localZetaID + 1) - at(neigNodeID)->getZetaTZero(localZetaID));
                        // Note: in SWI2 documentation is, BOUY is calculated by adding headdiff (would be out +=),
                        // MODFLOW code for headdiff is as implemented here (with out -=)
                    }
                    // second part of the flux correction term (vertical conductance *  BOUY)
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(neigNodeID));
                    out = verticalConductance *
                          (headdiff +
                           0.5 * (at(neigNodeID)->getZetasTZero().back() - getZetasTZero().front()) *
                           (at(neigNodeID)->getNusBot() + getNusTop()));
                    // Note in SWI2 documentation, BOUY is calculated with a - between NUBOT and NUTOP,
                    // in MODFLOW code there is a + in the calculation of QLEXTRA
                    //LOG(debug) << "headdiff: " << headdiff.value() << std::endl;
                }
                return out;
            }

            t_vol_t getVerticalFluxCorrections(){
                t_vol_t out = 0 * (si::cubic_meter / day);

                for (auto const &[neigPos, neigNodeID]: vertical_neighbours) {
                    if (neigPos == NeighbourPosition::TOP) {
                        out -= getVerticalFluxCorrection();
                    }
                    if (neigPos == NeighbourPosition::DOWN) {
                        out += at(neigNodeID)->getVerticalFluxCorrection();
                    }
                }
                return out;
            }

            /**
             * @brief Calculates the upward vertical flux correction for variable density flow
             * @return volume per time
             * @note in SWI2 code: qztop
             */
            t_vol_t getFluxTop() {
                t_vol_t out = 0 * (si::cubic_meter / day);
                if (neighbours.find(TOP) != neighbours.end()) {
                    large_num neigNodeID = neighbours[TOP];
                    t_vol_t fluxFromTopNode = getVerticalFluxCorrection();
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(neigNodeID));
                    out = (verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(neigNodeID))) - fluxFromTopNode;
                    //LOG(debug) << "getFluxTop: " << out.value();
                }
                return out;
            }

            /**
             * @brief Calculates the downward vertical flux correction for variable density flow
             * @return volume per time
             * @note in SWI2 code: qzbot
             */
            t_vol_t getFluxDown() {
                t_vol_t out = 0 * (si::cubic_meter / day);
                if (neighbours.find(DOWN) != neighbours.end()) {
                    large_num neigNodeID = neighbours[DOWN];
                    t_vol_t fluxFromDownNode = at(neigNodeID)->getVerticalFluxCorrection();
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(neigNodeID));
                    out = (verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(neigNodeID))) + fluxFromDownNode;
                    //LOG(debug) << "getFluxDown: " << out.value();
                }
                return out;
            }

            /**
             * @brief Updating top zeta surface height after the flow/head solution is found
             * @note Top zeta surface height is set to the new groundwater head. In SWI2: SSWI2_UPZ1
             */
            void setTopZetaToHead(){
                t_meter head = getHead();
                t_meter bottomOfNode = getBottom();
                t_meter topOfNode = getElevation();
                t_meter newHeight;
                // if groundwater head is ABOVE the top of the node
                if (head > topOfNode) {
                    // clip zeta to the top of the node
                    setZeta(0, topOfNode);
                // if groundwater head is BELOW the top of the node
                } else {
                    // if groundwater head is ABOVE the bottom of the node
                    if (head > bottomOfNode) {
                        newHeight = head;
                        // if groundwater head is BELOW OR EQUAL to the bottom of the node
                    } else { // head <= bottomOfNode
                        newHeight = bottomOfNode;
                    }
                    // update the first zeta surface
                    setZeta(0, newHeight);
                }
                // update all other zeta surfaces that are ABOVE the updated first zeta surface
                if (getZetas().size() > 1) {
                    for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                        if (getZeta(localZetaID) > getZeta(0)) {
                            setZeta(localZetaID, getZeta(0));
                        }
                    }
                }
                setZetasTZero(); // set ZetasTZero = Zetas
            }

            /**
             * @brief Vertical movement of zeta surfaces through top of this node. This function is required to move a
             * zeta surface with the height equal to the top of a lower node and bottom to an upper node. Without this
             * function, that zeta surface would be stuck there since the other only consider "active" surfaces (which
             * are between the top and bottom of a node)
             * @note in SWI2 code: SSWI2_VERTMOVE
             */
            void zetaMovementBetweenLayers() {
                // skip nodes where head is below the bottom of the node
                t_vol_t fluxCorrectionTop; // in SWI2: qztop
                t_s_meter_t verticalConductanceTop;
                t_meter deltaZeta;
                if (neighbours.find(TOP) != neighbours.end()) {
                    large_num topNodeID = neighbours[TOP];
                    // if head is above the bottom of the node, both in this node AND in the node above (SWI2 line 2325)
                    if (get<t_meter, Head>() >= getBottom() and
                        getAt<t_meter, Head>(topNodeID) >=
                                (getAt<t_meter, Elevation>(topNodeID) - getAt<t_meter, VerticalSize>(topNodeID))) {

                        for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                            // zeta only moves through the top of a node if there is a ZETA surface
                            // - at the top of the current node (in SWI2: IPLPOS_(i,j,k,n) = 1)
                            // - AND at the bottom of the top node (in SWI2: IPLPOS_(i,j,k-1,n) = 2)
                            if (isZetaAtTop(localZetaID) and at(topNodeID)->isZetaAtBottom(localZetaID)) {

                                // calculate flux through the top
                                fluxCorrectionTop = getFluxTop();
                                //LOG(debug) << "fluxCorrectionTop: " << fluxCorrectionTop.value() << std::endl;

                                // if vertical flux through the top of the node is positive...
                                if (fluxCorrectionTop.value() > 0 and at(topNodeID)->getEffectivePorosity().value() > 0) {
                                    deltaZeta = (fluxCorrectionTop * (day)) /
                                                (get<t_s_meter, Area>() * getAt<t_dim, EffectivePorosity>(topNodeID)); // * get<t_dim, StepSize>()
                                    // ...lift zeta height of the lowest zeta surface in top node
                                    t_meter zeta_back_top = at(topNodeID)->getZetas().back();
                                    at(topNodeID)->setZeta(localZetaID, zeta_back_top + deltaZeta);

                                    // if vertical flux through the top of the node is negative...
                                } else if (fluxCorrectionTop.value() < 0 and getEffectivePorosity().value() > 0) {
                                    deltaZeta = (fluxCorrectionTop * (day)) /
                                                (get<t_s_meter, Area>() * getEffectivePorosity()); // * get<t_dim, StepSize>()
                                    //LOG(debug) << "deltaZeta: " << deltaZeta.value() << std::endl;
                                    // ...lower zeta height of this zeta surface by delta zeta
                                    setZeta(localZetaID, getZetas().front() + deltaZeta);
                                }
                            }
                        }
                    }
                }
            }

            /**
             * @brief horizontal tracking of zeta surface tips and toes (at top and bottom of the current node)
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 documentation: "Tip and Toe Tracking", in SWI2 code: SSWI2_HORZMOVE
             */
            void horizontalZetaMovement(){
                t_meter maxDelta; // maximum "allowed" height difference of a zeta surface n between adjacent nodes
                t_meter delta_self; // potential zeta height adjustment for this node
                t_meter delta_neig; // potential zeta height adjustment for neighbour node
                t_meter delta_opp; // potential zeta height adjustment for opposite neighbour node
                if (getHead() < getBottom()) { return; }
                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                        if (isZetaActive(localZetaID)) {
                            // get max delta of zeta between nodes
                            if (at(neigNodeID)->isZetaAtBottom(localZetaID)) {
                                maxDelta = 0.5 * (getNodeLength(neigPos) + getLengthNeig(neigPos, neigNodeID)) * get<t_dim, MaxToeSlope>();
                            } else if (at(neigNodeID)->isZetaAtTop(localZetaID)) {
                                maxDelta = 0.5 * (getNodeLength(neigPos) + getLengthNeig(neigPos, neigNodeID)) * get<t_dim, MaxTipSlope>();
                            }
                            //LOG(debug) << "maxDelta: " << maxDelta.value() << std::endl;

                            // if tracking tip/toe: raise/lower this zeta surface in this node by:
                            delta_self = get<t_dim, SlopeAdjFactor>() * maxDelta *
                                         ((at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID)) /
                                          ((getEffectivePorosity() * getNodeLength(neigPos)) +
                                           (at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID))));
                            // if tracking tip/toe: lower/raise this zeta surface in neighbouring node by:
                            delta_neig = get<t_dim, SlopeAdjFactor>() * maxDelta *
                                         ((getEffectivePorosity() * getNodeLength(neigPos)) /
                                          ((getEffectivePorosity() * getNodeLength(neigPos)) +
                                           (at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID))));

                            if (at(neigNodeID)->isZetaAtBottom(localZetaID)) {
                                //%% Toe tracking %%
                                t_meter zetaDif = getZeta(localZetaID) - at(neigNodeID)->getZetas().back();
                                //LOG(debug) << "zetaDif (toe): " << zetaDif.value() << std::endl;
                                if (zetaDif > maxDelta) {
                                    setZeta(localZetaID, getZeta(localZetaID) - delta_self);
                                    //LOG(debug) << "delta_self (toe): " << delta_self.value() << std::endl;
                                    t_meter zeta_back_neig = at(neigNodeID)->getZetas().back();
                                    at(neigNodeID)->setZeta(localZetaID, zeta_back_neig + delta_neig);
                                    //LOG(debug) << "delta_neig (toe): " << delta_neig.value() << std::endl;
                                }
                            } else if (at(neigNodeID)->isZetaAtTop(localZetaID)) {
                                //%% Tip tracking %%
                                t_meter zetaDif = at(neigNodeID)->getZetas().front() - getZeta(localZetaID);
                                //LOG(debug) << "zetaDif (tip): " << zetaDif.value() << std::endl;
                                if (zetaDif > maxDelta) {
                                    setZeta(localZetaID, getZeta(localZetaID) + delta_self);
                                    t_meter zeta_front_neig = at(neigNodeID)->getZetas().front();
                                    at(neigNodeID)->setZeta(localZetaID, zeta_front_neig - delta_neig);
                                }
                            }

                            if ((getZeta(localZetaID) - getZetas().back()) < (get<t_dim, MinDepthFactor>() * delta_neig)) {
                                NeighbourPosition oppNeigPos = getOppositePosition(neigPos);
                                if (neighbours.find(getOppositePosition(neigPos)) != neighbours.end()){
                                    large_num nodeIDOppNeig = neighbours[getOppositePosition(neigPos)];
                                    if (at(nodeIDOppNeig)->isZetaActive(localZetaID)) {
                                        // change zeta in other direction neighbour
                                        delta_opp = ((getZeta(localZetaID) - getZetas().back()) *
                                                     (getNodeLength(oppNeigPos) * getEffectivePorosity()) /
                                                     (getLengthNeig(oppNeigPos, nodeIDOppNeig) *
                                                      at(nodeIDOppNeig)->getEffectivePorosity()));
                                        t_meter zeta_opp = at(nodeIDOppNeig)->getZeta(localZetaID);
                                        at(nodeIDOppNeig)->setZeta(localZetaID, zeta_opp + delta_opp);
                                        setZeta(localZetaID, getZetas().back());
                                        //LOG(debug) << "delta_opp (toe): " << delta_opp.value() << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            NeighbourPosition getOppositePosition(NeighbourPosition position) {
                if (position == NeighbourPosition::BACK or position == NeighbourPosition::BACKBACK or
                    position == NeighbourPosition::BACKLEFT or position == NeighbourPosition::BACKRIGHT or
                    position == NeighbourPosition::BACKBACKLEFT or position == NeighbourPosition::BACKBACKRIGHT) {
                    return NeighbourPosition::FRONT;
                } else if (position == NeighbourPosition::FRONT or position == NeighbourPosition::FRONTFRONT or
                           position == NeighbourPosition::FRONTLEFT or position == NeighbourPosition::FRONTRIGHT or
                           position == NeighbourPosition::FRONTFRONTLEFT or position == NeighbourPosition::FRONTFRONTRIGHT) {
                    return NeighbourPosition::BACK;
                } else if (position == NeighbourPosition::LEFT or position == NeighbourPosition::LEFTLEFT or
                           position == NeighbourPosition::LEFTFRONT or position == NeighbourPosition::LEFTBACK or
                           position == NeighbourPosition::LEFTLEFTFRONT or position == NeighbourPosition::LEFTLEFTBACK) {
                    return NeighbourPosition::RIGHT;
                } else if (position == NeighbourPosition::RIGHT or position == NeighbourPosition::RIGHTRIGHT or
                           position == NeighbourPosition::RIGHTFRONT or position == NeighbourPosition::RIGHTBACK or
                            position == NeighbourPosition::RIGHTRIGHTFRONT or position == NeighbourPosition::RIGHTRIGHTBACK) {
                    return NeighbourPosition::LEFT;
                } else {
                    throw "Position unavailable for function getOppositePosition";
                }
            }

            /**
             * @brief clip zeta surfaces that are outside of the bounds (bounds: first and last zeta surface)
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ZETACLIP
             */
            void clipInnerZetas() {
                for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                    if (isZetaActive(localZetaID)) {
                        if (getZeta(localZetaID) < getZetas().back()) { setZeta(localZetaID, getZetas().back()); }
                        if (getZeta(localZetaID) > getZetas().front()) { setZeta(localZetaID, getZetas().front()); }
                    }
                }
            }

            /**
             * @brief Adjust zeta surfaces if they have crossed or their height difference is below a minimum value
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ZETACROSS
             *  localZetaID | x       | n
             *  1           | -       | -
             *  2           | 1       | 1,2,3
             *  3           | 1,2     | 2,3,4; 1,2,3,4
             *  4           | 1,2,3   | 3,4,5; 2,3,4,5; 1,2,3,4,5
             */
            void correctCrossingZetas(){
                t_meter zetaDifferenceCap = 0.001 * si::meter; // todo move to config
                t_meter zetaAverage;
                t_meter zetaSum;
                int zetaRef;
                int counter;
                for (int localZetaID = 1; localZetaID < getZetas().size() - 2; ++localZetaID) {
                    // if zeta surface n is very close to or lower than the zeta surface that SHOULD be below (n+1)
                    if (getZeta(localZetaID) - getZeta(localZetaID + 1) < zetaDifferenceCap) {
                        // make the zeta height of both surfaces their average
                        zetaAverage = 0.5 * (getZeta(localZetaID) + getZeta(localZetaID + 1));
                        setZeta(localZetaID, zetaAverage);
                        setZeta(localZetaID + 1, zetaAverage);
                        // if there are zeta surfaces above that could possibly be crossing
                        for (int x = 1; x < localZetaID; ++x) {
                            zetaSum = 0 * si::meter;
                            // create a range from n_min (n-x) to n_max (n+1)
                            zetaRef = localZetaID - x;
                            // if surface zetaRef crosses (or is very close to) the zeta surface at localZetaID+1
                            if (getZeta(zetaRef) - getZeta(localZetaID + 1) < zetaDifferenceCap) {
                                // sum heights from zetaRef surface to localZetaID+1
                                counter = 0;
                                for (int n = zetaRef; n <= localZetaID + 1; ++n) {
                                    zetaSum += getZeta(n);
                                    ++counter;
                                }
                                if (counter != 0){
                                    // calculate the average height from surface zetaRef to localZetaID+1
                                    zetaAverage = zetaSum / (counter * si::si_dimensionless);
                                    // set all zeta surfaces from zetaRef to localZetaID+1 to that average
                                    for (int n = zetaRef; n <= localZetaID + 1; ++n) { setZeta(n, zetaAverage); }
                                }
                            }
                        }
                    }
                }
            }

            /**
             * @brief Avoid zeta surface locking
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ANTILOCKMIN
             */
            void preventZetaLocking(){
                t_meter maxDelta;
                t_meter maxAdjustment_self;
                t_meter maxAdjustment_neig;
                t_meter delta_self;
                t_meter delta_neig;
                t_meter zeta_neig;
                for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                    if (isZetaAtTop(localZetaID) or isZetaAtBottom(localZetaID)) {
                        // iterate through horizontal neighbours
                        for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {

                            // determine max delta zeta
                            maxAdjustment_self = get<t_meter, VerticalSize>() * get<t_dim, SlopeAdjFactor>();
                            maxAdjustment_neig = getAt<t_meter, VerticalSize>(neigNodeID) * get<t_dim, SlopeAdjFactor>();
                            if (get<t_meter, VDFLock>() > maxAdjustment_self or
                                get<t_meter, VDFLock>() > maxAdjustment_neig) {
                                maxDelta = std::min(maxAdjustment_self, maxAdjustment_neig);
                            } else {
                                maxDelta = get<t_meter, VDFLock>();
                            }

                            if (getEffectivePorosity().value() == 0 and
                                at(neigNodeID)->getEffectivePorosity().value() == 0){
                                continue;
                            }
                            // calculate delta zeta of node and neighbour
                            delta_self = maxDelta * (at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID)) /
                                         (getEffectivePorosity() * getNodeLength(neigPos) +
                                          at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID));
                            delta_neig = maxDelta * (getEffectivePorosity() * getNodeLength(neigPos)) /
                                         (getEffectivePorosity() * getNodeLength(neigPos) +
                                          at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID));

                            // if a zeta surface is at the BOTTOM of this node and at the TOP of the neighbour
                            // else if a zeta surface is at the TOP of this node and at the BOTTOM of the neighbour
                            if (isZetaAtBottom(localZetaID) && // IPLPOS_self = 2
                                at(neigNodeID)->isZetaAtTop(localZetaID)) { // IPLPOS_neig = 1
                                // adjust zeta surface heights
                                setZeta(localZetaID, getZeta(localZetaID) + delta_self);
                                zeta_neig = at(neigNodeID)->getZeta(localZetaID);
                                at(neigNodeID)->setZeta(localZetaID, zeta_neig - delta_neig);
                            } else if (isZetaAtTop(localZetaID) && // IPLPOS_self = 1
                                       at(neigNodeID)->isZetaAtBottom(localZetaID)) {// IPLPOS_neig = 2
                                // adjust zeta surface heights
                                setZeta(localZetaID, getZeta(localZetaID) - delta_self);
                                zeta_neig = at(neigNodeID)->getZeta(localZetaID);
                                at(neigNodeID)->setZeta(localZetaID, zeta_neig + delta_neig);
                            }
                        }
                    }
                }
            }

            /**
             * @brief The length of the neighbouring node, in direction to the neighbouring node
             * @param got neighbour node
             * @return meter
             */
            t_meter getLengthNeig(NeighbourPosition neigPos, large_num neigNodeID){
                if ( isLeftOrRight(neigPos) ) {
                    return getAt<t_meter, EdgeLengthLeftRight>(neigNodeID);
                } else if ( isFrontOrBack(neigPos) ) {
                    return getAt<t_meter, EdgeLengthFrontBack>(neigNodeID);
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            /**
             * @brief The length of this node, in direction to the neighbouring node
             * @param got neighbour node
             * @return meter
             */
            t_meter getNodeLength(NeighbourPosition neigPos){
                if ( isLeftOrRight(neigPos) ){
                    return get<t_meter, EdgeLengthLeftRight>();
                } else if ( isFrontOrBack(neigPos) ){
                    return get<t_meter, EdgeLengthFrontBack>();
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            /**
             * @brief The width of this node, perpendicular to the direction to the neighbouring node
             * @param got neighbour node
             * @return meter
             */
            t_meter getNodeWidth(NeighbourPosition neigPos) {
                if ( isLeftOrRight(neigPos) ){
                    return get<t_meter, EdgeLengthFrontBack>();
                } else if ( isFrontOrBack(neigPos)  ){
                    return get<t_meter, EdgeLengthLeftRight>();
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            bool isLeftOrRight(NeighbourPosition neig) {
                std::unordered_map<NeighbourPosition, bool> leftRightMap = {
                        {LEFT, true}, {RIGHT, true},
                        {LEFTFRONT, true}, {LEFTBACK, true}, {RIGHTFRONT, true}, {RIGHTBACK, true},
                        {LEFTLEFTFRONT, true}, {LEFTLEFTBACK, true}, {RIGHTRIGHTFRONT, true}, {RIGHTRIGHTBACK, true},
                        {LEFTLEFT, true}, {RIGHTRIGHT, true}
                };
                try{
                    return leftRightMap.at(neig);
                } catch (const std::out_of_range &e) { return false;}
                //return (neig == LEFT or neig == RIGHT or
                //        neig == LEFTFRONT or neig == LEFTBACK or neig == RIGHTFRONT or neig == RIGHTBACK or);
            }

            bool isFrontOrBack(NeighbourPosition neig) {
                std::unordered_map<NeighbourPosition, bool> frontBackMap = {
                        {FRONT, true}, {BACK, true},
                        {FRONTLEFT, true}, {FRONTRIGHT, true}, {BACKLEFT, true}, {BACKRIGHT, true},
                        {FRONTFRONTLEFT, true}, {FRONTFRONTRIGHT, true}, {BACKBACKLEFT, true}, {BACKBACKRIGHT, true},
                        {FRONTFRONT, true}, {BACKBACK, true}
                };

                try{
                    return frontBackMap.at(neig);
                } catch (const std::out_of_range &e) { return false;}
                //return (neig == FRONT or neig == BACK or
                //       neig == FRONTLEFT or neig == FRONTRIGHT or neig == BACKLEFT or neig == BACKRIGHT);
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
                t_vol_t eqFlow = getEqFlow();
                for (const auto &flow : externalFlows) {
                    if (is(flow.second.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                        // if flow depends on head
                        if (flow.second.isFlowHeadDependent(head)) {
                            out += flow.second.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>();
                        }
                    } else { // RECHARGE, GENERAL_HEAD_BOUNDARY
                        out += flow.second.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>();
                    }
                }
                return out;
            }

            /**
             * @brief Get P part of source part of RHS of current density zone
             * @return volume over time
             */
            t_s_meter_t getP(int localZetaID) noexcept {
                int densityZone = localZetaID; // the denity zone has the same ID as the zeta interface
                t_s_meter_t out = 0.0 * (si::square_meter / day);
                t_s_meter_t P = 0.0 * (si::square_meter / day);

                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>(); // need the current head
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me
                }
                t_vol_t eqFlow = getEqFlow();
                for (const auto &[flowType, extFlow] : externalFlows) { // todo at all loops through external flows
                    if (is(flowType).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                        if (extFlow.isFlowHeadDependent(head)) {
                            P = extFlow.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>();
                            // if this is zone 0
                            //if (densityZone == 0) { out += P; }
                            out += P * getFlowFractionOfZone(localZetaID, extFlow);
                        }
                    } else {
                        if (flowType == RECHARGE) { // getP for recharge is 0
                            continue;
                        } else if (flowType == GENERAL_HEAD_BOUNDARY) {
                            P = extFlow.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>();
                            // if GHB flow is into node (gw_head < ghb_elevation) and ...
                            if (getHead() < extFlow.getFlowHead()) {
                                // ... this is the zone of GHB inflow
                                if (densityZone == getSourceZoneGHB()) {
                                    out += P;
                                }
                            // if GHB flow is out of node (gw_head > ghb_elevation) and ...
                            } else if (getHead() > extFlow.getFlowHead()) {
                                // ... interface above is inactive (-> no zeta above active)
                                if (!isZetaActive(localZetaID-1)) {
                                    out += P;
                                }
                            }
                            //LOG(debug) << "GHB P at node " << getID() << " is " << P.value();
                        }
                    }
                }
                return out;
            }


            double getFlowFractionOfZone(int localZetaID, ExternalFlow extFlow) {
                t_meter bottomOfFlowInZone;
                t_meter flowHeightInZone;

                // if gw_head is below the flow head (-> flow is into node)
                if (getHead() < extFlow.getFlowHead()) {
                    // todo currently, rivers provide freshwater only.
                    //  If they should be capable of adding saline water, GHB elevation should be taken into account
                    return 0; // rivers provide freshwater only, hence do not add to any density zone
                    // gw_head is above flow head (-> flow is out of node)
                } else {
                    // height of flow is: gw_head - flow_bottom
                    t_meter flowHeight = getHead() - extFlow.getBottom();

                    // if this zeta is below the flow bottom
                    if (getZetaTZero(localZetaID) < extFlow.getBottom()) {
                        // no water from its zone enters the SWB
                        return 0;
                        // this zeta is above the flow bottom
                    } else {
                        // if the zeta below is higher than the river bottom ...
                        if (getZetaTZero(localZetaID + 1) > extFlow.getBottom()) {
                            // ... the bottom of flow into/out of this zone is the lower zeta
                            bottomOfFlowInZone = getZetaTZero(localZetaID + 1);
                            // the zeta below is not higher than the river bottom, ...
                        } else {
                            // ... so the bottom of flow into/out of this zone is the flow bottom
                            bottomOfFlowInZone = extFlow.getBottom();
                        }
                        flowHeightInZone = getZetaTZero(localZetaID) - bottomOfFlowInZone;
                        if (flowHeightInZone.value() < 0) {
                            LOG(userinfo)
                                << "getFlowFractionOfZone: flowHeightInZone is below 0. Should never happen!";
                            return 0;
                        }
                        return flowHeightInZone / flowHeight;
                    }
                }
            }

            /**
             * @brief Get flow which is not groundwater head dependent
             * @return volume over time
             * Flow can be added to constant flows on right side of the equations
             * If head is above river bottom for example
             */
            t_vol_t calculateNotHeadDependentFlows() noexcept {
                t_meter eq_head = get<t_meter, EQHead>();
                t_meter head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                try {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                } catch (const std::out_of_range &e) {//ignore me
                }
                t_vol_t eqFlow = getEqFlow();
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                //Q part is already subtracted in RHS
                for (const auto &flow : externalFlows) {
                    if (is(flow.second.getType()).in(RIVER, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND, DRAIN)) {
                        if (not flow.second.isFlowHeadDependent(head)) {
                            out += flow.second.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepSize>() *
                                   flow.second.getBottom();
                        }
                    }
                }
                return out;
            }

            /**
             * @brief The matrix entry for the cell
             * @return map <CellID,Conductance>
             * The left hand side of the equation
             */
            std::unordered_map<large_num, t_s_meter_t> getMatrixEntries() {
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(neighbours.size()+1);
                t_s_meter_t conduct;

                //Get all conductances from neighbouring cells
                for (const auto &neighbour: neighbours) {
                    auto neigPos = neighbour.first;
                    auto neigNodeID = neighbour.second;
                    //LOG(debug) << "neigNodeID: " << neigNodeID;
                    conduct = 0 * si::square_meter / day;
                    if (neigPos == TOP or neigPos == DOWN) {
                        conduct = mechanics.calculateVerticalConductance(createDataTuple(neigNodeID));
                        //LOG(debug) << "vertical conductance: " << conduct.value();
                    } else {
                        if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                            conduct = mechanics.calculateEFoldingConductance(createDataTuple<Head>(neigPos, neigNodeID),
                                                                             get<t_meter, EFolding>(),
                                                                             getAt<t_meter, EFolding>(neigNodeID));
                            //LOG(debug) << "horizontal conductance (using e-folding): " << conduct.value();
                        } else {
                            conduct = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                            //LOG(debug) << "horizontal conductance: " << conduct.value();
                        }
                    }
                    NANChecker(conduct.value(), "Conductances");
                    // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                    // (used to solve for the head at the neighbouring node)
                    out[neigNodeID] = std::move(conduct);
                }

                t_s_meter_t conductNode = 0 * si::square_meter / day;

                // To solve for the head at this node, the conductances to neighbours and HCOF are used
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &[neigNodeID, conductNeig]: out) { conductNode -= conductNeig; }
                // add HCOF
                t_s_meter_t hcof = getP() - (getStorageCapacity() / (day * get<t_dim, StepSize>()));
                conductNode += hcof;
                //LOG(debug) << "conductNode: " << conductNode.value();
                NANChecker(conductNode.value(), "conductNode");

                // add resulting conductance to solve for the head at this node to out
                out[getID()] = std::move(conductNode);
                return out;
            };

            /**
             * @brief The matrix entry for the left hand side of the zeta surface equation
             * @param localZetaID zeta surface id in this node
             * @return map <CellID,Conductance>
             */
            std::unordered_map<large_num, t_s_meter_t> getMatrixEntries(int localZetaID) {
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(horizontal_neighbours.size()+1);
                auto delnus = get<std::vector<t_dim>, Delnus>();
                std::vector<t_s_meter_t> zoneConductances;
                t_s_meter_t zoneConductanceCum;
                t_s_meter_t zetaMovementConductance;

                // if this zeta interface is not active: return directly
                if (!isZetaActive(localZetaID)){ return out; }

                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    zetaMovementConductance = 0 * (si::square_meter / day);
                    if (at(neigNodeID)->isZetaActive(localZetaID)) {
                        // on left hand side of equation -> need updated zeta heights
                        std::vector<t_meter> zoneThicknesses = calculateZoneThicknesses(neigPos, neigNodeID);
                        zoneConductances = getZoneConductances(neigPos, neigNodeID, zoneThicknesses);
                        zoneConductanceCum = getZoneConductanceCum(localZetaID, zoneConductances);
                        zetaMovementConductance += delnus[localZetaID] * zoneConductanceCum; // in SWI2: SWISOLCC/R
                        //LOG(debug) << "zoneConductanceCum = " << zoneConductanceCum.value() << std::endl;
                        NANChecker(zetaMovementConductance.value(), "zetaMovementConductance");
                        // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                        out[neigNodeID] = zetaMovementConductance;
                    }
                }

                // To solve for zeta in this node, the conductances to neighbours and porosity term are used
                t_s_meter_t conductNode = 0 * (si::square_meter / day); // SWI_HCOF
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &c: out) { conductNode -= c.second; }
                // subtract effective porosity term
                conductNode -= getEffectivePorosityTerm(); // subtracting effective porosity term
                //LOG(debug) << "effectivePorosityTerm = " << getEffectivePorosityTerm().value() << std::endl;

                // add conductance of this node to out, the key in the unordered map is the ID of this node
                out[get<large_num, ID>()] = conductNode;
                return out;
            };

            /**
             * @brief The right hand side of the flow equation
             * @return volume per time
             */
            t_vol_t getRHS() {
                t_vol_t externalSources = - getQ() - calculateNotHeadDependentFlows(); // e.g., recharge, river, lakes, wetlands
                //LOG(debug) << "externalSources: " << externalSources.value() << std::endl;
                t_vol_t dewateredFlow = calculateDewateredFlow(); // only if node has bottom neighbour
                //LOG(debug) << "dewateredFlow: " << dewateredFlow.value() << std::endl;
                t_vol_t storageFlow = -getStorageCapacity() * getHead_TZero() / (day * get<t_dim, StepSize>());
                //LOG(debug) << "storageFlow: " << storageFlow.value() << std::endl;
                t_vol_t gncFromNodes = getGNCFromNodes();
                //LOG(debug) << "gncFromNodes: " << gncFromNodes.value() << std::endl;
                t_vol_t gncToRefined = getGNCToRefinedNode();
                //LOG(debug) << "gncToRefined: " << gncToRefined.value() << std::endl;
                t_vol_t internalSources = dewateredFlow + storageFlow + gncFromNodes + gncToRefined;
                t_vol_t out = externalSources + internalSources;
                //LOG(debug) << "RHS constant density: " << out.value() << std::endl;
                NANChecker(out.value(), "RHS constant density");

                if (get<bool, IsDensityVariable>()) {
                    // save constant density RHS (without variable density terms) for calculation of zeta movement
                    //set<t_vol_t, InternalSources>(internalSources);

                    // calculate Pseudo-Source Flow
                    t_vol_t pseudoSourceNode = getPseudoSourceNode();
                    // calculate Vertical Flux Correction (from top neighbour)
                    t_vol_t verticalFluxCorrections = getVerticalFluxCorrections();
                    out += pseudoSourceNode + verticalFluxCorrections;
                }
                //LOG(debug) << "getRHS (nodeID: " << getID() << "): " << out.value() << std::endl;
                NANChecker(out.value(), "RHS");
                return out;
            }

            /**
             * @brief calculate the right hand side for zeta surface equation (b_zetas)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             */
            t_vol_t getRHS(int localZetaID){
                if (!isZetaActive(localZetaID)) { throw "getRHS(int localzetaID) was called for inactive zeta";} // (line 3571)

                t_vol_t porosityTerm = - getEffectivePorosityTerm() * getZetaTZero(localZetaID);
                //LOG(debug) << "porosityTerm: " << porosityTerm.value();
                t_vol_t pseudoSourceBelowZeta = getPseudoSourceBelowZeta(localZetaID); // in SWI2 code: SSWI2_SD and SSWI2_SR
                //LOG(debug) << "pseudoSourceBelowZeta: " << pseudoSourceBelowZeta.value() << std::endl;
                t_vol_t sources = - getSources(localZetaID); // in SWI2 code: part of BRHS; in SWI2 doc: G or known source term below zeta
                //LOG(debug) << "sources: " << sources.value();
                t_vol_t tipToeFlow = getTipToeFlow(localZetaID); // in SWI2 code: SSWI2_QR and SSWI2_QC
                //LOG(debug) << "tipToeFlow: " << tipToeFlow.value();

                t_vol_t out = porosityTerm + sources + tipToeFlow + pseudoSourceBelowZeta;
                //LOG(debug) << "out: " << out.value();

                NANChecker(out.value(), "getRHS(int localZetaID)");
                return out;
            }

            void setHeadAndHeadChange(t_meter change) noexcept {
                NANChecker(change.value(), "Set Head Change at nodeID = " + std::to_string(getID()));
                set<t_meter, HeadChange>(change);
                setHead(getHead() + change);
            }

            t_meter calcInitialHead(t_meter initialParam) noexcept { return __calcInitialHead(initialParam); }

            bool isStaticNode() noexcept { return __isStaticNode(); }

            PhysicalProperties &getProperties() { return fields; } // Question: rename to getNodeProperties?

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
                t_vol_t dewateredFlow = -calculateDewateredFlow();
                if (compare(dewateredFlow.value())) {
                    out += boost::units::abs(dewateredFlow);
                }
                return out;
            }

            /**
             * @brief Calculate the lateral flow velocity
             * @param pos
             * @return
             */
            quantity<Velocity> getVelocity(NeighbourPosition neigPos, large_num neigNodeID) {
                t_vol_t lateral_flow{0 * si::cubic_meter / day};
                auto vert_size = get<t_meter, VerticalSize>();
                t_s_meter_t conductance;
                if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                    conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(neigPos, neigNodeID),
                                                                                               get<t_meter, EFolding>(),
                                                                                               getAt<t_meter, EFolding>(neigNodeID));
                } else {
                    conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                }

                lateral_flow = conductance * (get<t_meter, Head>() - getAt<t_meter, Head>(neigNodeID));
                return lateral_flow / (vert_size * vert_size);
            }

            /**
             * @brief Calculate flow velocity for flow tracking
             * Vx and Vy represent the flow velocity in x and y direction.
             * A negative value represents a flow in the opposite direction.
             * @return Velocity vector (x,y)
             */
            std::pair<double, double> getVelocityVector() {
                quantity<Velocity> Vx{};
                quantity<Velocity> Vy{};

                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    if (neigPos == NeighbourPosition::LEFT) {
                        Vx += -getVelocity(neigPos, neigNodeID);
                    }
                    if (neigPos == NeighbourPosition::BACK) {
                        Vy += -getVelocity(neigPos, neigNodeID);
                    }

                    if (neigPos == NeighbourPosition::RIGHT) {
                        Vx += getVelocity(neigPos, neigNodeID);
                    }
                    if (neigPos == NeighbourPosition::FRONT) {
                        Vy += getVelocity(neigPos, neigNodeID);
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
                         large_num SpatID,
                         large_num ID,
                         t_vel K,
                         t_meter head,
                         double aquiferDepth,
                         double anisotropy,
                         double specificYield,
                         double specificStorage,
                         bool useEfolding,
                         bool confined,
                         large_num refID,
                         large_num maxRefinement,
                         bool isSteadyState,
                         bool isDensityVariable,
                         std::vector<t_dim> delnus,
                         std::vector<t_dim> nusInZones,
                         double effPorosity,
                         double maxTipSlope,
                         double maxToeSlope,
                         double minDepthFactor,
                         double slopeAdjFactor,
                         t_meter vdfLock,
                         int sourceZoneGHB,
                         int sourceZoneRecharge)
                    : NodeInterface(nodes, lat, lon, area, edgeLengthLeftRight, edgeLengthFrontBack, SpatID, ID, K,
                                    head, aquiferDepth, anisotropy, specificYield, specificStorage,
                                    useEfolding, confined, refID, maxRefinement, isSteadyState, isDensityVariable, delnus,
                                    nusInZones, effPorosity, maxTipSlope, maxToeSlope, minDepthFactor, slopeAdjFactor,
                                    vdfLock, sourceZoneGHB, sourceZoneRecharge) {}
        private:
            // implementation
            friend class NodeInterface;

            //Learning weight
            t_dim weight = 0.1;

            virtual t_meter
            __calcInitialHead(t_meter initialParam) {
                t_meter elevation = get<t_meter, TopElevation>();
                if (elevation >= initialParam) {
                    return elevation - initialParam; // todo check whether this is correct
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
                ar & initial_head;
                ar & simpleDistance;
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
                           t_meter edgeLengthFrontBack)
                    : NodeInterface(
                    nodes,
                    0,
                    0,
                    area,
                    edgeLengthLeftRight,
                    edgeLengthFrontBack,
                    ID,
                    ID,
                    0.3 * (si::meter / day), 1 * si::meter, 100, 10, 0.15,
                    0.000015, false, true, 0, 1, true, false, {0.0, 0.1}, {0.0, 0.1},
                    0.2, 0.2, 0.2, 0.1, 0.1, 0.001 * si::meter, 0, 0) {}

        private:
            friend class NodeInterface;

            virtual void __setHeadChange(t_meter change) {
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
                ar & initial_head;
                ar & simpleDistance;
                ar & fields;
                ar & cached;
                ar & eq_flow;
            }

        };

    }
}

#endif //NODE_HPP
