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
 */
        enum NeighbourPosition {
            TOP = 1,
            DOWN,               // 2
            RIGHT,              // 3
            LEFT,               // 4
            FRONT,              // 5
            BACK,               // 6
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
            t_dim efoldingFromData(t_meter verticalSize) {
                auto folding = get<t_meter, EFolding>();
                if (folding == 0.0 * si::meter)
                    return 1 * si::si_dimensionless;
                //Alter if a different size should be used and not full vertical size
                //depth = (depth / (2 * si::si_dimensionless));
                t_dim out = exp(-verticalSize / folding);
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
                t_vol_t storageFlow = calculateStorageFlow();
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

            /**
             * @brief Get the node identifier of this node
             * @return large_num
             */
            large_num getID() { return get<large_num, ID>(); }

            /**
             * @brief Get the spatial identifier of this node
             * @return int (for quick build of neighbouring)
             */
            int getSpatID() {return (int) get<large_num, SpatID>();}

            /**
             * @brief Get the longitude of the center of this node
             * @return double
             */
            double getLat() {return get<double, Lat>();}


            /**
             * @brief Get the latitude of the center of this node
             * @return double
             */
            double getLon() {return get<double, Lon>();}

            /**
             * @brief Get the surface area of this node
             * @return square meter
             */
            t_s_meter getArea(){return get<t_s_meter, Area>();}

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
                    t_vel out = get<t_vel, K>() * e_fold;
                    if (out < 1e-20 * si::meter / day) {
                        out = 1e-20 * si::meter / day;
                    }
                    return out;
                } else {
                    return get<t_vel, K>();
                }
            }

            /**
             * @brief Get hydraulic vertical conductivity scaled by anisotropy (and below first layer by e-folding)
             * @return meter per time
             */
            t_vel getK_vertical() noexcept { return (getK() / get<t_dim, Anisotropy>()); }

            /**
             * Calculate the equilibrium lateral flows
             * @return cubic meters per time
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
             * @return cubic meters per time
             */
            t_vol_t getLateralFlows() { return calcLateralFlows<Head>(false); }

            /**
             * Get the current lateral out flows
             * @return cubic meters per time
             */
            t_vol_t getLateralOutFlows() { return calcLateralFlows<Head>(true); }

            /**
             * @brief Get the effective porosity
             * @return dimensionless
             */
            t_dim getEffectivePorosity(){return get<t_dim, EffectivePorosity>();}

            /**
             * @brief Get the length of this node from left end to right end
             * @return meter
             */
            t_meter getEdgeLengthLeftRight(){return get<t_meter, EdgeLengthLeftRight>();}

            /**
             * @brief Get the length of this node from front end to back end
             * @return meter
             */
            t_meter getEdgeLengthFrontBack(){return get<t_meter, EdgeLengthFrontBack>();}

            /**
             * @brief Get the elevation of this node
             * @return
             */
            t_meter getElevation(){return get<t_meter, Elevation>();}

            /**
             * @brief Get the vertical size of this node
             * @return
             */
            t_meter getVerticalSize(){return get< t_meter, VerticalSize >(); }

            /**
             * @brief Get the height of the bottom of this node
             * @return meter
             */
            t_meter getBottom(){return get<t_meter, Elevation>() - get<t_meter, VerticalSize>();}

            /**
             * @brief Get the groundwater head
             * @return meter
             */
            t_meter getHead(){ return get<t_meter, Head>(); }

            /**
             * @brief Get the groundwater head of previous time step
             * @return meter
             */
            t_meter getHead_TZero() { return get<t_meter, Head_TZero>(); }

            /**
             * @brief Get the minimum of density interface closeness (are snapped together if closer)
             * @return meter
             */
            t_meter getVDFLock(){ return get<t_meter, VDFLock>(); }

            /**
             * @brief Get the layer
             * @return integer
             */
            int getLayer(){ return get<int, Layer>(); }

            /**
             * @brief Get the step size
             * @return dimensionless
             */
            t_dim getStepSize(){ return get<t_dim, StepSize>(); }

            /**
             * @brief Get all outflow since simulation start
             * @return cubic meter per time
             */
            t_vol_t getOUT() noexcept { return get<t_vol_t, OUT>(); }

            /**
             * @brief Get all inflow since simulation start
             * @return cubic meter per time
             */
            t_vol_t getIN() noexcept { return get<t_vol_t, IN>(); }

            /**
             * @brief Get all current IN flow
             * @return cubic meters per time
             */
            t_vol_t getCurrentIN() noexcept { return getFlow([](double a) -> bool { return a > 0; }); }

            /**
             * @brief Get all current OUT flow
             * @return cubic meters per time
             */
            t_vol_t getCurrentOUT() noexcept { return -getFlow([](double a) -> bool { return a < 0; }); }

            t_vol_t getVDF_IN() { return get<t_vol_t, VDF_IN>(); }

            t_vol_t getVDF_OUT() { return get<t_vol_t, VDF_OUT>(); }

            /**
             * @brief Get the zone change out of node
             * @return cubic meters per time
             */
            t_vol_t getVDF_ZCHG_OUT() { return get<t_vol_t, VDF_ZCHG_OUT>(); }

            /**
             * @brief Get the zone change into node
             * @return cubic meters per time
             */
            t_vol_t getVDF_ZCHG_IN() { return get<t_vol_t, VDF_ZCHG_IN>(); }

            t_vol_t getVDF_INST_OUT() { return get<t_vol_t, VDF_INST_OUT>(); }

            t_vol_t getVDF_INST_IN() { return get<t_vol_t, VDF_INST_IN>(); }

            t_vol_t getVDF_TTT_OUT() { return get<t_vol_t, VDF_TTT_OUT>(); }

            t_vol_t getVDF_TTT_IN() { return get<t_vol_t, VDF_TTT_IN>(); }

            /**
             * Calculate the lateral groundwater flow to (-)/from (+) a neighbouring node
             * @return unordered map of neighbor positions to flow values as double
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
                        lateralFlow = -conductance * (getHead() - getAt<t_meter, Head>(neigNodeID));
                        mapNeighboursToFlows[neigPos] = lateralFlow.value() ;
                    }


                }
                return mapNeighboursToFlows;
            }

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

            void setNewExternalFlowElevation(FlowType type, double head) {
                if (hasTypeOfExternalFlow(type)) {
                    auto conductance = externalFlows.at(type).getConductance().value();
                    auto bottom = externalFlows.at(type).getBottomElev();
                    removeExternalFlow(type);
                    addExternalFlow(type,
                                    head * Model::si::meter,
                                    conductance,
                                    bottom);
                }
            }

            /**
            * @brief Get bottom of an external flow by its FlowType
            * @param type The flow type
            * @return Bottom
            */
            t_meter getExternalFlowBottom(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type).getBottomElev();
                } else {
                    return 0 * si::meter;
                }
            }

            /**
             * @brief Get the number of neighbours
             * @return integer
             */
            int getNumOfNeighbours() { return (int) neighbours.size(); }

            /**
             * @brief Get the number of horizontal neighbours
             * @return integer
             */
            int getNumOfHorizontalNeighbours() { return (int) horizontal_neighbours.size(); }

            /**
             * @brief Get the list of neighbours (their position and identifier)
             * @return unordered map of neighbour position (e.g., left) and identifier
             */
            std::unordered_map<NeighbourPosition, large_num> getListOfNeighbours(){
                return neighbours;
            }

            /**
             * @brief Get the list of neighbours (their position and identifier)
             * @return unordered map of neighbour position (e.g., left) and identifier
             */
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
             * @brief The number of external flows
             * @return int
             */
            int getNumOfExternalFlows(){ return numOfExternalFlows; }

/*****************************************************************
Set Properties
******************************************************************/

            /**
             * @brief Set elevation on top layer and propagate to lower layers
             * @param elevation The top elevation (e.g. from DEM)
             */
            void setElevation_allLayers(const t_meter& elevation) {
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
             * @brief Set aquifer depth from data
             * @param aquifer_depth aquifer depth value in meters
             */
            void setVerticalSize(const t_meter& aquifer_depth) { set < t_meter, VerticalSize > (aquifer_depth); };

            /**
             * @brief Set e-folding factor from data on all layers
             * @param efold efolding value
             */
            void setEfold(double efold) {
                set < t_meter, EFolding > (efold * si::meter);
                applyToAllLayers([efold](NodeInterface *nodeInterface) {
                    try { nodeInterface->setEfold(efold); }
                    catch (...) {}
                });
            };

            /**
             * @brief Set the groundwater head
             * @param head groundwater head
             */
            void setHead(const t_meter& head) noexcept {
                NANChecker(head.value(), "Set Head");
                set<t_meter, Head>(head);
            }

            /**
             * @brief Set the same groundwater head in all aquifer layers
             * @param head groundwater head
             */
            void setHead_allLayers(t_meter head) noexcept {
                setHead(head);
                applyToAllLayers([&head](NodeInterface *nodeInterface) {
                    try {
                        nodeInterface->setHead(head);
                    }
                    catch (...) {}

                });
            }

            /**
             * @brief Set the equilibrium groundwater head using water table depth
             * @param wtd water table depth
             */
            void setEqHead_allLayers(const t_meter& wtd) {
                auto eqhead = get<t_meter, Elevation>() - wtd;
                set<t_meter, EQHead>(eqhead);
                applyToAllLayers([eqhead](NodeInterface *nodeInterface) {
                    try {
                        nodeInterface->set<t_meter, EQHead>(eqhead);
                    }
                    catch (...) {}
                });
            }

            /**
             * @brief Write equilibrium groundwater head onto "regular" groundwater head
             * @note Should be called after setting equilibrium groundwater head
             */
            void setHeadToEqHead_allLayers() {
                auto eqhead = get<t_meter, EQHead>();
                setHead(eqhead);
                applyToAllLayers([eqhead](NodeInterface *nodeInterface) {
                    try {
                        nodeInterface->setHead(eqhead);
                    }
                    catch (...) {}
                });
            }

            /**
             * @brief Write "regular" groundwater head onto equilibrium groundwater head
             * @note Should be called after setting "regular" groundwater head
             */
            void setEqHeadToHead_allLayers() {
                auto head = getHead();
                set<t_meter, EQHead>(head);
                applyToAllLayers([head](NodeInterface *nodeInterface) {
                    try {
                        nodeInterface->set<t_meter, EQHead>(head);
                    }
                    catch (...) {}
                });
            }

            /**
             * @brief Update the head change of previous time step
             * @note Should only be called at end of time step
             */
            void updateHeadChange_TZero() noexcept {
                set < t_meter, HeadChange_TZero > (getHead() - getHead_TZero());
            }

            /**
             * @brief Update the groundwater head of previous time step using the current groundwater head
             */
            void updateHead_TZero() noexcept { setHead_TZero(getHead()); }

            /**
             * @brief Set the groundwater head of previous time step using the current groundwater head
             * @param head groundwater head
             */
            void setHead_TZero(const t_meter& head) noexcept { set<t_meter, Head_TZero>(head); }

            /**
             * @brief Set the groundwater head of previous time step in all layers using the current groundwater head
             * @param head groundwater head
             */
            void setHead_TZero_allLayers(const t_meter& head) noexcept {
                setHead_TZero(head);
                applyToAllLayers([&head](NodeInterface *nodeInterface) {
                    nodeInterface->setHead_TZero(head);
                });
            }

            /**
             * @brief Set the density zone into which water from the general head boundary flows
             * @param sourceZoneGHB density zone into which water from the general head boundary flows
             */
            void setSourceZoneGHB(int sourceZoneGHB) { set<int, SourceZoneGHB>(sourceZoneGHB); }

            /**
             * @brief Modify hydraulic conductivity in all layers
             * @param conduct conductivity (if e-folding enabled scaled on below first layer)
             */
            void setK_allLayers(t_vel conduct) {
                setK(conduct);
                applyToAllLayers([&conduct](NodeInterface *nodeInterface) {
                    nodeInterface->setK(conduct);
                });
            }

            /**
             * @brief Set hydraulic conductivity
             * @param conduct hydraulic conductivity
             */
            void setK(const t_vel& conduct) { set < t_vel, K > (conduct); }

            /**
             * @brief Update whether time step is steady state or not
             * @param isSteadyState boolean defining if time step is steady state or not
             */
            void updateIsSteadyState(bool isSteadyState) { set<bool, IsSteadyState>(isSteadyState); }

            /**
             * @brief Update whether time step solves variable density flow or not
             * @param isSteadyState boolean defining if time step solves variable density flow or not
             */
            void updateIsDensityVariable(bool isDensityVariable) { set<bool, IsDensityVariable>(isDensityVariable); }

            /**
             * @brief Update step size of groundwater flow
             * @param stepSize time passing in one time step
             */
            void updateStepSize(double stepSize) { set<t_dim, StepSize>(stepSize * si::si_dimensionless); }

            /**
             * @brief Adapt step size of variable density evolution
             */
            void updateZetaStepSize(int vdf_steps_per_head_step) {
                set<t_dim, ZetaStepSize>( getStepSize() / vdf_steps_per_head_step );
            }

            /**
             * @brief Tell cell to save its flow budget
             */
            void saveMassBalance() noexcept {
                fields.addTo<t_vol_t, OUT>(getCurrentOUT());
                fields.addTo<t_vol_t, IN>(getCurrentIN());
            }

            void saveVDFMassBalance() noexcept {
                fields.addTo<t_vol_t, VDF_ZCHG_OUT>(get< t_vol_t, VDF_ZCHG_CUR_OUT>());
                fields.addTo<t_vol_t, VDF_INST_OUT>(calculateInstantaneousMixing(false));
                fields.addTo<t_vol_t, VDF_TTT_OUT>(calculateTipToeTrackingZoneChange(false));
                set<t_vol_t, VDF_OUT>(getVDF_ZCHG_OUT() + getVDF_INST_OUT() + getVDF_TTT_OUT());

                fields.addTo<t_vol_t, VDF_ZCHG_IN>(get< t_vol_t, VDF_ZCHG_CUR_IN>());
                fields.addTo<t_vol_t, VDF_INST_IN>(calculateInstantaneousMixing(true));
                fields.addTo<t_vol_t, VDF_TTT_IN>(calculateTipToeTrackingZoneChange(true));
                set<t_vol_t, VDF_IN>(getVDF_ZCHG_IN() + getVDF_INST_IN() + getVDF_TTT_IN());
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


            void resetVDFBudget() {
                set < t_vol_t, VDF_OUT > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_OUT > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_ZCHG_OUT > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_ZCHG_IN > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_INST_OUT > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_INST_IN > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_TTT_OUT > (0 * si::cubic_meter / day);
                set < t_vol_t, VDF_TTT_IN > (0 * si::cubic_meter / day);
            }


            /**
             * @brief Scales river conduct by 50%
             * @warning Should only be used in spin up phase
             */
            void scaleRiverConduct() {
                eq_flow = eq_flow * 1.5;
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
                    for(auto const& imap: externalFlows) {
                        LOG(debug) << " " << imap.first;
                    }
                    throw "Number of external flows don't match";
                }
            }


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
             * @brief Get neighbour positions at left, right, front and back
             * @return vector of neighbour positions
             */
            static std::vector<NeighbourPosition>
            getNeigPos_LRFB(){
                return {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                        NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
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
                return lateral_flow * get<t_dim, StepSize>();
            }

            /**
             * @brief Get flow budget based on head change
             * @return cubic meters per time
             * Note: Water entering storage is treated as an outflow (-), that is a loss of water from the flow system
             * while water released from storage is treated as inflow (+), that is a source of water to the flow system
             */
            t_vol_t calculateStorageFlow() noexcept {
                return -getStorageCapacity() * get<t_meter, HeadChange_TZero>() / day; // FIXME: rename "day" to "time"?
            }

            /**
             * @brief Get flow budget of a specific external flows
             * @param flow An external flow
             * @return cubic meters per time
             * Note: Water entering storage is treated as an outflow (-), that is a loss of water from the flow system
             * while water released from storage is treated as inflow (+), that is a source of water to the flow system
             */
            t_vol_t calculateExternalFlowVolume(const ExternalFlow &flow) {
                t_vol_t ex;
                auto eq_head = get<t_meter, EQHead>();
                auto head = get<t_meter, Head>();
                t_vol_t eqFlow = getEqFlow();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                if (hasTypeOfExternalFlow(RECHARGE)) {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                }
                if (is(flow.getType()).in(RECHARGE, NET_ABSTRACTION)) {
                    ex = recharge;
                } else if (is(flow.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                    if (flow.isFlowHeadDependent(head)) {
                        ex = (flow.getP(eq_head, head, recharge, eqFlow) * head +
                              flow.getQ(eq_head, head, recharge, eqFlow));
                    } else { // flow is not head dependent when the head is below the bottom of the simulated cell
                        ex = (flow.getP(eq_head, head, recharge, eqFlow) * flow.getBottomElev() +
                              flow.getQ(eq_head, head, recharge, eqFlow));
                    }
                } else {  // GENERAL_HEAD_BOUNDARY (Question: what about FLOODPLAIN_DRAIN, EVAPOTRANSPIRATION, FAST_SURFACE_RUNOFF)
                    ex = (flow.getP(eq_head, head, recharge, eqFlow) * head + // = - ghb_conductance * gw_head
                          flow.getQ(eq_head, head, recharge, eqFlow)); // = ghb_conductance * ghb_elevation (often = 0)
                }
                return ex * get<t_dim, StepSize>();
            }

            /**
             * @brief Calculate dewatered flow
             * @return cubic meters per time
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
                        out += conductance_below * (head_n - elev);
                    }
                }

                if (neighbours.find(TOP) != neighbours.end()) {
                    large_num nodeIDTop = neighbours.find(TOP)->first;
                    auto elev = getAt<t_meter, Elevation>(nodeIDTop);
                    auto head_n = getAt<t_meter, Head>(nodeIDTop);
                    //Check if a dewatered condition is present
                    if (get<t_meter, Head>() < get<t_meter, Elevation>() and head_n > elev) {
                        t_s_meter_t conductance_above =
                                mechanics.calculateVerticalConductance(createDataTuple(nodeIDTop));
                        out += conductance_above * (get<t_meter, Elevation>() - get<t_meter, Head>());
                    }
                }
                NANChecker(out.value(), "Dewatered flow");
                return out;
            }

            /**
             * @brief Calculates zone change budget of a node at zetaID (for variable density flow budget)
             * @return cubic meters per time
             * @note in SWI2 code: SSWI2_ZCHG
             */
            t_vol_t calculateZoneChange(large_num zetaID) {
                t_vol_t out = 0.0 * si::cubic_meter / day;

                t_meter zeta = getZeta(zetaID);
                t_meter zetaBelow = getZeta(zetaID + 1);
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

                t_meter zetaOld = getZeta_TZero(zetaID);
                t_meter zetaBelowOld = getZeta_TZero(zetaID + 1);
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
                        ((zeta - zetaBelow) - (zetaOld - zetaBelowOld)) / day;
                return out;
            }

            /**
             * @brief Instantaneous mixing of water, as described in SWI2 documentation under "Vertical Leakage Between
             * Aquifers": (2) when freshwater leaks up into an aquifer containing only saline water, that freshwater is
             *                added as saline water
             *            (4) when saline water leaks down into an aquifer containing only freshwater, that saline water
             *                is added as freshwater
             * @param zetaID
             * @note in SWI2 code: SSWI2_IMIX, comments referring to lines (at each "if"), refer to lines in gwf2swi27.f
             */
            t_c_meter calcInstMix(large_num zetaID) {
                t_c_meter out = 0.0 * si::cubic_meter;
                // skip nodes that do not have a down neighbour
                if (neighbours.find(DOWN) == neighbours.end()) { return out; } // line 4142

                // skip if dimensionless density at bottom of this node is below or equal to
                // dimension less density at the top of down node
                large_num nodeIDDown = neighbours[DOWN];
                if (getNusBot() <= at(nodeIDDown)->getNusTop()){ return out; } // line 4150

                // skip if head in this or neighbour node is below node bottom
                if (getHead() < getBottom() or at(nodeIDDown)->getHead() < at(nodeIDDown)->getBottom()){ return out; } // lines 4140 and 4145

                // skip if zetaID is not involved
                auto nusInZones = getNusInZones();
                if (nusInZones[zetaID] != getNusBot()) { return out; } // line 4156
                if (nusInZones[zetaID] != at(nodeIDDown)->getNusTop()) { return out; } // line 4162

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

            /**
             * @brief Calculate mass change from instantaneous mixing in variable density routine
             * @param in boolean defining if flow in or out should be returned
             * @return cubic meters per time
             */
            t_vol_t calculateInstantaneousMixing(bool in) noexcept {
                t_c_meter iMix_in;
                t_c_meter iMix_out;
                for (large_num zetaID = 0; zetaID < zetasSize() - 1; ++zetaID) {
                    t_c_meter iMix = calcInstMix(zetaID);
                    if (iMix.value() > 0) { iMix_in += iMix; } else { iMix_out += iMix; }
                }
                if (in) { return iMix_in / day; } else { return iMix_out / day; }
            }

            /**
             * @brief Save volumetric density zone change between last and new time step
             */
            void saveCurrentZoneChange() noexcept {
                set<t_vol_t, VDF_ZCHG_CUR_IN>(getZoneChange(true));
                set<t_vol_t, VDF_ZCHG_CUR_OUT>(getZoneChange(false));
            }

            /**
             * @brief Calculate mass change due to density zone change
             * @param in boolean defining if flow in or out should be returned
             * @return cubic meters per time
             */
            t_vol_t getZoneChange(bool in) noexcept {
                t_vol_t zoneChange_in = 0 * si::cubic_meters / day;
                t_vol_t zoneChange_out = 0 * si::cubic_meters / day;
                for (large_num zetaID = 0; zetaID < zetasSize() - 1; ++zetaID) {
                    t_vol_t zoneChange = calculateZoneChange(zetaID);
                    if (zoneChange.value() > 0) { zoneChange_in += zoneChange; } else { zoneChange_out += zoneChange; }
                }
                if (in) { return zoneChange_in; } else { return zoneChange_out; }
            }

            /**
             * @brief Calculate mass change due to tip and toe tracking (horizontal move of interface into neighbor node)
             * @param in boolean defining if flow in or out should be returned
             * @return cubic meters per time
             */
            t_vol_t calculateTipToeTrackingZoneChange(bool in) {
                t_vol_t result = 0 * si::cubic_meters / day;
                t_vol_t tttOut = getZoneChange(false) - get< t_vol_t, VDF_ZCHG_CUR_OUT>();
                t_vol_t tttIn = getZoneChange(true) - get< t_vol_t, VDF_ZCHG_CUR_IN>();
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

            class NodeNotFoundException : public std::exception {
                virtual const char *what() const throw() { return "Node does not exist"; }
            };

            /**
             * @brief Add an external flow to the cell
             * @param type The flow type
             * @param flowHead The flow head
             * @param cond The conductance
             * @param bottom The bottom of the flow (e.g river bottom)
             * @return Number assigned by cell to flow
             */
            int addExternalFlow(FlowType type, t_meter flowHead, double cond, const t_meter& bottom) {
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
                    for(auto const& imap: externalFlows) {
                        LOG(debug) << " " << imap.first;
                    }
                    throw "Number of external flows don't match";
                }
                return numOfExternalFlows;
            }


            void initializeZetas() {
                std::vector<t_meter> zetas;
                // add top zeta (at node elevation or gw head or node bottom)
                zetas.push_back(applyFrontZetaLimits());
                // add bottom zeta at node bottom
                zetas.push_back(getBottom());
                setZetas_direct(zetas);
            }

            /**
             * @brief Add a zeta surface to the cell (bounded by elevation at top and by cell bottom at bottom).
             * @param zetaID
             * @param zeta the zeta surface height in meters
             */
            void addZeta(large_num zetaID, t_meter zeta){
                NANChecker(zeta.value(), "zeta (in addZeta)");
                if (zetaID == 0) {
                    LOG(userinfo) << "zetaID should not be 0 when adding zetas";
                    throw "zetaID should not be 0 when adding zetas"; }
                auto zetas = getZetas();

                if (zeta > zetas.front() - getVDFLock()) {
                    zeta = zetas.front();
                    // limit zeta height to bottom zeta (in SWI2: lines 660-680)
                } else if (zeta < zetas.back() + getVDFLock()) {
                    zeta = zetas.back();
                }
                zetas.insert(zetas.begin() + zetaID , zeta);
                checkZetas();

                if (zetas.back() != getBottom()) {
                    LOG(userinfo) << "At nodeID " << getID() << ": zeta at back is below bottom of node!";
                    throw "Zeta at back must be at the bottom of the node!";
                }

                setZetas_direct(zetas);
                setZetas_TZero_direct(zetas); // give zetas before the iteration an initial value
            }

            void checkZetas() {
                auto zetas = getZetas();
                // check if zetas vector is sorted
                for (int id = 0; id < zetas.size()-1; ++id) {
                    if (zetas[id] < zetas[id+1]) {
                        LOG(userinfo) << "Zetas in wrong order at node(" << getID() << "): zeta(" << id << ") = "
                                      << zetas[id].value() << ", zeta(" << id+1 << ") = "
                                      << zetas[id+1].value();
                        throw "Vector of zetas needs to be sorted!";
                    }
                }

                if (zetas.front() != applyFrontZetaLimits()) {
                    LOG(userinfo) << "zeta.front() is out of bounds";
                }

                if (zetas.back() != getBottom()) {
                    LOG(userinfo) << "zeta.back() is out of bounds";
                }
            }

            /**
             * @brief Update zetas after one or multiple inner iteration
             * @param zetaID zeta surface id in this node
             * @param zeta zeta surface height in meter
             */
            void setZeta(large_num zetaID, const t_meter& zeta, bool setZetaTZero = false) {
                NANChecker(zeta.value(), "zeta (in setZeta)");
                auto zetas = getZetas();
                if (0 < zetaID < zetasSize()-1) {
                    zetas[zetaID] = applyInnerZetaLimits(zetaID, zeta);
                    setZetas_direct(zetas);
                    if (setZetaTZero) { setZetas_TZero_direct(zetas); }
                } else {
                    LOG(userinfo) << "zetaID too large in setZeta";
                    throw "zetaID too large in setZeta";
                }
            }

            void setZeta_direct(large_num zetaID, const t_meter& zeta, bool setZetaTZero = false) {
                NANChecker(zeta.value(), "zeta (in setZeta)");
                auto zetas = getZetas();
                if (0 < zetaID < zetasSize()-1) {
                    zetas[zetaID] = zeta;
                    setZetas_direct(zetas);
                    if (setZetaTZero) { setZetas_TZero_direct(zetas); }
                } else {
                    LOG(userinfo) << "zetaID too large in setZeta";
                    throw "zetaID too large in setZeta";
                }
            }

            t_meter applyFrontZetaLimits(){
                // the top zeta (zetaID = 0) is at GW head, limited by the top and bottom elevation of the node
                // info: the bottom zeta (set to the bottom of the node in "initializeZetas") is never changed
                if (getHead() > getBottom()) {
                    // if groundwater head is ABOVE the top of the node
                    if (getHead() > getElevation()) {
                        // top zeta is max at node elevation
                        return getElevation();
                        // if groundwater head is BELOW the top of the node
                    } else {
                        // top zeta is at gw head
                        return getHead();
                    }
                } else { // if head is below node bottom
                    return getBottom();
                }
            }

            /**
             * @brief Apply limits for inner zeta interfaces (e.g., needs to be below upper and above lower zeta)
             * @param zetaID zeta surface id in this node
             * @param zeta zeta surface height in meter
             * @return the limited zeta surface height
             */
            t_meter applyInnerZetaLimits(large_num zetaID, const t_meter& zeta) {
                // all inner zetas (except the last one) are limited by
                //  - the zeta height of the top and bottom zetas (in SWI2: lines 660-680)
                //  - the height of the GHB of this node (can not be higher than GHB)
                //  - the zeta height of the zetas above and below
                if (getEffectivePorosity().value() == 0) { // no saline water is in node if the effective porosity is 0
                    return getBottom();
                }
                auto zeta_tmp = zeta;
                auto zetas = getZetas();

                // upper limit of zetas at the coast: the elevation of the ocean (= GHB elevation)
                if (hasGHB()) {
                    auto ghbElevation = getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY);
                    if (zetas.back().value() < ghbElevation and ghbElevation < zeta_tmp.value()) {
                        zeta_tmp = ghbElevation * si::meter;
                    }
                } else {
                    // upper limit of zetas inland (i.e., not at the coast): the highest neighbouring zeta
                    // -> 1) find the highest neighbouring zeta
                    t_meter max_neig_zeta = zetas.back(); // lower limit of the new height: the lowest zeta (= node bottom)
                    t_meter neig_zeta;
                    for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                        neig_zeta = at(neigNodeID)->getZeta(zetaID);
                        if (neig_zeta > max_neig_zeta) {
                            max_neig_zeta = neig_zeta;
                        }
                    }
                    // -> 2) limit zeta
                    if (zeta_tmp > max_neig_zeta) { zeta_tmp = max_neig_zeta; }
                }

                // lower limit of the new height: the lowest zeta (= node bottom)
                if (zeta_tmp < zetas.back() + getVDFLock()) {
                    zeta_tmp = zetas.back();
                }
                // upper limit of the new height: the highest zeta (= node top or groundwater head)
                if (zeta_tmp > zetas.front() - getVDFLock()) {
                    zeta_tmp = zetas.front();
                }
                return zeta_tmp;
            }

            /**
             * @brief Set density interfaces directly (without applying limits)
             * @param zetas vector of density interfaces
             */
            void setZetas_direct(const std::vector<t_meter>& zetas) { set<std::vector<t_meter>, Zetas>(zetas); }

            void limitZetas(bool limitZetasTZero = false) {
                auto zetas = getZetas();
                for (large_num zetaID = 1; zetaID < zetasSize() - 1; ++zetaID) {
                    zetas[zetaID] = applyInnerZetaLimits(zetaID, zetas[zetaID]);
                }
                setZetas_direct(zetas);
                if (limitZetasTZero) { setZetas_TZero_direct(zetas); }
            }

            /**
             * @brief Adjust zeta surfaces if they have crossed or their height difference is below a minimum value
             * @param zetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ZETACROSS
             *  zetaID      | n_above | n
             *  1           | -       | -
             *  2           | 1       | 1,2,3
             *  3           | 1,2     | 2,3,4; 1,2,3,4
             *  4           | 1,2,3   | 3,4,5; 2,3,4,5; 1,2,3,4,5
             */
            void correctCrossingZetas(){
                t_meter zetaAverage;
                t_meter zetaSum;
                large_num zetaID_above;
                int counter;
                for (large_num zetaID = 1; zetaID < zetasSize() - 2; ++zetaID) {
                    // if zeta surface is very close to or lower than the zeta surface that SHOULD be below
                    if (getZeta(zetaID) - getZeta(zetaID + 1) < getVDFLock()) {
                        // make the zeta height of both surfaces their average
                        zetaAverage = 0.5 * (getZeta(zetaID) + getZeta(zetaID + 1));
                        setZeta(zetaID, zetaAverage); // set this zeta surface to average
                        setZeta(zetaID + 1, zetaAverage); // set zeta surface "below" to average
                        // if there are zeta surfaces above which could possibly be crossing
                        for (int n_above = 1; n_above < zetaID; ++n_above) {
                            zetaSum = 0 * si::meter;
                            // create a range from n_min (n-x) to n_max (n+1)
                            zetaID_above = zetaID - n_above;
                            // if surface zetaID_above crosses (or is very close to) this zeta surface
                            if (getZeta(zetaID_above) - getZeta(zetaID) < getVDFLock()) {
                                // sum heights from zetaID_above to zetaID
                                counter = 0;
                                for (int n = zetaID_above; n <= zetaID; ++n) {
                                    zetaSum += getZeta(n);
                                    ++counter;
                                }
                                if (counter != 0){
                                    // calculate the average height from surface zetaID_above to zetaID+1
                                    zetaAverage = zetaSum / (counter * si::si_dimensionless);
                                    // set all zeta surfaces from zetaRef to zetaID+1 to that average
                                    for (int n = zetaID_above; n <= zetaID; ++n) { setZeta(n, zetaAverage); }
                                }
                            }
                        }
                    }
                }
            }


            /**
             * @brief Get all zeta surface heights
             * @return vector<meter>
             */
            std::vector<t_meter> getZetas() { return get<std::vector<t_meter>, Zetas>();}

            large_num zetasSize() { return getZetas().size(); }

            /**
             * @brief get one zeta surface height
             * @param zetaID
             * @return meter
             */
            t_meter getZeta(large_num zetaID) {
                auto zetas = getZetas();
                if (zetaID < zetasSize()){
                    auto zeta = zetas[zetaID];
                    return zeta;
                } else {
                    LOG(userinfo) << "Not set at nodeID " + std::to_string(getID()) +
                                     ": Zetas[zetaID = " + std::to_string(zetaID) + "]";
                    throw "Not set at nodeID " + std::to_string(getID()) +
                          ": Zetas[zetaID = " + std::to_string(zetaID) + "]";
                }
            }

            std::vector<t_meter> getZetas_TZero() { return get<std::vector<t_meter>, Zetas_TZero>();}

            t_meter getZeta_TZero(large_num zetaID){
                if (zetaID < get<std::vector<t_meter>, Zetas_TZero>().size()){
                    return get<std::vector<t_meter>, Zetas_TZero>()[zetaID];
                } else {
                    LOG(userinfo) << "Not set at nodeID " + std::to_string(getID()) +
                                     ": Zetas_TZero[zetaID = " + std::to_string(zetaID) + "]";
                    throw "Not set at nodeID " + std::to_string(getID()) +
                          ": Zetas_TZero[zetaID = " + std::to_string(zetaID) + "]";
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
             * @param flow
             * @param lock
             */
            void updateUniqueFlow(double amount, FlowType flow = RECHARGE, bool lock = true) {
                if (lock and flow == RECHARGE) {
                    if (hasTypeOfExternalFlow(RIVER_MM)) {
                        //get current recharge and lock it before setting new recharge
                        //in arid regions recharge might be 0
                        t_vol_t recharge{0 * si::cubic_meter /day};
                        if(hasTypeOfExternalFlow(RECHARGE)){
                            recharge = getExternalFlowByName(RECHARGE).getRecharge();
                        }
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
                    LOG(userinfo) << "Number of external flows don't match";
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
                    t_meter bottom = getExternalFlowByName(type).getBottomElev();
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
                    RiverDepth = getExternalFlowByName(Model::FlowType::RIVER_MM).getFlowHead().value() - getExternalFlowByName(Model::FlowType::RIVER_MM).getBottomElev().value();
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
                    t_meter bottom = getExternalFlowByName(type).getBottomElev();
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
                    t_meter bottom = getExternalFlowByName(type).getBottomElev();
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
                    t_meter bottom{externalFlow.getBottomElev()};
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
                    t_meter bottom = getExternalFlowByName(type).getBottomElev() * amount; // TODO: if amount puts bottom upwards head might be under bottom
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
                auto eq_head = get<t_meter, EQHead>();
                auto head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                if (hasTypeOfExternalFlow(RECHARGE)) {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                }
                t_vol_t eqFlow = getEqFlow();
                t_vol_t ex = 0.0 * (si::cubic_meter / day);
                for (const auto &flow : externalFlows) {
                    ex += flow.second.getQ(eq_head, head, recharge, eqFlow);
                }
                return ex;
            }

            /**
             * @brief Get Q part (external sources) source part of RHS of current density zone
             * @return volume over time
             * @note The source zones of GHB and GWR are specified in config and/or in file, sink of GHB is the top zone
             */
            t_vol_t getQ(large_num zetaID) noexcept {
                large_num densityZone = zetaID; // the denity zone has the same ID as the zeta interface
                t_vol_t ex = 0.0 * (si::cubic_meter / day);
                t_vol_t Q = 0.0 * (si::cubic_meter / day);

                auto eq_head = get<t_meter, EQHead>();
                t_meter head = getHead(); // need the previous head
                t_vol_t recharge = 0 * si::cubic_meter / day;
                if (hasTypeOfExternalFlow(RECHARGE)) {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                }
                t_vol_t eqFlow = getEqFlow();
                for (const auto &[flowType, extFlow] : externalFlows) {
                    Q = extFlow.getQ(eq_head, head, recharge, eqFlow);
                    if (flowType == RECHARGE) {
                        // if this is the zone of recharge
                        if (densityZone == getSourceZoneRecharge()) {
                            //LOG(debug) << "adding gw recharge to sources below zeta";
                            ex += Q;
                        }
                    } else if (flowType == GENERAL_HEAD_BOUNDARY) {
                        // if GHB flow is into node (gw_head < ghb_elevation) and ...
                        if (getHead() < extFlow.getFlowHead()) {
                            // ... this is the zone of GHB inflow
                            if (densityZone == getSourceZoneGHB()) {
                                ex += Q;
                            }
                        // if GHB flow is out of node (gw_head > ghb_elevation) and ...
                        } else if (getHead() > extFlow.getFlowHead()) {
                            // ... interface above is not active (-> no zeta above active)
                            if (!isZetaActive(zetaID - 1)) {
                                ex -= Q;
                            }
                        }
                    }
                }
                return ex;
            }

            /**
             * @brief The effective porosity term for this node
             * @return square meter over time
             * @note In SSWI2: SWIHCOF (line 3744)
             */
            t_s_meter_t getEffectivePorosityTerm(){ // computed independent of steady or transient flow (see SWI2 doc "Using the SWI2 Package")
                t_s_meter_t out = (getEffectivePorosity() * get<t_s_meter, Area>()) / (day * get<t_dim, ZetaStepSize>());
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
                for (large_num zetaID = 1; zetaID < getZetas_TZero().size() - 1; zetaID++){
                    if (isZetaTZeroAtFront(zetaID) and getEffectivePorosity().value() > 0){
                        auto delnus = get<std::vector<t_dim>, Delnus>();
                        out += delnus[zetaID];
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
                for (large_num zetaID = 1; zetaID < getZetas_TZero().size() - 1; zetaID++){
                    if (isZetaTZeroAtBottom(zetaID) and getEffectivePorosity().value() > 0){
                        auto delnus = get<std::vector<t_dim>, Delnus>();
                        out -= delnus[zetaID];
                    }
                }
                return out;
            }

            /**
             * @brief Checking if zeta surface is at top zeta surface (and NOT at bottom)
             * @param zeta
             * @return bool
             */
            bool isZetaAtElevation(const large_num& zetaID){
                return (getZeta(zetaID) >= getElevation() - getVDFLock());
            }

            bool isZetaTZeroAtElevation(const large_num& zetaID){
                return ((getZeta_TZero(zetaID) >= getElevation() - getVDFLock()));
            }


            /**
             * @brief Checking if zeta surface is at top zeta surface (and NOT at bottom)
             * @param zeta
             * @return bool
             */
            bool isZetaAtFront(const large_num& zetaID){
                return (getZeta(zetaID) >= getZetas().front() - getVDFLock());
            }

            bool isZetaTZeroAtFront(const large_num& zetaID){
                return (getZeta_TZero(zetaID) >= getZetas_TZero().front() - getVDFLock());
            }

            /**
             * @brief Check if density interface is at node bottom
             * @param zetaID density interface identifier
             * @return boolean
             */
            bool isZetaAtBottom(const large_num& zetaID){
                return (getZeta(zetaID) <= getZetas().back() + getVDFLock());
            }

            /**
             * @brief Check if density interface was at node bottom in previous step
             * @param zetaID density interface identifier
             * @return boolean
             */
            bool isZetaTZeroAtBottom(const large_num& zetaID){
                return (getZeta_TZero(zetaID) <= getZetas_TZero().back() + getVDFLock());
            }

            /**
             * @brief Check if density interface is between first and last density interface
             * @param zetaID density interface identifier
             * @return bool
             */
            bool isZetaBetween(const large_num& zetaID){
                return (!isZetaAtFront(zetaID) and !isZetaAtBottom(zetaID));
            }

            /**
             * @brief Checking if density interface was between first and last density interface in previous time step
             * @param zetaID density interface identifier
             * @return bool
             */
            bool isZetaTZeroBetween(const large_num& zetaID){
                return (!isZetaTZeroAtFront(zetaID) and !isZetaTZeroAtBottom(zetaID));
            }

            /**
             * @brief Check if any zeta surface is active
             * @return bool
             */
            bool isAnyZetaActive(){
                for (large_num zetaID = 0; zetaID < zetasSize(); ++zetaID){
                    if (isZetaActive(zetaID)){ return true;}
                }
                return false;
            }

            /**
             * @brief Checking if zeta surface is active (between first and last interface and porosity above zero)
             * @param zetaID  density interface identifier
             * @return bool
             */
            bool isZetaActive(large_num zetaID){
                return (isZetaBetween(zetaID) and getEffectivePorosity().value() > 0);
            }

            /**
             * @brief Check if zeta surface is inactive (opposite of active)
             * @param zetaID  density interface identifier
             * @return bool
             */
            bool isZetaInactive(large_num zetaID) {
                return !isZetaActive(zetaID);
            }

            /**
             * @brief Check if zeta surface was active in previous time step
             * @param zetaID  density interface identifier
             * @return bool
             */
            bool isZetaTZeroActive(large_num zetaID) {
                return (isZetaTZeroBetween(zetaID) and getEffectivePorosity().value() > 0);
            }

            /**
             * @brief Check if zeta surface was inactive in previous time step
             * @param zetaID  density interface identifier
             * @return bool
             */
            bool isZetaTZeroInactive(large_num zetaID) {
                return !isZetaTZeroActive(zetaID);
            }

            /**
             * @brief Check if zeta surface qualifies as an output
             * @param zetaID  density interface identifier
             * @return bool
             */
            bool isZetaOutput(large_num zetaID){
                return (!isZetaAtBottom(zetaID) and getEffectivePorosity().value() > 0);
            }

            /**
             * @brief The source flow below a zeta surface (for the right hand side in the zeta equation)
             * @param zetaID zeta surface id in this node
             * @return cubic meters per time
             * @note like G in SWI2 doc, but without vertical leakage; in SWI2 code: lines 3523-3569
             * G = RHS (of flow, for constant density) - HCOF_(i,j,k,n)*h^(m)_(i,j,k) + (verticalLeakage_(i,j,k-1,n) - verticalLeakage_(i,j,k,n))
             */
            /*t_vol_t getSources(large_num zetaID){
                // We use the following sink/source concept:
                //  - the source zones of GHB and GWR are specified in config and/or in file, sink of GHB is the top zone
                //  - the source/sink zone of SWBs depend on interface heights

                // todo: compute BUFF with SSWI2_BDCH for constant head cells
                t_vol_t sources = 0.0 * (si::cubic_meter / day);
                // if the new groundwater head is above or equal to the node bottom
                if (getHead() >= getBottom()) { // lines 3532-3536
                    t_vol_t externalSources = - getQ(zetaID) - getP_aboveFlowBottom(zetaID) * getHead();
                    t_vol_t storageChange = getStorageCapacity() * (getHead() - getHead_TZero()) /
                                                    (day * get<t_dim, StepSize>());
                    sources = externalSources + storageChange * getZoneFraction(zetaID); // see line 3535-3536
                }
                return sources;
            }*/

            t_dim getZoneFraction(large_num zetaID) {
                if (zetaID >= getZetas_TZero().size() -1) {
                    LOG(userinfo) << "getZoneFraction(): zetaID too high";
                    throw "getZoneFraction(): zetaID too high";
                }
                // return: height of this zone / height of all zones
                return (getZeta_TZero(zetaID) - getZeta_TZero(zetaID + 1)) /
                       (getZetas_TZero().front() - getZetas_TZero().back());
            }

            /**
             * @brief The pseudo source for a zeta surface (for the right hand side in the zeta equation)
             * @param zetaID zeta surface id in this node
             * @return cubic meters per time
             * @note in SWI2 code lines 3574-3635, using SSWI2_SD and SSWI2_SR)
             */
            t_vol_t getPseudoSourceBelowZeta(large_num zetaID) {
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<std::vector<t_dim>, Delnus>();

                // pseudo source term calculation (in 2 parts)
                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    // on right side of equation -> need to use zeta interface heights of previous time step
                    std::vector<t_meter> zoneThicknesses = calculateZoneThicknesses(neigPos, neigNodeID, getZetas_TZero(),
                                                                                    at(neigNodeID)->getZetas_TZero());
                    // calculating zone conductances for pseudo source term calculation
                    std::vector<t_s_meter_t> zoneConductances = getZoneConductances(neigPos, neigNodeID,
                                                                                    zoneThicknesses);
                    t_s_meter_t zoneCondCum = getZoneConductanceCum(zetaID,zoneConductances);
                    if (isZetaTZeroActive(zetaID) and // if IPLPOS == 0 (line 3571)
                        at(neigNodeID)->isZetaTZeroActive(zetaID)) { // if neighbouring IPLPOS == 0 (e.g. line 2094)
                        //%% head part %%
                        t_vol_t head_part = zoneCondCum * (at(neigNodeID)->getHead() - getHead());
                        out -= head_part;
                        //LOG(debug) << "head_part (pseudo source): " << head_part.value() << std::endl;

                        t_s_meter_t zoneCondCumDelnus;
                        for (large_num zetaID_delnus = 0; zetaID_delnus < getZetas_TZero().size() - 1; zetaID_delnus++) {
                            //%% delnus part %%
                            if (zetaID_delnus < zetaID) {
                                zoneCondCumDelnus = zoneCondCum;
                            } else if (zetaID_delnus == zetaID) {
                                continue;
                            } else { //if (zetaID_delnus > zetaID) {
                                zoneCondCumDelnus = getZoneConductanceCum(zetaID_delnus, zoneConductances);
                            }
                            t_vol_t delnus_part = delnus[zetaID_delnus] * zoneCondCumDelnus *
                                                  (at(neigNodeID)->getZeta_TZero(zetaID_delnus) - getZeta_TZero(zetaID_delnus));
                            out -= delnus_part;
                            //LOG(debug) << "delnus_part (pseudo source) (zetaID_delnus = " << zetaID_delnus << "): " << delnus_part.value() << std::endl;
                        }
                    }
                }

                // add pseudo-source from GHB
                // (can be flow into or out of node, depends on the zeta height above the zone where GHB flows into)
                /*if (hasGHB() and getSourceZoneGHB() == zetaID) { // zoneID is the same as the zetaID
                    auto ghbElevation = externalFlows.at(GENERAL_HEAD_BOUNDARY).getFlowHead();
                    auto ghbConductance = externalFlows.at(GENERAL_HEAD_BOUNDARY).getConductance();
                    auto ghbZetaID = getSourceZoneGHB();

                    t_vol_t pseudoSourceGHB = delnus[zetaID] * ghbConductance * (ghbElevation - getZeta_TZero(zetaID));
                    out -= pseudoSourceGHB;
                }*/
                NANChecker(out.value(), "getPseudoSourceBelowZeta");
                return out;
            }

            /**
             * @brief Flow at tips and toes (on right hand side of zeta equation)
             * @param got neighbouring node
             * @param zetaID zeta surface id in this node
             * @return cubic meters per time
             */
            t_vol_t getFluxHorizontal(NeighbourPosition neigPos, large_num neigNodeID, large_num zetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<std::vector<t_dim>, Delnus>();

                // on right hand side of equation -> we need the zone conductances of previous time step
                std::vector<t_meter> zoneThicknesses = calculateZoneThicknesses(neigPos, neigNodeID, getZetas_TZero(),
                                                                                at(neigNodeID)->getZetas_TZero());
                std::vector<t_s_meter_t> zoneConductances = getZoneConductances(neigPos, neigNodeID,
                                                                                zoneThicknesses);

                t_s_meter_t zoneCondCum = getZoneConductanceCum(zetaID, zoneConductances);
                // %%head part %%
                t_vol_t head_part = zoneCondCum * (getAt<t_meter, Head>(neigNodeID) - get<t_meter, Head>());
                out += head_part;
                //LOG(debug) << "head_part (tip/toe): " << head_part.value() << std::endl;

                // %%delnus part %%
                t_s_meter_t zoneCondCumDelnus;
                for (large_num zetaID_delnus = 0; zetaID_delnus < getZetas_TZero().size() - 1; zetaID_delnus++) {
                    if (zetaID_delnus <= zetaID){
                        zoneCondCumDelnus = zoneCondCum;
                    } else {
                        zoneCondCumDelnus = getZoneConductanceCum(zetaID_delnus,zoneConductances);
                    }
                    t_vol_t delnus_part = delnus[zetaID_delnus] * zoneCondCumDelnus *
                                          (at(neigNodeID)->getZeta_TZero(zetaID_delnus) - getZeta_TZero(zetaID_delnus));
                    out += delnus_part;
                    //LOG(debug) << "delnus_part (tip/toe): " << delnus_part.value() << std::endl;
                }
                return out;
            }

            /**
             * @brief Specification of boundary condition at tips and toes
             * @param zetaID zeta surface id in this node
             * @return cubic meters per time
             * @note adapted from SWI2 code lines 3637-3703 (includes usage of SSWI2_QR and SSWI2_QC)
             */
            t_vol_t getTipToeFlow(large_num zetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);

                t_meter zeta_neig;
                // tip and toe flow calculation
                if (isZetaTZeroActive(zetaID)) { // if IPLPOS == 0 (line 3571)
                    for (auto const &[neigPos, neigNodeID]: neighbours) {
                        if (!at(neigNodeID)->isZetaTZeroActive(zetaID)) { // if IPLPOS == 0 (line 3571)
                            if (neigPos == NeighbourPosition::TOP) {
                                // vertical leakage to TOP neighbour
                                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                                if (getFluxTop() < 0 * (si::cubic_meter / day) and
                                    at(neigNodeID)->getNusBot() >= nusInZones[zetaID] and
                                    getNusBot() >= at(neigNodeID)->getNusBot()) { // IF ((qztop.LT.0).AND.(NUBOT(i,j,k-1).GE.NUS(iz)).AND.(NUBOT(j,i,k).GE.NUBOT(i,j,k-1))) THEN
                                    out += getFluxTop(); // in SWI2: qztop
                                }
                            } else if (neigPos == NeighbourPosition::DOWN) {
                                // vertical leakage to DOWN neighbour
                                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                                if (getFluxDown() < 0 * (si::cubic_meter / day) and
                                    at(neigNodeID)->getNusTop() < nusInZones[zetaID] and
                                    getNusTop() <= at(neigNodeID)->getNusTop()) { // IF ((qzbot.LT.0).AND.(NUTOP(i,j,k+1).LT.NUS(iz)).AND.(NUTOP(j,i,k).LE.NUTOP(i,j,k+1))) THEN
                                    continue;
                                } else{
                                    out += getFluxDown(); // in SWI2: qzbot
                                }
                            } else { // if neighbour is horizontal
                                out -= getFluxHorizontal(neigPos, neigNodeID, zetaID); // SSWI2_QR (left/right), SSWI2_QC (front/back)
                            }
                        }
                    }
                }
                return out;
            }

            /**
             * @brief The density zone conductance for variable density flow calculations (both flow and zeta equations)
             * @param neigPos position (LEFT, RIGHT...) of neighbouring node
             * @param neigNodeID nodeID of neigbouring node
             * @param zoneThicknesses vertical thicknesses of density zones
             * @return vector with entries in square meter per time
             */
            std::vector<t_s_meter_t> getZoneConductances(NeighbourPosition neigPos, large_num neigNodeID,
                                                         const std::vector<t_meter>& zoneThicknesses) {
                std::vector<t_s_meter_t> out;
                t_s_meter_t conductance;
                t_s_meter_t zoneConductance;

                // calculate the zone thickness sum
                t_meter sumOfZoneThicknesses = 0 * si::meter;
                for (const auto &zoneThickness : zoneThicknesses) {
                    sumOfZoneThicknesses += zoneThickness;
;                }

                // calculate the density zone conductances
                for (const auto &zoneThickness : zoneThicknesses) {
                    //LOG(debug) << "zoneThickness: " << zoneThickness.value();
                    zoneConductance = 0 * si::square_meter / day;
                    if (sumOfZoneThicknesses != (0 * si::meter)) { // adapted from SWI2 code line 1159
                        conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                        //LOG(debug) << "conductance:" << conductance.value();
                        zoneConductance = conductance * (zoneThickness / sumOfZoneThicknesses);
                        //LOG(debug) << "zoneConductance: " << zoneConductance.value();
                    }
                    out.push_back(zoneConductance);
                    NANChecker(zoneConductance.value(), "zoneConductance");
                }
                return out;
            }

            /**
             * @brief calculate thicknesses of density zones
             * @param neigPos position (LEFT, RIGHT...) of neighbouring node
             * @param neigNodeID nodeID of neighbouring node
             * @param zetas vector of density interface heights
             * @param zetas_neig vector of density interface heights in neighbouring node
             * @return
             */
            std::vector<t_meter> calculateZoneThicknesses(NeighbourPosition neigPos, large_num neigNodeID,
                                                          const std::vector<t_meter>& zetas,
                                                          const std::vector<t_meter>& zetas_neig) {
                std::vector<t_meter> zoneThicknesses;
                t_meter zoneThickness;
                t_meter deltaZeta;
                t_meter deltaZeta_neig;
                for (large_num zetaID = 0; zetaID < zetas.size() - 1; zetaID++) {
                    deltaZeta = zetas[zetaID] - zetas[zetaID + 1];
                    deltaZeta_neig = zetas_neig[zetaID] - zetas_neig[zetaID + 1];
                    //LOG(debug) << "nodeID: " << getID() << ", zetas[" << zetaID << "]: " << zetas[zetaID].value() <<
                    //           ", zetas[" << zetaID+1 << "]: " << zetas[zetaID+1].value();

                    //LOG(debug) << "neigNodeID: " << neigNodeID << ", zetas_neig[" << zetaID << "]: " << zetas_neig[zetaID].value() <<
                    //           ", zetas_neig[" << zetaID+1 << "]: " << zetas_neig[zetaID+1].value();

                    if (deltaZeta <= (0 * si::meter) or deltaZeta_neig <= (0 * si::meter)){ // adapted from SWI2 code line 1149
                        zoneThickness = 0 * si::meter;
                    } else {
                        zoneThickness = ((getLengthNeig(neigPos, neigNodeID) * deltaZeta) +
                                         (getNodeLength(neigPos) * deltaZeta_neig)) /
                                        (getLengthNeig(neigPos, neigNodeID) + getNodeLength(neigPos));
                    }
                    //LOG(debug) << "nodeID: " << getID() << ", zoneThickness: " << zoneThickness.value();
                    NANChecker(zoneThickness.value(), "zoneThicknessTZero");
                    zoneThicknesses.push_back(zoneThickness);
                }
                return zoneThicknesses;
            }

            /**
             * @brief Sum of density zone conductances below current local zeta id
             * @param zetaID identifier of density interface
             * @param zoneConductances density zone conductances
             * @return square meter per time
             */
            static t_s_meter_t getZoneConductanceCum(large_num zetaID, std::vector<t_s_meter_t> zoneConductances) {
                // calculate the sum of density zone conductances below a zeta surface n and add to vector out
                t_s_meter_t out = 0 * si::square_meter / day;

                for (int zoneID = zetaID; zoneID < zoneConductances.size() ; zoneID++) {
                    out += zoneConductances[zoneID];
                    //LOG(debug) << "zoneConductances[" << zoneID << "]: " << zoneConductances[zoneID].value();
                };
                //LOG(debug) << "zoneConductanceCum:" << out.value() << std::endl;
                NANChecker(out.value(), "getZoneConductanceCum");
                return out;
            }

            /**
             * @brief The pseudo-source for the flow equation, only used if variable density flow is active
             * @return cubic meters per time
             * @note This accounts for the effects of variable density flow (in SWI2 code: QREXTRA/QFEXTRA)
             */
            t_vol_t getPseudoSourceNode() {
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<std::vector<t_dim>, Delnus>();
                // pseudo-source calculation
                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    if (isAnyZetaActive() or at(neigNodeID)->isAnyZetaActive()) { // check if there are any active zeta surfaces
                        // in flow equation -> zeta interface heights are not updated between iterations
                        std::vector<t_meter> zoneThicknesses = calculateZoneThicknesses(neigPos, neigNodeID,
                                                                                        getZetas_TZero(),
                                                                                        at(neigNodeID)->getZetas_TZero());
                        std::vector<t_s_meter_t> zoneConductances = getZoneConductances(neigPos, neigNodeID,
                                                                                        zoneThicknesses);
                        for (large_num zetaID = 0; zetaID < zetasSize() - 1; zetaID++) {
                            t_s_meter_t zoneConductanceCum = getZoneConductanceCum(zetaID, zoneConductances);
                            //LOG(debug) << "zoneConductanceCum: " << zoneConductanceCum.value();
                            t_vol_t pseudoSource = delnus[zetaID] * zoneConductanceCum *
                                                   (at(neigNodeID)->getZeta(zetaID) - getZeta(zetaID));

                            out -= pseudoSource;
                        }
                    }
                }

                // not in SWI2: add pseudo-source from GHB
                /*if (hasGHB()) {
                    auto ghbElevation = externalFlows.at(GENERAL_HEAD_BOUNDARY).getFlowHead();
                    auto ghbCondcutance = externalFlows.at(GENERAL_HEAD_BOUNDARY).getConductance();
                    // get the zetaID at the top of the density zone into which GHB flows
                    auto ghbZetaID = getSourceZoneGHB(); // zetaID = id of zone below
                    t_vol_t pseudoSourceGHB = delnus[ghbZetaID] * ghbCondcutance * (ghbElevation - getZeta(ghbZetaID));
                    out -= pseudoSourceGHB;
                }*/
                return out;
            }

            /**
             * @brief The flux correction term in vertical direction
             * @return cubic meters per time
             * @note in SWI2 documentation: CV*BOUY; in code: QLEXTRA
             */
            t_vol_t getVerticalFluxCorrection(){
                auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                t_meter headdiff = 0 * si::meter;
                t_vol_t out = 0 * (si::cubic_meter / day);

                if (neighbours.find(TOP) != neighbours.end()){ //Current node has a top node
                    large_num topNodeID = neighbours[TOP];
                    // first part of the flux correction term
                    for (large_num zetaID = 0; zetaID < getZetas_TZero().size() - 1; zetaID++){
                        headdiff -= nusInZones[zetaID] *
                                    (at(topNodeID)->getZeta_TZero(zetaID + 1) - at(topNodeID)->getZeta_TZero(zetaID));
                        // Note: in SWI2 documentation is, BOUY is calculated by adding headdiff (would be out +=),
                        // MODFLOW code for headdiff is as implemented here (with out -=)
                    }
                    // second part of the flux correction term (vertical conductance *  BOUY)
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(topNodeID));
                    out = verticalConductance *
                          (headdiff +
                           0.5 * (at(topNodeID)->getZetas_TZero().back() - getZetas_TZero().front()) *
                           (at(topNodeID)->getNusBot() + getNusTop()));
                    // Note in SWI2 documentation, BOUY is calculated with a - between NUBOT and NUTOP,
                    // in MODFLOW code there is a + in the calculation of QLEXTRA
                    //LOG(debug) << "headdiff: " << headdiff.value() << std::endl;
                }
                return out;
            }

            /**
             * @brief Calculates the vertical flux correction
             * @return cubic meters per time
             */
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
             * @return cubic meters per time
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
             * @return cubic meters per time
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
             * @brief Updating top zeta height after the head solution is found
             * @note Top zeta surface height is set to the new groundwater head. In SWI2: SSWI2_UPZ1
             */
            void clipFrontZeta(){
                t_meter frontZeta = applyFrontZetaLimits();
                // update the first zeta surface
                auto zetas = getZetas();
                zetas[0] = frontZeta;
                setZetas_direct(zetas);
                // update all other zeta surfaces if they are ABOVE the new height of the first zeta surface
                bool limitZetasTZero{true};
                limitZetas(limitZetasTZero);
            }

            /**
             * @brief save the density interfaces from previous time step
             * @param zetas vector of density interface heights
             */
            void setZetas_TZero_direct(std::vector<t_meter>& zetas) { set<std::vector<t_meter>, Zetas_TZero>(zetas); };

            /**
             * @brief Vertical movement of zeta surfaces through top of this node. This function is required to move a
             * zeta surface with the height equal to the top of a lower node and bottom to an upper node. Without this
             * function, that zeta surface would be stuck there since the other only consider "active" surfaces (which
             * are between the top and bottom of a node)
             * @note in SWI2 code: SSWI2_VERTMOVE
             */
            void zetaMovementBetweenLayers() {
                t_vol_t fluxCorrectionTop; // in SWI2: qztop
                t_meter deltaZeta;
                // calculate flux through the top (checks if there is a neighbour at the top)
                fluxCorrectionTop = getFluxTop();
                //LOG(debug) << "fluxCorrectionTop: " << fluxCorrectionTop.value() << std::endl;
                if (fluxCorrectionTop.value() != 0) {
                    large_num topNodeID = neighbours[TOP];

                    if (getEffectivePorosity().value() <= 0 or at(topNodeID)->getEffectivePorosity().value() <= 0){
                        return;
                    }

                    // if head is above the bottom of the node, both in this node AND in the node above (SWI2 line 2325)
                    if (get<t_meter, Head>() >= getBottom() and
                        getAt<t_meter, Head>(topNodeID) >= at(topNodeID)->getBottom()) {

                        for (large_num zetaID = 1; zetaID < zetasSize() - 1; zetaID++) {
                            // zeta only moves through the top of a node if there is a ZETA surface
                            // - at the top of the current node (in SWI2: IPLPOS_(i,j,k,n) = 1)
                            // - AND at the bottom of the top node (in SWI2: IPLPOS_(i,j,k-1,n) = 2)
                            if (isZetaAtElevation(zetaID) and at(topNodeID)->isZetaAtBottom(zetaID) ) {

                                // if vertical flux through the top of the node is positive...
                                if (fluxCorrectionTop.value() > 0) {
                                    deltaZeta = (fluxCorrectionTop * (day * get<t_dim, StepSize>())) /
                                                (get<t_s_meter, Area>() * getAt<t_dim, EffectivePorosity>(topNodeID));
                                    // ...lift zeta height of the lowest zeta surface in top node
                                    t_meter zeta_back_top = at(topNodeID)->getZetas().back();
                                    at(topNodeID)->setZeta(zetaID, zeta_back_top + deltaZeta);
                                // if vertical flux through the top of the node is negative...
                                } else if (fluxCorrectionTop.value() < 0) {
                                    deltaZeta = (fluxCorrectionTop * (day * get<t_dim, StepSize>())) /
                                                (get<t_s_meter, Area>() * getEffectivePorosity());
                                    //LOG(debug) << "deltaZeta: " << deltaZeta.value() << std::endl;
                                    // ...lower zeta height of this zeta surface by delta zeta
                                    setZeta(zetaID, getZetas().front() + deltaZeta);
                                }
                            }
                        }
                    }
                }
            }

            /**
             * @brief Horizontal tracking of density interface movement between horizontal neighbour nodes
             * @note in SWI2 documentation: "Tip and Toe Tracking", in SWI2 code: SSWI2_HORZMOVE
             */
            void horizontalZetaMovement(){
                t_meter maxDelta; // maximum "allowed" height difference of a zeta surface n between adjacent nodes
                t_dim factor_self; // potential zeta height factor for this node
                t_dim factor_neig; // potential zeta height factor for neighbour node
                t_meter delta_opp; // potential zeta height adjustment for opposite neighbour node
                t_meter length_self;
                t_meter length_neig;
                t_meter zeta_self;
                t_meter zeta_neig;
                if (getHead() < getBottom()) { return; }

                for (large_num zetaID = 1; zetaID < zetasSize() - 1; ++zetaID) {
                    if (isZetaActive(zetaID)) {
                        for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                            if (at(neigNodeID)->isZetaInactive(zetaID)) {
                                if (at(neigNodeID)->getEffectivePorosity().value() > 0) {
                                    length_self = getNodeLength(neigPos);
                                    length_neig = getLengthNeig(neigPos, neigNodeID);
                                    // if tracking tip/toe: raise/lower this zeta surface in this node by:
                                    factor_self = get<t_dim, SlopeAdjFactor>() *
                                                  ((at(neigNodeID)->getEffectivePorosity() * length_neig) /
                                                   ((getEffectivePorosity() * length_self) +
                                                    (at(neigNodeID)->getEffectivePorosity() * length_neig)));
                                    // if tracking tip/toe: lower/raise this zeta surface in neighbouring node by:
                                    factor_neig = get<t_dim, SlopeAdjFactor>() *
                                                  ((getEffectivePorosity() * length_self) /
                                                   ((getEffectivePorosity() * length_self) +
                                                    (at(neigNodeID)->getEffectivePorosity() * length_neig)));

                                    zeta_self = getZeta(zetaID);
                                    zeta_neig = at(neigNodeID)->getZeta(zetaID);

                                    if (at(neigNodeID)->isZetaAtBottom(zetaID) and
                                        neighbours.find(DOWN) == neighbours.end()) {
                                        //%% Toe tracking %%
                                        maxDelta = 0.5 * (length_self + length_neig) * get<t_dim, MaxToeSlope>();
                                        if (zeta_self > zeta_neig + maxDelta) {
                                            setZeta(zetaID, zeta_self - (factor_self * maxDelta));
                                            at(neigNodeID)->setZeta(zetaID, zeta_neig + (factor_neig * maxDelta));
                                        }
                                    } else if (at(neigNodeID)->isZetaAtFront(zetaID) and
                                               neighbours.find(TOP) == neighbours.end()) {
                                        //%% Tip tracking %%
                                        maxDelta = 0.5 * (length_self + length_neig) * get<t_dim, MaxTipSlope>();
                                        if (zeta_neig > zeta_self + maxDelta) {
                                            setZeta(zetaID, zeta_self + (factor_self * maxDelta));
                                            at(neigNodeID)->setZeta(zetaID, zeta_neig - (factor_neig * maxDelta));
                                        }
                                    }

                                    maxDelta = 0.5 * (length_self + length_neig) * get<t_dim, MaxToeSlope>();
                                    if ((zeta_self - getZetas().back()) < (get<t_dim, MinDepthFactor>() * factor_self * maxDelta)) {
                                        NeighbourPosition oppNeigPos = getOppositePosition(neigPos);
                                        if (neighbours.find(oppNeigPos) != neighbours.end()){
                                            large_num nodeIDOppNeig = neighbours[oppNeigPos];
                                            if (at(nodeIDOppNeig)->isZetaActive(zetaID)) {
                                                // change zeta in other direction neighbour
                                                delta_opp = ((zeta_self - getZetas().back()) *
                                                             (getNodeLength(oppNeigPos) * getEffectivePorosity()) /
                                                             (getLengthNeig(oppNeigPos, nodeIDOppNeig) *
                                                              at(nodeIDOppNeig)->getEffectivePorosity()));
                                                t_meter zeta_opp = at(nodeIDOppNeig)->getZeta(zetaID);
                                                at(nodeIDOppNeig)->setZeta(zetaID, zeta_opp + delta_opp);
                                                setZeta(zetaID, getZetas().back());
                                                //LOG(debug) << "delta_opp (toe): " << delta_opp.value() << std::endl;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    int zoneID = zetaID;
                    // allow intrusion at the coast if...
                    if (hasGHB()){ // (1) node has GHB
                        if (getSourceZoneGHB() == zoneID) { // (2) source zone of GHB is this zone
                            if (isZetaAtBottom(zetaID)) { // (3) zeta with current zetaID is at bottom
                                if (getEffectivePorosity().value() > 0) { // (4) effective porosity is above 0
                                    zeta_self = getZeta(zetaID);
                                    maxDelta = 0.5 * std::sqrt(getArea().value()) * si::meter * get<t_dim, MaxTipSlope>();
                                    if (externalFlows.at(GENERAL_HEAD_BOUNDARY).getFlowHead() > zeta_self + maxDelta) {
                                        factor_self = get<t_dim, SlopeAdjFactor>();
                                        // raise this zeta surface in this node by:
                                        t_meter newZeta = getZetas().back() + (factor_self * maxDelta);
                                        //LOG(debug) << "zetaBack = " << getZetas().back().value() << ", newZeta = " << newZeta.value();
                                        setZeta(zetaID, newZeta);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /**
             * @brief Get the neighbour position of the opposite neighbour (e.g., RIGHT -> LEFT)
             * @param position
             * @return
             */
            NeighbourPosition getOppositePosition(NeighbourPosition position) {
                if (position == NeighbourPosition::BACK) {
                    return NeighbourPosition::FRONT;
                } else if (position == NeighbourPosition::FRONT) {
                    return NeighbourPosition::BACK;
                } else if (position == NeighbourPosition::LEFT) {
                    return NeighbourPosition::RIGHT;
                } else if (position == NeighbourPosition::RIGHT) {
                    return NeighbourPosition::LEFT;
                } else {
                    LOG(userinfo) << "Position unavailable for function getOppositePosition";
                    throw "Position unavailable for function getOppositePosition";
                }
            }

            /**
             * @brief clip zeta surfaces that are outside of the bounds (bounds: first and last zeta surface)
             * @param zetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ZETACLIP
             */
            void clipInnerZetas() {
                for (large_num zetaID = 1; zetaID < zetasSize() - 1; zetaID++) {
                    if (isZetaActive(zetaID)) {
                        if (getZeta(zetaID) < getZetas().back()) { setZeta(zetaID, getZetas().back()); }
                        if (getZeta(zetaID) > getZetas().front()) { setZeta(zetaID, getZetas().front()); }
                    }
                }
            }

            /**
             * @brief Avoid zeta surface locking
             * @param zetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ANTILOCKMIN
             */
            void preventZetaLocking(){
                t_meter maxDelta;
                t_meter maxAdjustment_self;
                t_meter maxAdjustment_neig;
                t_meter delta_self;
                t_meter delta_neig;
                t_meter zeta_self;
                t_meter zeta_neig;
                t_meter zeta_self_new;
                t_meter zeta_neig_new;
                for (large_num zetaID = 1; zetaID < zetasSize() - 1; zetaID++) {
                    if (isZetaInactive(zetaID) and getEffectivePorosity().value() > 0) {
                        // iterate through horizontal neighbours
                        for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                            if (at(neigNodeID)->getEffectivePorosity().value() > 0) {
                                // determine max delta zeta
                                maxAdjustment_self = get<t_meter, VerticalSize>() * get<t_dim, SlopeAdjFactor>();
                                maxAdjustment_neig = getAt<t_meter, VerticalSize>(neigNodeID) * get<t_dim, SlopeAdjFactor>();
                                maxDelta = std::min(std::min(maxAdjustment_self, maxAdjustment_neig), getVDFLock());

                                // calculate delta zeta of node and neighbour
                                delta_self = maxDelta * (at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID)) /
                                             (getEffectivePorosity() * getNodeLength(neigPos) +
                                              at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID));
                                delta_neig = maxDelta * (getEffectivePorosity() * getNodeLength(neigPos)) /
                                             (getEffectivePorosity() * getNodeLength(neigPos) +
                                              at(neigNodeID)->getEffectivePorosity() * getLengthNeig(neigPos, neigNodeID));

                                zeta_self = getZeta(zetaID);
                                zeta_neig = at(neigNodeID)->getZeta(zetaID);
                                // if a zeta surface is at the BOTTOM of this node and at the TOP of the neighbour
                                if (isZetaAtBottom(zetaID) and
                                        (at(neigNodeID)->isZetaAtFront(zetaID) and
                                         !at(neigNodeID)->isZetaAtBottom(zetaID))) {
                                    zeta_self_new = zeta_self + delta_self;
                                    zeta_neig_new = zeta_neig - delta_neig;
                                    if (zeta_self_new < zeta_neig_new) {
                                        // adjust zeta surface heights
                                        setZeta(zetaID, zeta_self_new);
                                        at(neigNodeID)->setZeta(zetaID, zeta_neig_new);
                                    }

                                // else if a zeta surface is at the TOP of this node and at the BOTTOM of the neighbour
                                } else if (at(neigNodeID)->isZetaAtBottom(zetaID) and
                                            (isZetaAtFront(zetaID) and
                                             !isZetaAtBottom(zetaID))) {
                                    zeta_self_new = zeta_self - delta_self;
                                    zeta_neig_new = zeta_neig + delta_neig;
                                    if (zeta_self_new > zeta_neig_new) {
                                        // adjust zeta surface heights
                                        setZeta(zetaID, zeta_self_new);
                                        at(neigNodeID)->setZeta(zetaID, zeta_neig_new);
                                    }
                                }
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
                    LOG(userinfo) << "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
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
                    LOG(userinfo) << "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
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
                    LOG(userinfo) << "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            bool isLeftOrRight(NeighbourPosition neig) {
                std::unordered_map<NeighbourPosition, bool> leftRightMap = {
                        {LEFT, true}, {RIGHT, true}
                };
                try{
                    return leftRightMap.at(neig);
                } catch (const std::out_of_range &e) { return false;}
            }

            bool isFrontOrBack(NeighbourPosition neig) {
                std::unordered_map<NeighbourPosition, bool> frontBackMap = {
                        {FRONT, true}, {BACK, true}
                };

                try{
                    return frontBackMap.at(neig);
                } catch (const std::out_of_range &e) { return false;}
            }

            /**
             * @brief Get P part of flow equations (for left side of the equation)
             * @return volume over time
             */
            t_s_meter_t getP_aboveFlowBottom() noexcept {
                t_s_meter_t ex = 0.0 * (si::square_meter / day);
                auto eq_head = get<t_meter, EQHead>();
                auto head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                if (hasTypeOfExternalFlow(RECHARGE)) {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                }
                t_vol_t eqFlow = getEqFlow();
                for (const auto &extFlow : externalFlows) {
                    if (is(extFlow.second.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                        // if flow depends on head
                        if (extFlow.second.isFlowHeadDependent(head)) {
                            ex += extFlow.second.getP(eq_head, head, recharge, eqFlow);
                        }
                    // extFlow.second.getP() with flowType RECHARGE returns 0
                    } else if (extFlow.second.getType() == GENERAL_HEAD_BOUNDARY){
                        ex += extFlow.second.getP(eq_head, head, recharge, eqFlow);
                    }
                }
                return ex;
            }

            /**
             * @brief Get P part of RHS
             * @return volume over time
             * @note The source zones of GHB and GWR are specified in config and/or in file, sink of GHB is the top zone
             */
            t_s_meter_t getP_aboveFlowBottom(large_num zetaID) noexcept {
                large_num densityZone = zetaID; // the density zone has the same ID as the zeta interface
                t_s_meter_t ex = 0.0 * (si::square_meter / day);
                t_s_meter_t P = 0.0 * (si::square_meter / day);

                auto eq_head = get<t_meter, EQHead>();
                auto head = get<t_meter, Head>(); // need the current head
                t_vol_t recharge = 0 * si::cubic_meter / day;
                if (hasTypeOfExternalFlow(RECHARGE)) {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                }
                t_vol_t eqFlow = getEqFlow();
                for (const auto &[flowType, extFlow] : externalFlows) {
                    // extFlow.getP with flowType RECHARGE returns 0
                    if (flowType == GENERAL_HEAD_BOUNDARY) {
                        P = extFlow.getP(eq_head, head, recharge, eqFlow);
                        // if GHB flow is into node (gw_head < ghb_elevation) and ...
                        if (getHead() < extFlow.getFlowHead()) {
                            // ... this is the zone of GHB inflow
                            if (densityZone == getSourceZoneGHB()) {
                                ex += P;
                            }
                        // else if GHB flow is out of node (gw_head > ghb_elevation) and ...
                        } else if (getHead() > extFlow.getFlowHead()) {
                            // ... interface above is at front
                            if (!isZetaActive(zetaID - 1)) {
                                ex -= P;
                            }
                        }
                    }
                }
                return ex;
            }


            /*double getFlowFractionOfZone(large_num zetaID, ExternalFlow extFlow) {
                t_meter bottomOfFlowInZone;
                t_meter flowHeightInZone;

                if (getHead() < extFlow.getFlowHead()) { // gw_head is below flow head (-> flow is into node)
                    // if this node has a general head boundary,...
                    if (hasGHB()) {
                        // ... and the river head is below the general head boundary, ...
                        if (extFlow.getFlowHead() < externalFlows.at(GENERAL_HEAD_BOUNDARY).getFlowHead()) {
                            // ... then we assume that the river is fully saline (-> water goes to lowest density zone)
                            // hence, if this is the zeta surface at the top of the most saline zone, ...
                            if (zetaID == getZetas_TZero().size() - 2) {
                                // ... it gets all of that saline water
                                return 1;
                            }
                        }
                    }
                } else { // gw_head is above flow head (-> flow is out of node)
                    t_meter flowHeight = getHead() - extFlow.getBottomElev(); // height of flow is: gw_head - flow_bottom

                    // if this zeta is above the flow bottom (-> water from its zone enters the SWB)
                    if (getZeta_TZero(zetaID) > extFlow.getBottomElev()) {
                        // if the zeta below is higher than the river bottom ...
                        if (getZeta_TZero(zetaID + 1) > extFlow.getBottomElev()) {
                            // ... the bottom of flow into/out of this zone is the lower zeta
                            bottomOfFlowInZone = getZeta_TZero(zetaID + 1);
                            // the zeta below is not higher than the river bottom, ...
                        } else {
                            // ... so the bottom of flow into/out of this zone is the flow bottom
                            bottomOfFlowInZone = extFlow.getBottomElev();
                        }
                        flowHeightInZone = getZeta_TZero(zetaID) - bottomOfFlowInZone;
                        if (flowHeightInZone.value() < 0) {
                            flowHeightInZone = 0 * si::meter;
                        }
                        return flowHeightInZone / flowHeight;
                    }
                }
                return 0;
            }*/

            /**
             * @brief Get flow which is not groundwater head dependent
             * @return volume over time
             * Flow can be added to constant flows on right side of the equations
             * e.g., when GW head is below external flow bottom, we expect drainage through the bottom of external flow
             */
            t_vol_t getP_belowFlowBottom() noexcept {
                auto eq_head = get<t_meter, EQHead>();
                auto head = get<t_meter, Head>();
                t_vol_t recharge = 0 * si::cubic_meter / day;
                if (hasTypeOfExternalFlow(RECHARGE)) {
                    recharge = getExternalFlowByName(RECHARGE).getRecharge();
                }
                t_vol_t eqFlow = getEqFlow();
                t_vol_t ex = 0.0 * (si::cubic_meter / day);
                //Q part is already subtracted in RHS
                for (const auto &flow : externalFlows) {
                    if (is(flow.second.getType()).in(RIVER, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND, DRAIN)) {
                        if (not flow.second.isFlowHeadDependent(head)) {
                            ex += flow.second.getP(eq_head, head, recharge, eqFlow) * flow.second.getBottomElev();
                        }
                    }

                }
                return ex;
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
                NeighbourPosition neigPos;
                large_num neigNodeID;

                //Get all conductances from neighbouring cells
                for (const auto &neighbour: neighbours) {
                    neigPos = neighbour.first;
                    neigNodeID = neighbour.second;
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
                        } else {
                            conduct = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(neigPos, neigNodeID));
                            if (conduct.value() <= 0) {
                                LOG(debug) << "conduct = " << conduct.value();
                            }
                        }
                    }
                    NANChecker(conduct.value(), "Conductances");
                    // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                    // (used to solve for the head at the neighbouring node)
                    out[neigNodeID] = conduct;
                }

                t_s_meter_t conductNode = 0 * si::square_meter / day;

                // To solve for the head at this node, the conductances to neighbours and HCOF are used
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &[_, conductNeig]: out) { conductNode -= conductNeig; }
                // add HCOF
                t_s_meter_t hcof = getP_aboveFlowBottom() - (getStorageCapacity() / (day * get<t_dim, StepSize>()));
                conductNode += hcof;
                //LOG(debug) << "conductNode: " << conductNode.value();
                NANChecker(conductNode.value(), "conductNode");

                // add resulting conductance to solve for the head at this node to out
                out[getID()] = std::move(conductNode);
                return out;
            };

            /**
             * @brief The matrix entry for the left hand side of the zeta surface equation
             * @param zetaID zeta surface id in this node
             * @return map <CellID,Conductance>
             */
            std::unordered_map<large_num, t_s_meter_t> getMatrixEntries(large_num zetaID) {
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(horizontal_neighbours.size()+1);
                auto delnus = get<std::vector<t_dim>, Delnus>();
                std::vector<t_s_meter_t> zoneConductances;
                t_s_meter_t zoneConductanceCum;
                t_s_meter_t zetaMovementConductance;
                std::vector<t_meter> zoneThicknesses;

                // if this zeta interface is not active: return directly
                /*if (isZetaInactive(zetaID)){
                    LOG(userinfo) << "Asking for matrix entry at inactive zeta, causing matrix row of zeros, and unsolvable equation";
                    throw "Asking for matrix entry at inactive zeta, causing matrix row of zeros, and unsolvable equation";
                }*/

                for (auto const &[neigPos, neigNodeID]: horizontal_neighbours) {
                    zetaMovementConductance = 0 * (si::square_meter / day);
                    if (at(neigNodeID)->isZetaTZeroActive(zetaID)) {
                        // on left hand side of equation -> need updated zeta heights
                        zoneThicknesses = calculateZoneThicknesses(neigPos, neigNodeID, getZetas(),
                                                                   at(neigNodeID)->getZetas());
                        zoneConductances = getZoneConductances(neigPos, neigNodeID, zoneThicknesses);
                        zoneConductanceCum = getZoneConductanceCum(zetaID, zoneConductances);
                        zetaMovementConductance += delnus[zetaID] * zoneConductanceCum; // in SWI2: SWISOLCC/R
                        NANChecker(zetaMovementConductance.value(), "zetaMovementConductance");
                        // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                        out[neigNodeID] = zetaMovementConductance;
                    }
                }

                // To solve for zeta in this node, the conductances to neighbours and porosity term are used
                t_s_meter_t conductNode = 0 * (si::square_meter / day);
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &c: out) { conductNode -= c.second; }
                // subtract effective porosity term
                conductNode -= getEffectivePorosityTerm(); // subtracting effective porosity term (SWIHCOF)
                NANChecker(conductNode.value(), "conductNode");
                // add conductance of this node to out, the key in the unordered map is the ID of this node
                out[get<large_num, ID>()] = conductNode;
                return out;
            };

            /**
             * @brief The right hand side of the flow equation
             * @return cubic meters per time
             */
            t_vol_t getRHS() {
                t_vol_t externalSources = - getQ() - getP_belowFlowBottom(); // e.g., recharge, river, lakes, wetlands
                //LOG(debug) << "externalSources: " << externalSources.value() << std::endl;
                t_vol_t dewateredFlow = calculateDewateredFlow(); // only if node has bottom neighbour
                //LOG(debug) << "dewateredFlow: " << dewateredFlow.value() << std::endl;
                t_vol_t storageFlow = -getStorageCapacity() * getHead_TZero() / (day * get<t_dim, StepSize>());
                //LOG(debug) << "storageFlow: " << storageFlow.value() << std::endl;
                t_vol_t internalSources = dewateredFlow + storageFlow;
                t_vol_t out = externalSources + internalSources;
                //LOG(debug) << "RHS constant density: " << out.value() << std::endl;
                NANChecker(out.value(), "RHS constant density");

                if (get<bool, IsDensityVariable>()) {
                    // save constant density RHS (without variable density terms) for calculation of zeta movement
                    //set<t_vol_t, InternalSources>(internalSources);

                    // calculate Pseudo-Source Flow
                    t_vol_t pseudoSourceNode = getPseudoSourceNode();
                    //LOG(debug) << "nodeID: " << getID() << ", pseudoSourceNode: " << pseudoSourceNode.value();
                    // calculate Vertical Flux Correction (from top neighbour)
                    t_vol_t verticalFluxCorrections = getVerticalFluxCorrections();
                    //LOG(debug) << "nodeID: " << getID() << ", verticalFluxCorrections: " << verticalFluxCorrections.value();
                    out += pseudoSourceNode + verticalFluxCorrections;
                }
                //LOG(debug) << "getRHS (nodeID: " << getID() << "): " << out.value() << std::endl;
                NANChecker(out.value(), "RHS");
                return out;
            }

            /**
             * @brief calculate the right hand side for zeta surface equation (b_zetas)
             * @param zetaID zeta surface id in this node
             * @return cubic meters per time
             */
            t_vol_t getRHS(large_num zetaID){
                t_vol_t out = 0 * si::cubic_meter / day;
                if (isZetaTZeroInactive(zetaID)) { return out; } // (line 3571)

                t_vol_t porosityTerm = getEffectivePorosityTerm() * getZeta_TZero(zetaID);
                //LOG(debug) << "porosityTerm: " << porosityTerm.value();
                t_vol_t pseudoSourceBelowZeta = getPseudoSourceBelowZeta(zetaID); // in SWI2 code: SSWI2_SD and SSWI2_SR
                //LOG(debug) << "pseudoSourceBelowZeta: " << pseudoSourceBelowZeta.value();

                // in SWI2 code: part of BRHS; in SWI2 doc: G or known source term below zeta
                t_vol_t externalFlow = 0.0 * (si::cubic_meter / day);
                t_vol_t storageChange = 0.0 * (si::cubic_meter / day);
                // if the new groundwater head is above or equal to the node bottom
                if (getHead() >= getBottom()) { // lines 3532-3536
                    externalFlow = - getQ(zetaID) - getP_aboveFlowBottom(zetaID) * getHead();
                    //LOG(debug) << "externalFlow: " << externalFlow.value();
                    storageChange = getZoneFraction(zetaID) * getStorageCapacity() * (getHead() - getHead_TZero()) /
                                            (day * get<t_dim, StepSize>());
                    //LOG(debug) << "storageChange: " << storageChange.value();

                }
                t_vol_t tipToeFlow = getTipToeFlow(zetaID); // in SWI2 code: SSWI2_QR and SSWI2_QC
                //LOG(debug) << "tipToeFlow: " << tipToeFlow.value();

                out = - porosityTerm + pseudoSourceBelowZeta + externalFlow + storageChange + tipToeFlow;
                NANChecker(out.value(), "getRHS(large_num zetaID)");
                return out;
            }

            /**
             * @brief set both head and the head change resulting from a solver iteration
             * @param change current head change
             */
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
                t_vol_t dewateredFlow = -calculateDewateredFlow() * get<t_dim, StepSize>();
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
                                    useEfolding, confined, isSteadyState, isDensityVariable, delnus,
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
                           const t_s_meter& area,
                           const t_meter& edgeLengthLeftRight,
                           const t_meter& edgeLengthFrontBack)
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
                    0.000015, false, true, true, false, {0.0, 0.1}, {0.0, 0.1},
                    0.2, 0.2, 0.2, 0.1, 0.1, 0.001 * si::meter, 0, 0) {}

        private:
            friend class NodeInterface;

            virtual void __setHeadChange(const t_meter& change) {
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