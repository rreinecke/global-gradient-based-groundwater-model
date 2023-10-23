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
            virtual void
            __setHeadChange(t_meter change) = 0;

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
            std::unordered_map<NeighbourPosition, large_num> neighbours;
            std::unordered_map<FlowType, ExternalFlow, FlowTypeHash> externalFlows;
            int numOfExternalFlows{0};
            bool initial_head{true};
            bool simpleDistance{false};
            bool steadyState{false};
            int zoneOfSinks; // individual for each node
            int zoneOfSources; // individual for each node

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

            p_node & at(map_itter pos) { return nodes->at(pos->second); }

            template<typename T, typename F>
            T getAt(map_itter pos) {
                return at(pos)->get<T, F>();
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
                t_vol_t storageFlow = 0.0 * si::cubic_meter / day;
                if (not steadyState) {
                    storageFlow = getTotalStorageFlow();
                    if (compare(storageFlow.value()))
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
                ar & steadyState;
                ar & fields;
                ar & cached;
                ar & eq_flow;
            }

            /**
             * @brief Overload ostream operator "<<" to return various node properties
             */
            friend std::ostream & operator<< (std::ostream & stream, const p_node & pNode) {
                stream << "Node properties:"
                       << "\nID [large_num]: "<< pNode->get<large_num, ID>()
                       << "\nSpatID [large_num]: " << pNode->get<large_num, SpatID>()
                       << "\nRefID [int]: " << pNode->get<large_num, RefID>()<< "\nLat [double]: " << pNode->get<double, Lat>()
                       << "\nLon [double]: " << pNode->get<double, Lon>()
                       << "\nLayer [int]: " << pNode->get<int, Layer>()
                       << "\nStepModifier [-]: " << pNode->get<t_dim, StepModifier>().value()
                       << "\nArea [m²]: " << pNode->get<t_s_meter, Area>().value()
                       << "\nVerticalSize [m]: " << pNode->get<t_meter, VerticalSize>().value()
                       << "\nElevation [m]: " << pNode->get<t_meter, Elevation>().value()
                       << "\nTopElevation [m]: " << pNode->get<t_meter, TopElevation>().value()
                       << "\nEFolding [m]: " << pNode->get<t_meter, EFolding>().value()
                       << "\nConfinement [bool]: " << pNode->get<bool, Confinement>()
                       << "\nK [m/s]: " << pNode->get<t_vel, K>().value()
                       << "\nAnisotropy [-]: " << pNode->get<t_dim, Anisotropy>().value()
                       << "\nOUT [m³]: " << pNode->get<t_c_meter, OUT>().value()
                       << "\nIN [m³]: " << pNode->get<t_c_meter, IN>().value()
                       << "\nHead [m]: " << pNode->get<t_meter, Head>().value()
                       << "\nEQHead [m]: " << pNode->get<t_meter, EQHead>().value()
                       << "\nSpecificYield [-]: " << pNode->get<t_dim, SpecificYield>().value()
                       << "\nSpecificStorage [perUnit]: " << pNode->get<quantity < perUnit>, SpecificStorage>().value()
                        << "\nEdgeLengthLeftRight [m]: " << pNode->get<t_meter, EdgeLengthLeftRight>().value()
                        << "\nEdgeLengthFrontBack [m]: " << pNode->get<t_meter, EdgeLengthFrontBack>().value()
                        << "\nSurfaceLeftRight [m²]: " << pNode->get<t_s_meter, SurfaceLeftRight>().value()
                        << "\nSurfaceFrontBack [m²]: " << pNode->get<t_s_meter, SurfaceFrontBack>().value()
                        << "\nEffectivePorosity [-]: " << pNode->get<t_dim, EffectivePorosity>().value()
                        << "\nDensityVariable [-]: " << pNode->get<bool, DensityVariable>()
                        << "\nMaxTipSlope [-]: " << pNode->get<t_dim, MaxTipSlope>().value()
                        << "\nMaxToeSlope [-]: " << pNode->get<t_dim, MaxToeSlope>().value()
                        << "\nSlopeAdjFactor [-]: " << pNode->get<t_dim, SlopeAdjFactor>().value()
                        << "\nVDFLock [-]: " << pNode->get<t_meter, VDFLock>().value();

                std::unordered_map<NeighbourPosition, large_num> neighbourList = pNode->getListOfNeighbours();
                if(neighbourList.find(DOWN) != neighbourList.end()) {
                    stream << "\nDOWN neighbour lat [double]: " << pNode->getNeighbour(DOWN)->get<double, Lat>()
                           << "\nDOWN neighbour lon [double]: " << pNode->getNeighbour(DOWN)->get<double, Lon>();
                }
                if(neighbourList.find(TOP) != neighbourList.end()) {
                    stream << "\nTOP neighbour lat [double]: " << pNode->getNeighbour(TOP)->get<double, Lat>()
                           << "\nTOP neighbour lon [double]: " << pNode->getNeighbour(TOP)->get<double, Lon>();
                }
                if(neighbourList.find(LEFT) != neighbourList.end()) {
                    stream << "\nLEFT neighbour lat [double]: " << pNode->getNeighbour(LEFT)->get<double, Lat>()
                           << "\nLEFT neighbour lon [double]: " << pNode->getNeighbour(LEFT)->get<double, Lon>();
                }
                if(neighbourList.find(RIGHT) != neighbourList.end()) {
                    stream << "\nRIGHT neighbour lat [double]: " << pNode->getNeighbour(RIGHT)->get<double, Lat>()
                           << "\nRIGHT neighbour lon [double]: " << pNode->getNeighbour(RIGHT)->get<double, Lon>();
                }
                if(neighbourList.find(FRONT) != neighbourList.end()) {
                    stream << "\nFRONT neighbour lat [double]: " << pNode->getNeighbour(FRONT)->get<double, Lat>()
                           << "\nFRONT neighbour lon [double]: " << pNode->getNeighbour(FRONT)->get<double, Lon>();
                }
                if(neighbourList.find(BACK) != neighbourList.end()) {
                    stream << "\nBACK neighbour lat [double]: " << pNode->getNeighbour(BACK)->get<double, Lat>()
                           << "\nBACK neighbour lon [double]: " << pNode->getNeighbour(BACK)->get<double, Lon>();
                }

                if (pNode->hasGHB()) {
                    stream << "\nGHB conductance: "
                           << pNode->getExternalFlowConductance(Model::GENERAL_HEAD_BOUNDARY).value();
                    stream << "\nGHB elevation: "
                           << pNode->getExternalFlowElevation(Model::GENERAL_HEAD_BOUNDARY).value();
                }

                if (pNode->hasTypeOfExternalFlow(RECHARGE)){
                    stream << "\nRecharge: "
                           << pNode->getExternalFlowByName(Model::RECHARGE).getRecharge().value();
                }

                stream << "\nZone of sinks: " << pNode->zoneOfSinks;
                stream << "\nZone of sources: " << pNode->zoneOfSources;

                return stream;
            };

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
                          bool hasRefinedNeighbour,
                          bool densityVariable,
                          std::vector<t_dim> delnus,
                          std::vector<t_dim> nusInZones,
                          double effPorosity,
                          double maxTipSlope,
                          double maxToeSlope,
                          double minDepthFactor,
                          double slopeAdjFactor,
                          t_meter vdfLock);

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
                    t_vel out = get<t_vel, K>() * e_fold * get<t_dim, StepModifier>();
                    if (out < 1e-20 * si::meter / day) {
                        out = 1e-20 * si::meter / day;
                    }
                    return out;
                } else {
                    return get<t_vel, K>() * get<t_dim, StepModifier>();
                }
            }

/*****************************************************************
Set Properties
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
             * @brief Update the current head change (in comparison to last time step)
             * @note Should only be called at end of time step
             */
            void updateHeadChange() noexcept {
                set < t_meter, HeadChange_TZero > (get<t_meter, Head>() - get<t_meter, Head_TZero>());
                set < t_meter, Head_TZero > (get<t_meter, Head>());
            }

            void initHead_t0() noexcept { set < t_meter, Head_TZero > (get<t_meter, Head>()); }

            void setHead(t_meter head) noexcept {
                NANChecker(head.value(), "Set Head");
                set<t_meter, Head>(head);
            }

            void setHead_direct(double head) noexcept { set < t_meter, Head > (head * si::meter); }

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
            FlowInputHor createDataTuple(map_itter got) {
                return std::make_tuple(at(got)->getK(),
                                       getK(),
                                       at(got)->getNodeLength(got), // length of neighbour node (parallel to direction)
                                       getNodeLength(got), // length of this node (parallel to direction)
                                       std::min(getNodeWidth(got), at(got)->getNodeWidth(got)), // width of smaller node (perpendicular to direction)
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
                                       getAt<t_meter, VerticalSize>(got),
                                       get<t_meter, Head>(),
                                       getAt<t_meter, Head>(got),
                                       get<t_meter, Elevation>(),
                                       getAt<t_meter, Elevation>(got),
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
                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();

                for (const auto &position: possible_neighbours) {
                    auto got = neighbours.find(position);
                    if (got != neighbours.end()) { //There is a neighbour node
                        t_s_meter_t conductance;

                        if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                            conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(got), get<t_meter, EFolding>(), getAt<t_meter, EFolding>(got));
                        } else {
                            conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got));
                        }

                        t_vol_t flow = conductance * (get<t_meter, HeadType>() - getAt<t_meter, HeadType>(got));

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
             * @return lateral flows
             */
            t_vol_t getLateralFlows() {
                return calcLateralFlows<Head>(false) * get<t_dim, StepModifier>();
            }

            /**
             * Get the current lateral out flows
             * @return lateral outflows
             */
            t_vol_t getLateralOutFlows() {
                return calcLateralFlows<Head>(true) * get<t_dim, StepModifier>();
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
            void setK(t_vel conduct) {
                setK_direct(conduct);
                applyToAllLayers([&conduct](NodeInterface *nodeInterface) {
                    nodeInterface->getProperties().set<t_vel, K>(conduct);
                });
            }

            /**
             * @brief Modify hydraulic conductivity (no e-folding, no layers)
             * @param conduct conductivity
             */
            void setK_direct(t_vel conduct) { set < t_vel, K > (conduct); }


            int getSpatID() {return (int) get<large_num, SpatID>();}

            double getLat() {return get<double, Lat>();}

            double getLon() {return get<double, Lon>();}

            large_num getRefID() {return get<large_num, RefID>(); }

            t_s_meter getArea(){return get<t_s_meter, Area>();}

            t_dim getEffectivePorosity(){return get<t_dim, EffectivePorosity>();}

            t_meter getEdgeLengthLeftRight(){return get<t_meter, EdgeLengthLeftRight>();}

            t_meter getEdgeLengthFrontBack(){return get<t_meter, EdgeLengthFrontBack>();}

            t_meter getElevation(){return get<t_meter, Elevation>();}

            t_meter getBottom(){return get<t_meter, Elevation>() - get<t_meter, VerticalSize>();}

            t_meter getHead(){ return get<t_meter, Head>(); }

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

            /**
             * @brief Toggle steady state simulation
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
            t_s_meter_t getExternalFlowConductance(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type).getConductance();
                } else {
                    return 0 * si::square_meter / day;
                }
            }

            /**
             * @brief Get head of an external flow by its FlowType
             * @param type The flow type
             * @return Head
             */
            t_meter getExternalFlowElevation(FlowType type) {
                if (hasTypeOfExternalFlow(type)) {
                    return externalFlows.at(type).getFlowHead();
                } else {
                    return 0 * si::meter;
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
                t_vol_t eqFlow = getEqFlow();
                if (is(flow.getType()).in(RIVER, DRAIN, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND)) {
                    if (flow.flowIsHeadDependent(head)) {
                        ex = (flow.getP(eq_head, head, recharge, eqFlow) * head +
                              flow.getQ(eq_head, head, recharge, eqFlow)) * get<t_dim, StepModifier>();
                    } else { // flow is not head dependent when the head is below the bottom of the simulated cell
                        ex = (flow.getP(eq_head, head, recharge, eqFlow) * flow.getBottom() +
                              flow.getQ(eq_head, head, recharge, eqFlow)) * get<t_dim, StepModifier>();
                    }
                } else {  // GENERAL_HEAD_BOUNDARY (Question: what about FLOODPLAIN_DRAIN, EVAPOTRANSPIRATION, FAST_SURFACE_RUNOFF)
                    ex = (flow.getP(eq_head, head, recharge, eqFlow) * head +
                          flow.getQ(eq_head, head, recharge, eqFlow)) * get<t_dim, StepModifier>();
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
                auto hasDown = neighbours.find(DOWN);
                auto hasUp = neighbours.find(TOP);
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
            t_c_meter getCurrentIN() noexcept { return getFlow([](double a) -> bool { return a > 0; }) * day;}

            /**
             * @brief Get all current OUT flow
             * @return Flow volume
             */
            t_c_meter getCurrentOUT() noexcept { return -getFlow([](double a) -> bool { return a < 0; }) * day; }

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
                if (get<t_meter, Head_TZero>() < getBottom()) {
                    // set zetas to bottom of node (zone thickness = zero)
                    zetaOld = getBottom();
                    zetaBelowOld = getBottom();
                }

                // calculate zone change ((zeta-zetaBelow) * Area is the current zone volume, with old zetas old vol)
                // not dividing by (day * step size) since we want to get the volumetric budget over full time step
                out += ((getEffectivePorosity() * get<t_s_meter, Area>()) *
                        ((zeta - zetaBelow) - (zetaOld - zetaBelowOld)));
                return out;
            }


            /** // todo debug
             * @brief Instantaneous mixing of water, as described in SWI2 documentation under "Vertical Leakage Between
             * Aquifers": (2) when freshwater leaks up into an aquifer containing only saline water, that freshwater is
             *                added as saline water
             *            (4) when saline water leaks down into an aquifer containing only freshwater, that saline water
             *                is added as freshwater
             * @note in SWI2 code: SSWI2_IMIX
             */
            t_c_meter calculateInstantaneousMixing(int localZetaID) {
                t_c_meter out = 0.0 * si::cubic_meter;
                auto down = neighbours.find(NeighbourPosition::DOWN);

                // skip nodes that do not have a down neighbour
                if (down != neighbours.end()) {
                    // skip if dimensionless density at bottom of this node is below or equal to
                    // dimension less density at the top of down node
                    if (getNusBot() <= at(down)->getNusTop()){ return out; } // return 0

                    // skip if head in this or neighbour node is below node bottom
                    if (getHead() < getBottom() or at(down)->getHead() < at(down)->getBottom()){ return out; } // return 0

                    // skip if localZetaID is not involved
                    // todo find the zone number at the bottom of the top node and the top of the down node
                    // if (localZetaID != getZetaID(getNusBot())) { return out; }
                    // if (localZetaID != getZetaID(at(down)->getNusTop())) { return out; }

                    // calculate the flux down
                    t_c_meter fluxDown = getFluxDown().value() * si::cubic_meter;

                    // if flux down is positive
                    if (fluxDown > (0 * si::cubic_meter)) {
                        // skip if dim-less density at bottom of down neighbour is greater or equal to
                        // dim-less density at the bottom of this node
                        if (at(down)->getNusBot() >= getNusBot()) { return out; } // return 0

                        // skip if dim-less density at top of this node is smaller or equal to
                        // dim-less density at the top of neighbour node
                        if (getNusTop() <= at(down)->getNusTop()){ return out; } // return 0
                    }
                    if (isZetaAtBottom(localZetaID) and !isZetaAtBottom(localZetaID - 1)) { out += fluxDown; }
                }
                return out;
            }

            /**
             * @brief save volumetric density zone change between last and new time step
             */
            void saveZoneChange() noexcept {
                set<t_c_meter, ZCHG_IN>(getZoneChange(true));
                set<t_c_meter, ZCHG_OUT>(getZoneChange(false));
            }

            t_c_meter getZoneChange(bool in) noexcept {
                t_c_meter zoneChange_in;
                t_c_meter zoneChange_out;
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; ++localZetaID) {
                    t_c_meter zoneChange = calculateZoneChange(localZetaID);
                    if (zoneChange.value() > 0) { zoneChange_in += zoneChange; } else { zoneChange_out += zoneChange; }
                }
                if (in) { return zoneChange_in; } else { return zoneChange_out; }
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

            t_c_meter getTipToeTrackingZoneChange(bool in) {
                t_c_meter result = 0 * si::cubic_meters;
                t_c_meter tttOut = getZoneChange(false) - get<t_c_meter, ZCHG_OUT>();
                t_c_meter tttIn = getZoneChange(true) - get<t_c_meter, ZCHG_IN>();
                if (in) {
                    if (tttOut.value() > 0) { result += tttOut; }
                    if (tttIn.value() > 0) { result += tttIn; }
                } else {
                    if (tttOut.value() < 0) { result += tttOut; }
                    if (tttIn.value() < 0) { result += tttIn; }
                }
                return result;
            }

            /**
             * @brief save variable density flow mass balance in node property (called before and after adjustment)
             * @param
             */
            void saveVDFMassBalance() noexcept {
                t_c_meter vdfOut = getVDF_OUT();
                t_c_meter vdfIn = getVDF_IN();

                // add current zone change before tip toe tracking
                vdfOut += getZCHG_OUT();
                vdfIn += getZCHG_IN();

                // add current instantaneous mixing budget
                vdfOut += getInstantaneousMixing(false);
                vdfIn += getInstantaneousMixing(true);

                vdfOut += getTipToeTrackingZoneChange(false);
                vdfIn += getTipToeTrackingZoneChange(true);

                set<t_c_meter, VDF_OUT>(vdfOut);
                set<t_c_meter, VDF_IN>(vdfIn);
            }

            t_c_meter getZCHG_OUT() { return get<t_c_meter, ZCHG_OUT>(); }

            t_c_meter getZCHG_IN() { return get<t_c_meter, ZCHG_IN>(); }

            t_c_meter getVDF_OUT() { return get<t_c_meter, VDF_OUT>(); }

            t_c_meter getVDF_IN() { return get<t_c_meter, VDF_IN>(); }

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
            void setNeighbour(large_num ID, NeighbourPosition neighbourPosition) { neighbours[neighbourPosition] = ID; }

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
                    neighbours[neighbourPosition] = nodeIDs_neig.at(0);
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
                            neighbours[refNeigPos] = nodeID_neig.second;
                        } catch (const std::out_of_range &ex) {
                            continue;
                        }
                    }

                    set<bool, HasRefinedNeighbour>(true);
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

            int getNumofNeighbours() { return (int) neighbours.size(); }

            class NodeNotFoundException : public std::exception {
                virtual const char *what() const throw() { return "Node does not exist"; }
            };
            std::unordered_map<NeighbourPosition, large_num> getListOfNeighbours(){
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
                    //LOG(debug) << "Adding flow " << std::to_string(type) << " that already existed at nodeID: " << getID();
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
                return -calculateGhostNodeCorrection(getPossibleRefinedNeighbours());
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
                std::forward_list<NeighbourPosition> neigPos_LRFB = getNeigPos_LRFB();

                for (const auto &neigPos: neigPos_LRFB) {
                    auto neig = neighbours.find(neigPos);
                    if (neig != neighbours.end()) {
                        // if neighbour is unrefined
                        if (at(neig)->getRefinedInto() == 1) { // if neig is not refined
                            // get the refined neighbour position this node has relative to that unrefined neighbour
                            NeighbourPosition refNeigPos = getRefNeigPosToUnrefNeig(get<large_num, RefID>(), neigPos,
                                                                                    getRefinedInto());
                            //LOG(debug) << "getGNCToRefinedNode, nodeID: " << getID() << ", neig nodeID: " << at(neig)->getID();

                            out += at(neig)->calculateGhostNodeCorrection({refNeigPos});
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
                                        {NeighbourPosition::RIGHT,  NeighbourPosition::LEFTFRONT}} },
                            {5, { {NeighbourPosition::LEFT,   NeighbourPosition::RIGHTRIGHTFRONT} } },
                            {6, { {} } }, // has only refined neighbours
                            {7, { {} } }, // has only refined neighbours
                            {8, { {NeighbourPosition::RIGHT,  NeighbourPosition::LEFTLEFTFRONT} } },
                            {9, { {NeighbourPosition::LEFT,   NeighbourPosition::RIGHTRIGHTBACK} } },
                            {10, { {} } }, // has only refined neighbours
                            {11, { {} } }, // has only refined neighbours
                            {12, { {NeighbourPosition::RIGHT, NeighbourPosition::LEFTLEFTBACK} } },
                            {13, { {NeighbourPosition::BACK, NeighbourPosition::FRONTLEFT},
                                        {NeighbourPosition::LEFT, NeighbourPosition::RIGHTBACK}} },
                            {14, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTFRONTLEFT} } },
                            {15, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTFRONTRIGHT} } },
                            {16, { {NeighbourPosition::BACK,  NeighbourPosition::FRONTRIGHT},
                                        {NeighbourPosition::RIGHT, NeighbourPosition::LEFTBACK} } }};
                }
                return mapRefIdToNeigToRefNeig.at(refID).at(neigPos);
            }

            t_vol_t calculateGhostNodeCorrection(const std::forward_list<NeighbourPosition>& possibleRefNeigPos){
                t_vol_t gnc;
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                t_meter head = getHead();
                t_dim multiplierContributor{};
                t_dim multiplierNodeInner{};
                t_dim multiplierNodeOuter{};
                t_s_meter_t transmissivitySelf = get<t_meter, VerticalSize>() * getK();
                t_s_meter_t contributorConductance = 0.0 * (si::square_meter / day);
                t_s_meter_t nodeConductance = 0.0 * (si::square_meter / day);
                std::forward_list<NeighbourPosition> potentialContributors;

                int neigRefInto{};
                int contrRefInto{};
                for (const auto &position: possibleRefNeigPos) {
                    auto refinedNeig = neighbours.find(position);
                    if (refinedNeig != neighbours.end()) {
                        potentialContributors = getPotentialContributors(position);
                        for (const auto &potContrPos: potentialContributors) {
                            auto contributor = neighbours.find(potContrPos); // contributor is named "j" in USG doc
                            if (contributor == neighbours.end()) { // at model boundary: no contributor
                                continue; // todo implement impact of GHB (not required if all nodes at coast/GHB are refined)
                            }

                            contrRefInto = (int) at(contributor)->getRefinedInto();
                            multiplierContributor = ( 1.0 / (2.0 * sqrt(contrRefInto) )) * si::si_dimensionless;

                            neigRefInto = (int) at(refinedNeig)->getRefinedInto();
                            if (neigRefInto == 4) {
                                multiplierNodeInner = (1.0 / 4.0) * si::si_dimensionless;
                                multiplierNodeOuter = (1.0 / 4.0) * si::si_dimensionless;
                            } else if (neigRefInto == 9) {
                                multiplierNodeInner = (2.0 / 6.0) * si::si_dimensionless;
                                multiplierNodeOuter = (1.0 / 6.0) * si::si_dimensionless;
                            } else if (neigRefInto == 16) {
                                if (refinedNeig->first > 18) {
                                    multiplierNodeInner = (1.0 / 8.0) * si::si_dimensionless;
                                    multiplierNodeOuter = (3.0 / 8.0) * si::si_dimensionless;
                                } else {
                                    multiplierNodeInner = (3.0 / 8.0) * si::si_dimensionless;
                                    multiplierNodeOuter = (1.0 / 8.0) * si::si_dimensionless;
                                }
                            }

                            t_s_meter_t transmissivityNeig =
                                    getAt<t_meter, VerticalSize>(contributor) * at(contributor)->getK();

                            if (transmissivityNeig != 0 * si::square_meter / day and
                                transmissivitySelf != 0 * si::square_meter / day) {
                                t_meter nodeWidth = std::min(getNodeWidth(contributor),
                                                             at(contributor)->getNodeWidth(contributor));
                                // conductance from contributor node to ghost node
                                // if refined into four, multiplierContributor is 0.5, multiplierNodeOuter is 0.25
                                contributorConductance = nodeWidth *
                                                         ((transmissivitySelf * transmissivityNeig)
                                                          / (transmissivitySelf *
                                                             at(contributor)->getNodeLength(contributor) *
                                                             multiplierContributor +
                                                             transmissivityNeig * getNodeLength(contributor) *
                                                             multiplierNodeOuter));
                            }
                            // conductance from this node's center to ghost node inside this node
                            // if refined into four, multiplierNodeInner is 0.25
                            nodeConductance = getNodeWidth(contributor) *
                                              (transmissivitySelf / (getNodeLength(contributor) * multiplierNodeInner));

                            // the alpha coefficient is used to weigh influence on ghost node height difference
                            t_dim alpha = contributorConductance / (contributorConductance + nodeConductance);
                            //LOG(debug) << "alpha = " << alpha << ", contributorConductance = " << contributorConductance.value();

                            // conductance to the refined neighbour node
                            t_s_meter_t conductance =
                                    mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(refinedNeig));

                            // calculate ghost node correction
                            gnc = conductance * (alpha * (getHead() - at(contributor)->getHead()));
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

            static std::forward_list<NeighbourPosition>
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
                    std::forward_list<NeighbourPosition> emptyList{};
                    return emptyList;
                }
            }

            static std::forward_list<NeighbourPosition>
            getPossibleNeighbours_vertical() {
                return {NeighbourPosition::TOP, NeighbourPosition::DOWN};
            }

            std::forward_list<NeighbourPosition>
            getPossibleRefinedNeighbours() {
                if (hasRefinedNeighbour()) {
                    return {NeighbourPosition::FRONTLEFT,  NeighbourPosition::FRONTFRONT, NeighbourPosition::FRONTRIGHT,
                            NeighbourPosition::FRONTFRONTLEFT,  NeighbourPosition::FRONTFRONTRIGHT,
                            NeighbourPosition::BACKLEFT,   NeighbourPosition::BACKBACK,  NeighbourPosition::BACKRIGHT,
                            NeighbourPosition::BACKBACKLEFT,   NeighbourPosition::BACKBACKRIGHT,
                            NeighbourPosition::LEFTFRONT,  NeighbourPosition::LEFTLEFT,  NeighbourPosition::LEFTBACK,
                            NeighbourPosition::LEFTLEFTFRONT,  NeighbourPosition::LEFTLEFTBACK,
                            NeighbourPosition::RIGHTFRONT, NeighbourPosition::RIGHTRIGHT, NeighbourPosition::RIGHTBACK,
                            NeighbourPosition::RIGHTRIGHTFRONT, NeighbourPosition::RIGHTRIGHTBACK};
                } else {
                    std::forward_list<NeighbourPosition> emptyList{};
                    return emptyList;
                }
            }

            static std::forward_list<NeighbourPosition>
            getNeigPos_LRFB(){
                return {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                        NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
            }

            std::forward_list<NeighbourPosition>
            getPossibleNeighbours_horizontal(){
                std::forward_list<NeighbourPosition> possiblePositions = getNeigPos_LRFB();
                if (hasRefinedNeighbour()){
                    possiblePositions.merge(getPossibleRefinedNeighbours());
                }
                return possiblePositions;
            }

            std::forward_list<NeighbourPosition>
            getPossibleNeighbours(){
                std::forward_list<NeighbourPosition> possiblePositions = getPossibleNeighbours_vertical();
                possiblePositions.merge(getPossibleNeighbours_horizontal());
                return possiblePositions;
            }
            /**
             * @brief Add a zeta surface to the cell (bounded by elevation at top and by cell bottom at bottom).
             * @param localZetaID
             * @param zeta the zeta surface height in meters
             */
            void addZeta(int localZetaID, t_meter zeta){
                NANChecker(zeta.value(), "zeta (in addZeta)");
                // limit zeta value to node top and bottom (in SWI2: lines 660-680)
                if (zeta > getElevation() - get<t_meter, VDFLock>()) {
                    zeta = getElevation();
                } else if (zeta < getBottom() + get<t_meter, VDFLock>()) {
                    zeta = getBottom();
                }

                // add initial values to vectors "zetas" and "zetasChange"
                std::vector<t_meter> zetas;
                std::vector<t_meter> zetasChange;
                if (localZetaID == 0) {
                    zetas.push_back(zeta);
                    zetasChange.push_back(0 * si::meter);
                } else {
                    zetas = getZetas();
                    zetasChange = getZetasChange();
                    zetas.insert(zetas.begin() + localZetaID, zeta);
                    zetasChange.insert(zetasChange.begin() + localZetaID, 0 * si::meter);
                }

                // check if zetas vector is sorted
                for (int zetaID = 0; zetaID < zetas.size() - 1; ++zetaID) {
                    if (zetas[zetaID] < zetas[zetaID+1]) {
                        LOG(userinfo) << "At nodeID " << getID() << ": zeta at ID= " << zetaID << ":" <<
                                      zetas[zetaID].value() << ", zeta at ID+1: " <<  zetas[zetaID+1].value();
                        throw "Vector of zetas needs to be sorted!";
                    }
                }

                setZetas(zetas);
                set<std::vector<t_meter>,ZetasChange>(zetas);
                //LOG(debug) << "nodeID: " << getID() << ", localZetaID: " << localZetaID << ", zeta: " << zeta.value() << ", getZeta: " << getZeta(localZetaID).value();
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
                    // limit zeta value to node top and bottom (in SWI2: lines 660-680)
                    auto topOfNode = get<t_meter, Elevation>();
                    if (zeta > topOfNode - get<t_meter, VDFLock>()) {
                        zetas[localZetaID] = topOfNode;
                    } else if (zeta < getBottom() + get<t_meter, VDFLock>()) {
                        zetas[localZetaID] = getBottom();
                    } else {
                        zetas[localZetaID] = zeta;
                    }
                }
                setZetas(zetas);
                //LOG(debug) << "nodeID: " << getID() << ", setZeta[" << localZetaID << "] = " << getZeta(localZetaID).value();
            }

            /**
             * @brief Update zetas after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            void setZetas(std::vector<t_meter> zetas) { set<std::vector<t_meter>, Zetas>(zetas); }

            void initZetasTZero() noexcept { set<std::vector<t_meter>, Zetas_TZero>(get<std::vector<t_meter>, Zetas>()); }

            /**
             * @brief Update zeta change after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            void setZetaChange(int localZetaID, t_meter zeta){
                NANChecker(zeta.value(), "zeta (in setZetaChange)");
                auto zetasChange = getZetasChange();
                if (localZetaID < zetasChange.size()) {
                    zetasChange[localZetaID] = zeta - getZeta(localZetaID);
                }
                set<std::vector<t_meter>, ZetasChange>(zetasChange);
            }

            void updateZetasTZero(){ set < std::vector<t_meter>, Zetas_TZero > (get<std::vector<t_meter>, Zetas>()); }

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

                if (localZetaID < get<std::vector<t_meter>, Zetas>().size()){
                    auto zetas = getZetas();
                    auto zeta = zetas[localZetaID];
                    return zeta;
                } else {
                    throw "Not set at nodeID " + std::to_string(getID()) +
                          ": Zetas[localZetaID = " + std::to_string(localZetaID) + "]";
                }
            }

            /**
             * @brief get one zeta surface height
             * @param localZetaID
             * @return meter
             */
            t_meter getZetaIfActive(int localZetaID) {

                if (localZetaID < get<std::vector<t_meter>, Zetas>().size()){
                    if (isZetaActive(localZetaID)) {
                        return get<std::vector<t_meter>, Zetas>()[localZetaID];
                    } else {
                        return 0 * si::meter;
                    }
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
                if (localZetaID < get<std::vector<t_meter>, ZetasChange>().size()){
                    return get<std::vector<t_meter>, ZetasChange>()[localZetaID];
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

            /**
             *  @brief set the zone(s) of sources and sinks
             *  @param zoneSinks
             *  @param zoneSources
             *  @param numOfZones
             *  @note from SWI2 doc:
             * - ISOURCE > 0: sources and sinks are of the same type as water in zone ISOURCE. If such a zone is not present in the cell, sources and sinks interact with the zone at the top of the aquifer.
             * - ISOURCE = 0: sources and sinks are of the same type as water at top of aquifer.
             * - ISOURCE < 0: for submarine groundwater discharge (SGD)
             *      sources are of the same type as water in zone abs(ISOURCE) (e.g. salt water from ocean)
             *      sinks are the same water type as the water at top of aquifer (e.g. fresh water from "upstream")
             * note: in SWI2 code, ISOURCE is IZONENR
             *
             * Here we use zoneOfSinks and zoneOfSources (containing values between 0 and number of density zones).
             * Thus, sources and sinks are associated to the respective zone. Rule: zoneOfSinks <= zoneOfSources
             * For simulation of submarine groundwater discharge:
             * - zoneOfSinks: 0 (fresh water)
             * - zoneOfSources: an integer between 1 and the number of zones (brackish/saline water)
             */
            void setZoneOfSinksAndSources(int zoneSinks, int zoneSources, int numOfZones) {
                if (zoneSinks > zoneSources) {
                    throw "Zone number of sinks cannot be larger than zone number of sources";
                }

                if (zoneSinks < zoneSources and zoneSinks > 0) {
                    throw "If zone number of sinks and sources are not equal (for simulating SGD), zone number of sinks must be 0";
                }

                if (zoneSinks < numOfZones and zoneSinks >= 0) {
                    zoneOfSinks = zoneSinks;
                } else {
                    throw "Zone number of sinks must be larger or equal to 0, and below the number of density zones";
                }

                if (zoneSources < numOfZones and zoneSources >= 0) {
                    zoneOfSources = zoneSources;
                } else {
                    throw "Zone number of sources must be larger or equal to 0, and below the number of density zones";
                }
            }

            int getZoneOfSources(){ return zoneOfSources;}

            int getZoneOfSinks(){ return zoneOfSinks;}

            /**
             * @brief Set effective porosity (applied to all layers below)
             * @param effectivePorosity effective porosity
             */
            void setEffectivePorosity(t_dim effectivePorosity) {
                setEffectivePorosity_direct(effectivePorosity);
                applyToAllLayers([&effectivePorosity](NodeInterface *nodeInterface) {
                    nodeInterface->getProperties().set<t_dim, EffectivePorosity>(effectivePorosity);
                });
            }

            /**
             * @brief Set effective porosity
             * @param effectivePorosity effective porosity in node
             */
            void setEffectivePorosity_direct(t_dim effectivePorosity) { set<t_dim, EffectivePorosity>(effectivePorosity); }

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
                        getExternalFlowByName(RIVER_MM).setLockRecharge(recharge); //TODO: never used; in calcERC read but not used
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
            * @brief Check if node has any refined neighbour
            * @return bool
            */
            bool hasRefinedNeighbour() { return get<bool, HasRefinedNeighbour>(); }

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
                } catch (const std::out_of_range &e) {//ignore me cell has no special_flow
                }
                t_vol_t eqFlow = getEqFlow();
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                for (const auto &flow : externalFlows) {
                    out += flow.second.getQ(eq_head, head, recharge, eqFlow) * get<t_dim, StepModifier>();
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
                                  (day * get<t_dim, StepModifier>()); // * get<t_dim, StepModifier>()
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
             * @note using ZetasTZero since they are the zetas from the beginning of the current iterative equation solving process
             */
            bool isZetaAtTop(int localZetaID){
                return (getZetaTZero(localZetaID) >= getZetasTZero().front());
            }

            /**
             * @brief Checking if zeta surface is at node bottom
             * @param localZetaID
             * @return bool
             * @note using ZetasTZero since they are the zetas from the beginning of the current iterative equation solving process
             */
            bool isZetaAtBottom(int localZetaID){
                return (getZetaTZero(localZetaID) <= getZetasTZero().back() or
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
             * @brief The source flow below a zeta surface (for the right hand side in the zeta equation)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note like G in SWI2 doc, but without vertical leakage; in SWI2 code: lines 3523-3569
             * G = RHS (of flow, for constant density) - HCOF_(i,j,k,n)*h^(m)_(i,j,k) + (verticalLeakage_(i,j,k-1,n) - verticalLeakage_(i,j,k,n))
             */
            t_vol_t getSources(int localZetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                /* if nodes can be inactive: return 0 at inactive nodes
                if (nodeInactive) { return out; }
                 */
                // todo: compute BUFF with SSWI2_BDCH for constant head cells
                t_vol_t sources = 0.0 * (si::cubic_meter / day);
                //LOG(debug) << "zoneToUse: " << zoneToUse << std::endl;

                if (isZetaActive(localZetaID)) { // if "iz.NE.1" and IPLPOS == 0 (line 3570-3571)
                    // if the new groundwater head is above or equal to the node bottom
                    if (get<t_meter, Head>() >= getBottom()) { // lines 3532-3536
                        // get RHS of flow equation without VDF terms (pseudo source term and flux correction)
                        auto RHSConstantDensity = get<t_vol_t, RHSConstantDensity_TZero>(); // in SWI2 code: RHSPRESWI

                        // get HCOF (= P-((SS_i,j,k * DELR_j * DELC_i)/(t^m - t^m-1))) of NEW time step
                        t_s_meter_t hcof = mechanics.getHCOF(steadyState,
                                                             get<t_dim, StepModifier>(),
                                                             getStorageCapacity(),
                                                             getP());
                        // calculate the boundary flux
                        sources = RHSConstantDensity - hcof * get<t_meter, Head>(); // see line 3535-3536
                        //LOG(debug) << "RHSConstantDensity: " << RHSConstantDensity.value() << std::endl;
                        //LOG(debug) << "hcof: " << hcof.value() << std::endl;
                        //LOG(debug) << "boundaryFlux: " << boundaryFlux.value() << std::endl;

                    }

                    int zoneToUse = zoneOfSources; // zone of flow sources in node
                    if (zoneOfSources > zoneOfSinks and // if we intend to simulate submarine groundwater discharge
                        sources > 0 * (si::cubic_meter / day)) { // and boundary flux is positive (-> out of node)
                        zoneToUse = zoneOfSinks; // use the zone of sinks (= top of aquifer = fresh water)
                    }
                    // Question: add 3539-3540 with "if (ibound < 0){q=-BUFF}"?
                    t_dim factor = 1 * si::si_dimensionless;
                    if (localZetaID <= zoneToUse and sources != 0 * (si::cubic_meter / day)) {
                        if (zoneToUse > 100){ // lines 3548-3553
                            if (sources > 0 * (si::cubic_meter / day) and
                                (getZetas().front() - getZetas().back()) > 0 * si::meter) {
                                factor = (getZeta(localZetaID) - getZetas().back()) / (getZetas().front() - getZetas().back());
                            }
                        }
                        out = sources * factor;
                    }
                }

                NANChecker(out.value(), "getSources");
                return out;
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
                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();

                // pseudo source term calculation (in 2 parts)
                for (const auto &position: possible_neighbours) {
                    auto got = neighbours.find(position);
                    if (got != neighbours.end()) {
                        // calculating zone conductances for pseudo source term calculation
                        std::vector<t_s_meter_t> zoneConductances = getZoneConductances(got);

                        if ((isZetaActive(localZetaID) and // if "iz.NE.1" and IPLPOS == 0 (line 3570-3571)
                             at(got)->isZetaActive(localZetaID))) { // if neighbouring IPLPOS == 0 (e.g. line 2094)
                            //%% head part %%
                            t_s_meter_t zoneCondCumHead = getZoneConductanceCum(localZetaID,zoneConductances);
                            t_vol_t head_part = -zoneCondCumHead * (getAt<t_meter, Head>(got) - get<t_meter, Head>());
                            out += head_part;
                            //LOG(debug) << "head_part: " << head_part.value() << std::endl;

                            t_s_meter_t zoneCondCumDelnus;
                            for (int zetaID = 0; zetaID < getZetas().size() - 1; zetaID++) {
                                //%% delnus part %%
                                if (zetaID < localZetaID) {
                                    zoneCondCumDelnus = getZoneConductanceCum(localZetaID, zoneConductances);
                                } else if (zetaID == localZetaID) {
                                    zoneCondCumDelnus = 0;
                                } else if (zetaID > localZetaID) {
                                    zoneCondCumDelnus = getZoneConductanceCum(zetaID, zoneConductances);
                                }
                                t_vol_t delnus_part = -delnus[zetaID] * zoneCondCumDelnus *
                                                      (at(got)->getZeta(zetaID) - getZeta(zetaID));
                                out += delnus_part;
                                //LOG(debug) << "delnus_part (zetaID = " << zetaID << "): " << delnus_part.value() << std::endl;
                            }
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
            t_vol_t getFluxHorizontal(map_itter got, int localZetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<std::vector<t_dim>, Delnus>();

                std::vector<t_s_meter_t> zoneConductances = getZoneConductances(got);
                t_s_meter_t zoneCondCum = getZoneConductanceCum(localZetaID, zoneConductances);
                // %%head part %% for left/back neighbour
                t_vol_t head_part = zoneCondCum * (getAt<t_meter, Head>(got) - get<t_meter, Head>());
                out += head_part;
                //LOG(debug) << "head_part (tip/toe): " << head_part.value() << std::endl;

                // %%delnus part %% for left/back neighbour
                t_s_meter_t zoneCondCumZeta;
                for (int zetaID = 0; zetaID < getZetas().size() - 1; zetaID++) {
                    if (zetaID <= localZetaID){
                        zoneCondCumZeta = getZoneConductanceCum(localZetaID,zoneConductances);
                    } else {
                        zoneCondCumZeta = getZoneConductanceCum(zetaID,zoneConductances);
                    }
                    t_vol_t delnus_part = delnus[zetaID] * zoneCondCumZeta *
                                          (at(got)->getZeta(zetaID) - getZeta(zetaID));
                    out += delnus_part;
                    //LOG(debug) << "delnus_part: " << delnus_part.value() << std::endl;
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

                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours();

                // tip and toe flow calculation (in 8 parts)
                if (isZetaActive(localZetaID)) { // if "iz.NE.1" and IPLPOS == 0 (line 3570-3571)
                    for (const auto &position: possible_neighbours) {
                        auto got = neighbours.find(position);
                        if (got != neighbours.end()) {
                            if (!at(got)->isZetaActive(localZetaID)) { // if "iz.NE.1" and IPLPOS == 0 (line 3570-3571)
                                if ((position == NeighbourPosition::LEFT) or
                                    (position == NeighbourPosition::BACK) or
                                    (position == NeighbourPosition::RIGHT) or
                                    (position == NeighbourPosition::FRONT)) {
                                    out -= getFluxHorizontal(got, localZetaID); // SSWI2_QR (left/right), SSWI2_QC (front/back)
                                } else if (position == NeighbourPosition::TOP) {
                                    // vertical leakage to TOP neighbour
                                    auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                                    if (getFluxTop() < 0 * (si::cubic_meter / day) and
                                        at(got)->getNusBot() >= nusInZones[localZetaID] and
                                        getNusBot() >= at(got)->getNusBot()) { // IF ((qztop.LT.0).AND.(NUBOT(i,j,k-1).GE.NUS(iz)).AND.(NUBOT(j,i,k).GE.NUBOT(i,j,k-1))) THEN
                                        out += getFluxTop(); // in SWI2: qztop
                                    }
                                } else if (position == NeighbourPosition::DOWN) {
                                    // vertical leakage to DOWN neighbour
                                    auto nusInZones = get<std::vector<t_dim>, NusInZones>();
                                    if(getFluxDown() < 0 * (si::cubic_meter / day) and
                                       at(got)->getNusTop() < nusInZones[localZetaID] and
                                       getNusTop() <= at(got)->getNusTop()) { // IF ((qzbot.LT.0).AND.(NUTOP(i,j,k+1).LT.NUS(iz)).AND.(NUTOP(j,i,k).LE.NUTOP(i,j,k+1))) THEN
                                        continue;
                                    } else{
                                        out += getFluxDown(); // in SWI2: qzbot
                                    }
                                }
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
            std::vector<t_s_meter_t> getZoneConductances(map_itter got) {
                std::vector<t_s_meter_t> out;
                t_s_meter_t conductance;
                t_s_meter_t zoneConductance;
                std::vector<t_meter> zoneThicknesses;
                t_meter zoneThickness;
                t_meter sumOfZoneThicknesses = 0 * si::meter;
                t_meter deltaZeta;
                t_meter deltaZeta_neig;

                // calculate the zone thicknesses and their sum
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; localZetaID++) {
                    deltaZeta = getZeta(localZetaID) - getZeta(localZetaID + 1);
                    deltaZeta_neig = at(got)->getZeta(localZetaID) - at(got)->getZeta(localZetaID + 1);
                    if (deltaZeta <= (0 * si::meter) or deltaZeta_neig <= (0 * si::meter)){ // adapted from SWI2 code line 1149
                        zoneThickness = 0 * si::meter;
                    } else {
                        zoneThickness = ((getLengthNeig(got) * deltaZeta) + (getNodeLength(got) * deltaZeta_neig)) /
                                        (getLengthNeig(got) + getNodeLength(got));
                        sumOfZoneThicknesses += zoneThickness;
                    }
                    NANChecker(zoneThickness.value(), "zoneThickness");
                    zoneThicknesses.push_back(zoneThickness);
                }

                // calculate the density zone conductances
                for (int localZetaID = 0; localZetaID < getZetas().size() - 1; localZetaID++) {
                    //LOG(debug) << "zoneThicknesses[" << localZetaID << "]:" << zoneThicknesses[localZetaID].value();
                    zoneConductance = 0 * si::square_meter / day;
                    if (sumOfZoneThicknesses != (0 * si::meter)) { // adapted from SWI2 code line 1159
                        conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got));
                        //LOG(debug) << "conductance:" << conductance.value();
                        zoneConductance = conductance * (zoneThicknesses[localZetaID] / sumOfZoneThicknesses);
                        //LOG(debug) << "zoneConductance[" << localZetaID << "] :" << zoneConductance.value();
                    }
                    /* if nodes can be inactive:
                    else {
                    if (nodeInactive()){
                        if (1 < localZetaID < (getZetas().size() - 2) {
                            if((getZetaPosInNode(localZetaID) != "between") or // adapted from SWI2 code lines 1212-1222
                               (at(got)->getZetaPosInNode(localZetaID) != "between") or
                               (getZetaPosInNode(localZetaID + 1) != "between") or
                               (at(got)->getZetaPosInNode(localZetaID + 1) != "between"))) {
                                // this section is reached if getZetas().size() >= 4 and numberOfZones >= 3
                                zoneConductance = 0 * si::square_meter / day; // zone is inactive or across layers
                            }
                        }
                    } */
                    //LOG(debug) << "zoneConductance[id = " << localZetaID << "] : " << zoneConductance.value() << std::endl;
                    out.push_back(zoneConductance);
                    NANChecker(zoneConductance.value(), "zoneConductance");
                }
                return out;
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
                /* if nodes can be inactive: return 0 at inactive nodes
                if (nodeInactive) { return out; }
                 */
                auto delnus = get<std::vector<t_dim>, Delnus>();
                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();

                // pseudo source term calculation
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got != neighbours.end()) {
                        if (isAnyZetaBetween() or
                            at(got)->isAnyZetaBetween()) { // check if there are any active zeta surfaces
                            std::vector<t_s_meter_t> zoneConductances = getZoneConductances(got);
                            for (int zetaID = 0; zetaID < getZetas().size() - 1; zetaID++) {
                                t_s_meter_t zoneConductanceCum = getZoneConductanceCum(zetaID, zoneConductances);
                                t_vol_t pseudoSource = delnus[zetaID] * zoneConductanceCum *
                                                       (at(got)->getZeta(zetaID) - getZeta(zetaID));
                                out -= pseudoSource;
                            }
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

                auto got = neighbours.find(NeighbourPosition::TOP);
                if (got != neighbours.end()){ //Current node has a top node
                    // first part of the flux correction term
                    for (int localZetaID = 0; localZetaID < getZetas().size() - 1; localZetaID++){
                        headdiff -= nusInZones[localZetaID] *
                                    (at(got)->getZeta(localZetaID + 1) - at(got)->getZeta(localZetaID));
                        // Note: in SWI2 documentation is, BOUY is calculated by adding headdiff (would be out +=),
                        // MODFLOW code for headdiff is as implemented here (with out -=)
                    }
                    // second part of the flux correction term (vertical conductance *  BOUY)
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                    out = verticalConductance *
                          (headdiff +
                           0.5 * (at(got)->getZetas().back() - getZetas().front()) *
                           (at(got)->getNusBot() + getNusTop()));
                    // Note in SWI2 documentation, BOUY is calculated with a - between NUBOT and NUTOP,
                    // in MODFLOW code there is a + in the calculation of QLEXTRA
                    //LOG(debug) << "headdiff: " << headdiff.value() << std::endl;
                }
                return out;
            }

            t_vol_t getVerticalFluxCorrections(){
                t_vol_t out = 0 * (si::cubic_meter / day);

                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_vertical();

                for (const auto &position : possible_neighbours) {
                    auto got = neighbours.find(position); // todo: enhance (this check is done multiple times for top)
                    if (got != neighbours.end()) {//Current node has a top or down node
                        if (position == NeighbourPosition::TOP) {
                            out -= getVerticalFluxCorrection();
                        }
                        if (position == NeighbourPosition::DOWN) {
                            out += at(got)->getVerticalFluxCorrection();
                        }
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
                auto got = neighbours.find(NeighbourPosition::TOP);
                if (got != neighbours.end()) {
                    t_vol_t fluxFromTopNode = getVerticalFluxCorrection();
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                    out = (verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(got))) - fluxFromTopNode;
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
                auto got = neighbours.find(NeighbourPosition::DOWN);
                if (got != neighbours.end()) {
                    t_vol_t fluxFromDownNode = at(got)->getVerticalFluxCorrection();
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                    out = (verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(got))) + fluxFromDownNode;
                    //LOG(debug) << "getFluxDown: " << out.value();
                }
                return out;
            }

            /**
             * @brief Updating top zeta surface height after the flow/head solution is found
             * @note Top zeta surface height is set to the new groundwater head, only if unconfined. In SWI2: SSWI2_UPZ1
             */
            void setTopZetaToHead(){
                if (get<bool, Confinement>()){
                    return;
                } else { // only clip if node is unconfined
                    auto head = get<t_meter, Head>();
                    t_meter bottomOfNode = getBottom();
                    auto topOfNode = getElevation();
                    t_meter newHeight;
                    // if groundwater head is BELOW the top of the node
                    if (head < topOfNode) {
                        // if groundwater head is ABOVE the bottom of the node
                        if (head > bottomOfNode) {
                            newHeight = head;
                            // if groundwater head is BELOW OR EQUAL to the bottom of the node
                        } else { // head <= bottomOfNode
                            newHeight = bottomOfNode;
                        }
                        // update the first zeta surface
                        setZeta(0, newHeight);

                        // update all other zeta surfaces that are ABOVE the updated first zeta surface
                        for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                            if (getZeta(localZetaID) > newHeight) {
                                setZeta(localZetaID, newHeight);
                            }
                        }
                        // if groundwater head is ABOVE OR EQUAL to the top of the node
                    } else { // head >= topOfNode
                        // clip zeta to the top of the node
                        setZeta(0, topOfNode);
                    }
                }
            }

            /**
             * @brief Vertical movement of zeta surfaces through top of this node. This function is required to move a
             * zeta surface with the height equal to the top of a lower node and bottom to an upper node. Without this
             * function, that zeta surface would be stuck there since the other only consider "active" surfaces (which
             * are between the top and bottom of a node)
             * @note in SWI2 code: SSWI2_VERTMOVE
             */
            void verticalZetaMovement() {
                // skip nodes where head is below the bottom of the node
                t_vol_t fluxCorrectionTop; // in SWI2: qztop
                t_s_meter_t verticalConductanceTop;
                t_meter deltaZeta;
                auto top = neighbours.find(NeighbourPosition::TOP);
                // skip nodes that do not have a top neighbour
                if (top == neighbours.end()) { // no neighbour at position
                } else {
                    // if head is above the bottom of the node, both in this node AND in the node above (SWI2 line 2325)
                    if (get<t_meter, Head>() >= getBottom() and
                        getAt<t_meter, Head>(top) >= (getAt<t_meter, Elevation>(top) - getAt<t_meter, VerticalSize>(top))) {

                        for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                            // zeta only moves through the top of a node if there is a ZETA surface
                            // - at the top of the current node (in SWI2: IPLPOS_(i,j,k,n) = 1)
                            // - AND at the bottom of the top node (in SWI2: IPLPOS_(i,j,k-1,n) = 2)
                            if (isZetaAtTop(localZetaID) and at(top)->isZetaAtBottom(localZetaID)) {

                                // calculate flux through the top
                                fluxCorrectionTop = getFluxTop();
                                //LOG(debug) << "fluxCorrectionTop: " << fluxCorrectionTop.value() << std::endl;

                                // if vertical flux through the top of the node is positive...
                                if (fluxCorrectionTop.value() > 0 and at(top)->getEffectivePorosity().value() > 0) {
                                    deltaZeta = (fluxCorrectionTop * (day)) /
                                                (get<t_s_meter, Area>() * getAt<t_dim, EffectivePorosity>(top) * get<t_dim, StepModifier>()); // * get<t_dim, StepModifier>()
                                    // ...lift zeta height of the lowest zeta surface in top node
                                    t_meter zeta_back_top = at(top)->getZetas().back();
                                    at(top)->setZeta(localZetaID, zeta_back_top + deltaZeta);

                                    // if vertical flux through the top of the node is negative...
                                } else if (fluxCorrectionTop.value() < 0 and getEffectivePorosity().value() > 0) {
                                    deltaZeta = (fluxCorrectionTop * (day)) /
                                                (get<t_s_meter, Area>() * getEffectivePorosity() * get<t_dim, StepModifier>()); // * get<t_dim, StepModifier>()
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
                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();

                for (const auto &position : possible_neighbours) {
                    auto got = neighbours.find(position);
                    auto got_opp = neighbours.find(getOppositePosition(position));
                    if (got != neighbours.end()) {
                        for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                            if (isZetaActive(localZetaID)) {
                                //if () // todo

                                // get max delta of zeta between nodes
                                if (at(got)->isZetaAtBottom(localZetaID)) {
                                    maxDelta = 0.5 * (getNodeLength(got) + getLengthNeig(got)) * get<t_dim, MaxToeSlope>();
                                } else if (at(got)->isZetaAtTop(localZetaID)) {
                                    maxDelta = 0.5 * (getNodeLength(got) + getLengthNeig(got)) * get<t_dim, MaxTipSlope>();
                                }
                                //LOG(debug) << "maxDelta: " << maxDelta.value() << std::endl;


                                // if tracking tip/toe: raise/lower this zeta surface in this node by:
                                delta_self = get<t_dim, SlopeAdjFactor>() * maxDelta *
                                             ((at(got)->getEffectivePorosity() * getLengthNeig(got)) /
                                              ((getEffectivePorosity() * getNodeLength(got)) +
                                               (at(got)->getEffectivePorosity() * getLengthNeig(got))));
                                // if tracking tip/toe: lower/raise this zeta surface in neighbouring node by:
                                delta_neig = get<t_dim, SlopeAdjFactor>() * maxDelta *
                                             ((getEffectivePorosity() * getNodeLength(got)) /
                                              ((getEffectivePorosity() * getNodeLength(got)) +
                                               (at(got)->getEffectivePorosity() * getLengthNeig(got))));

                                if (at(got)->isZetaAtBottom(localZetaID)) {
                                    //%% Toe tracking %%
                                    t_meter zetaDif = getZeta(localZetaID) - at(got)->getZetas().back();
                                    //LOG(debug) << "zetaDif (toe): " << zetaDif.value() << std::endl;
                                    if (zetaDif > maxDelta) {
                                        setZeta(localZetaID, getZeta(localZetaID) - delta_self);
                                        //LOG(debug) << "delta_self (toe): " << delta_self.value() << std::endl;
                                        t_meter zeta_back_neig = at(got)->getZetas().back();
                                        at(got)->setZeta(localZetaID, zeta_back_neig + delta_neig);
                                        //LOG(debug) << "delta_neig (toe): " << delta_neig.value() << std::endl;
                                    }
                                } else if (at(got)->isZetaAtTop(localZetaID)) {
                                    //%% Tip tracking %%
                                    t_meter zetaDif = at(got)->getZetas().front() - getZeta(localZetaID);
                                    //LOG(debug) << "zetaDif (tip): " << zetaDif.value() << std::endl;
                                    if (zetaDif > maxDelta) {
                                        setZeta(localZetaID, getZeta(localZetaID) + delta_self);
                                        //LOG(debug) << "zetaChange_self (tip): " << delta_self.value() << std::endl;
                                        t_meter zeta_front_neig = at(got)->getZetas().front();
                                        at(got)->setZeta(localZetaID, zeta_front_neig - delta_neig);
                                        //LOG(debug) << "zetaChange_neig (tip): " << delta_neig.value() << std::endl;
                                    }
                                }

                                if ((getZeta(localZetaID) - getZetas().back()) < (get<t_dim, MinDepthFactor>() * delta_neig)) {
                                    if (got_opp != neighbours.end()){
                                        if (at(got_opp)->isZetaActive(localZetaID)) {
                                            // change zeta in other direction neighbour
                                            delta_opp = ((getZeta(localZetaID) - getZetas().back()) *
                                                         (getNodeLength(got) * getEffectivePorosity()) /
                                                         (getLengthNeig(got_opp) *
                                                          at(got_opp)->getEffectivePorosity()));
                                            t_meter zeta_opp = at(got_opp)->getZeta(localZetaID);
                                            at(got_opp)->setZeta(localZetaID, zeta_opp + delta_opp);
                                            setZeta(localZetaID, getZetas().back());
                                            //LOG(debug) << "delta_opp (toe): " << delta_opp.value() << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            NeighbourPosition getOppositePosition(NeighbourPosition position) {
                if (position == NeighbourPosition::BACK or
                    position == NeighbourPosition::BACKLEFT or position == NeighbourPosition::BACKRIGHT) {
                    return NeighbourPosition::FRONT;
                } else if (position == NeighbourPosition::FRONT or
                           position == NeighbourPosition::FRONTLEFT or position == NeighbourPosition::FRONTRIGHT) {
                    return NeighbourPosition::BACK;
                } else if (position == NeighbourPosition::LEFT or
                           position == NeighbourPosition::LEFTFRONT or position == NeighbourPosition::LEFTBACK) {
                    return NeighbourPosition::RIGHT;
                } else if (position == NeighbourPosition::RIGHT or
                           position == NeighbourPosition::RIGHTFRONT or position == NeighbourPosition::RIGHTBACK) {
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
                /* if nodes can be inactive: return at inactive nodes
                if (nodeInactive) { return; }
                 */
                t_meter maxDelta;
                t_meter maxAdjustment_self;
                t_meter maxAdjustment_neig;
                t_meter delta_self;
                t_meter delta_neig;
                t_meter zeta_neig;
                for (int localZetaID = 1; localZetaID < getZetas().size() - 1; localZetaID++) {
                    if (isZetaAtTop(localZetaID) or isZetaAtBottom(localZetaID)) {
                        // iterate through horizontal neighbours
                        std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();
                        for (const auto &position: possible_neighbours) {
                            auto got = neighbours.find(position);
                            if (got != neighbours.end()) {
                                /* if nodes can be inactive: return at inactive nodes
                                if (at(got)->nodeInactive) { return; }
                                 */

                                // determine max delta zeta
                                maxAdjustment_self = get<t_meter, VerticalSize>() * get<t_dim, SlopeAdjFactor>();
                                maxAdjustment_neig = getAt<t_meter, VerticalSize>(got) * get<t_dim, SlopeAdjFactor>();
                                if (get<t_meter, VDFLock>() > maxAdjustment_self or
                                    get<t_meter, VDFLock>() > maxAdjustment_neig) {
                                    maxDelta = std::min(maxAdjustment_self, maxAdjustment_neig);
                                } else {
                                    maxDelta = get<t_meter, VDFLock>();
                                }

                                if (getEffectivePorosity().value() == 0 and
                                    at(got)->getEffectivePorosity().value() == 0){
                                    continue;
                                }
                                // calculate delta zeta of node and neighbour
                                delta_self = maxDelta * (at(got)->getEffectivePorosity() * getLengthNeig(got)) /
                                             (getEffectivePorosity() * getNodeLength(got) +
                                              at(got)->getEffectivePorosity() * getLengthNeig(got));
                                delta_neig = maxDelta * (getEffectivePorosity() * getNodeLength(got)) /
                                             (getEffectivePorosity() * getNodeLength(got) +
                                              at(got)->getEffectivePorosity() * getLengthNeig(got));

                                // if a zeta surface is at the BOTTOM of this node and at the TOP of the neighbour
                                // else if a zeta surface is at the TOP of this node and at the BOTTOM of the neighbour
                                if (isZetaAtBottom(localZetaID) && // IPLPOS_self = 2
                                    at(got)->isZetaAtTop(localZetaID)) { // IPLPOS_neig = 1
                                    // adjust zeta surface heights
                                    setZeta(localZetaID, getZeta(localZetaID) + delta_self);
                                    zeta_neig = at(got)->getZeta(localZetaID);
                                    at(got)->setZeta(localZetaID, zeta_neig - delta_neig);
                                } else if (isZetaAtTop(localZetaID) && // IPLPOS_self = 1
                                           at(got)->isZetaAtBottom(localZetaID)) {// IPLPOS_neig = 2
                                    // adjust zeta surface heights
                                    setZeta(localZetaID, getZeta(localZetaID) - delta_self);
                                    zeta_neig = at(got)->getZeta(localZetaID);
                                    at(got)->setZeta(localZetaID, zeta_neig + delta_neig);
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
            t_meter getLengthNeig(map_itter got){
                if ( isLeftOrRight(got->first) ) {
                    return getAt<t_meter, EdgeLengthLeftRight>(got);
                } else if ( isFrontOrBack(got->first) ) {
                    return getAt<t_meter, EdgeLengthFrontBack>(got);
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            /**
             * @brief The length of this node, in direction to the neighbouring node
             * @param got neighbour node
             * @return meter
             */
            t_meter getNodeLength(map_itter got){
                if ( isLeftOrRight(got->first) ){
                    return get<t_meter, EdgeLengthLeftRight>();
                } else if ( isFrontOrBack(got->first) ){
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
            t_meter getNodeWidth(map_itter got) {
                if ( isLeftOrRight(got->first) ){
                    return get<t_meter, EdgeLengthFrontBack>();
                } else if ( isFrontOrBack(got->first)  ){
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
                        if (flow.second.flowIsHeadDependent(head)) {
                            out += flow.second.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepModifier>();
                        }
                    } else { // RECHARGE, ...
                        out += flow.second.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepModifier>();
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
                    if (is(flow.second.getType()).in(RIVER, RIVER_MM, LAKE, GLOBAL_LAKE, WETLAND, GLOBAL_WETLAND, DRAIN)) { // Question: GLOBAL_LAKE?
                        if (not flow.second.flowIsHeadDependent(head)) {
                            out += flow.second.getP(eq_head, head, recharge, eqFlow) * get<t_dim, StepModifier>() *
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
                size_t numC = getNumofNeighbours() + 1; // matrix needs 1 entry per neighbour + 1 entry for this node
                //LOG(debug) << "numC: " << numC;
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(numC);
                t_s_meter_t conduct;

                //Get all conductances from neighbouring cells
                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours();

                for (const auto &position: possible_neighbours) {
                    auto got = neighbours.find(position);
                    conduct = 0 * si::square_meter / day;
                    if (got != neighbours.end()) { //There is a neighbour node
                        if (got->first == TOP or got->first == DOWN) {
                            conduct = mechanics.calculateVerticalConductance(createDataTuple(got));
                            //LOG(debug) << "vertical conductance: " << conduct.value();
                        } else {
                            if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                                conduct = mechanics.calculateEFoldingConductance(createDataTuple<Head>(got),
                                                                                 get<t_meter, EFolding>(),
                                                                                 getAt<t_meter, EFolding>(got));
                                //LOG(debug) << "horizontal conductance (using e-folding): " << conduct.value();
                            } else {
                                conduct = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got));
                                //LOG(debug) << "horizontal conductance: " << conduct.value();
                            }
                        }
                        NANChecker(conduct.value(), "Conductances");

                        // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                        // (used to solve for the head at the neighbouring node)
                        out[nodes->at(got->second)->get<large_num, ID>()] = conduct;
                    }
                }

                // To solve for the head at this node, the conductances to neighbours and HCOF are used
                t_s_meter_t conductNode = 0 * si::square_meter / day;
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &c : out) { conductNode = conductNode - c.second; }
                // add HCOF
                t_s_meter_t hcof = mechanics.getHCOF(steadyState,
                                                     get<t_dim, StepModifier>(),
                                                     getStorageCapacity(),
                                                     getP());
                conductNode = conductNode + hcof;
                //LOG(debug) << "conductNode: " << conductNode.value();

                // check for nan
                NANChecker(conductNode.value(), "conductNode");

                // add resulting conductance to solve for the head at this node to out
                out[get<large_num, ID>()] = conductNode;
                return out;
            };

            /**
             * @brief The matrix entry for the left hand side of the zeta surface equation
             * @param localZetaID zeta surface id in this node
             * @return map <CellID,Conductance>
             */
            std::unordered_map<large_num, t_s_meter_t> getMatrixEntries(int localZetaID) {
                size_t numC = 5;
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(numC);
                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();
                auto delnus = get<std::vector<t_dim>, Delnus>();

                std::vector<t_s_meter_t> zoneConductances;
                t_s_meter_t zoneConductanceCum;
                t_s_meter_t zetaMovementConductance;

                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    zetaMovementConductance = 0 * (si::square_meter / day);
                    if (got != neighbours.end()) { // there is a neighbour at position
                        if ((isZetaActive(localZetaID) and
                             at(got)->isZetaActive(localZetaID))) {
                            zoneConductances = getZoneConductances(got);
                            zoneConductanceCum = getZoneConductanceCum(localZetaID, zoneConductances);
                            zetaMovementConductance += delnus[localZetaID] * zoneConductanceCum; // in SWI2: SWISOLCC/R
                            //LOG(debug) << "zoneConductanceCum = " << zoneConductanceCum.value() << std::endl;
                        }
                        NANChecker(zetaMovementConductance.value(), "zetaMovementConductance");
                        // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                        out[nodes->at(got->second)->get<large_num, ID>()] = zetaMovementConductance;
                    }
                }

                // To solve for zeta in this node, the conductances to neighbours and porosity term are used
                t_s_meter_t conductNode = 0 * (si::square_meter / day); // SWI_HCOF
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &c: out) { conductNode = conductNode - c.second; }
                // subtract effective porosity term
                conductNode = conductNode - getEffectivePorosityTerm(); // subtracting SWIHCOF
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
                t_vol_t extFlows = -getQ(); // e.g., recharge, river, lakes, wetlands
                //LOG(debug) << "extFlows: " << extFlows.value() << std::endl;
                t_vol_t dewateredFlow = calculateDewateredFlow(); // only if node has bottom neighbour
                //LOG(debug) << "dewateredFlow: " << dewateredFlow.value() << std::endl;
                t_vol_t notHeadDependentFlows = calculateNotHeadDependentFlows(); // e.g., rivers, lakes, wetlands //
                //LOG(debug) << "notHeadDependentFlows: " << notHeadDependentFlows.value() << std::endl;
                t_vol_t storageFlow =
                        getStorageCapacity() * (get<t_meter, Head_TZero>() / (day * get<t_dim, StepModifier>()));
                if (steadyState) {
                    storageFlow = 0 * (si::cubic_meter / day);
                }
                //LOG(debug) << "storageFlow: " << storageFlow.value() << std::endl;

                t_vol_t gncFromNodes = 0 * (si::cubic_meter / day);
                //gncFromNodes = getGNCFromNodes();
                //LOG(debug) << "gncFromNodes: " << gncFromNodes.value() << std::endl;

                t_vol_t gncToRefined = 0 * (si::cubic_meter / day);
                //gncToRefined = getGNCToRefinedNode();
                //LOG(debug) << "gncToRefined: " << gncToRefined.value() << std::endl;

                t_vol_t out = extFlows + dewateredFlow - notHeadDependentFlows - storageFlow + gncFromNodes +
                              gncToRefined;
                //LOG(debug) << "RHS constant density: " << out.value() << std::endl;
                NANChecker(out.value(), "RHS constant density");

                if (get<bool, DensityVariable>()) {
                    // save constant density RHS (without variable density terms) for calculation of zeta movement
                    set<t_vol_t, RHSConstantDensity_TZero>(out);

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
                t_vol_t porosityTerm = 0 * (si::cubic_meter / day);
                if (isZetaActive(localZetaID)) { // if "iz.NE.1" and IPLPOS == 0 (line 3570-3571)
                    porosityTerm = getEffectivePorosityTerm() * getZeta(localZetaID);
                }
                //LOG(debug) << "porosityTerm: " << porosityTerm.value() << std::endl;
                t_vol_t sources = getSources(localZetaID); // in SWI2 code: part of BRHS; in SWI2 doc: G or known source term below zeta
                //LOG(debug) << "sources: " << sources.value() << std::endl;
                t_vol_t tipToeFlow = getTipToeFlow(localZetaID); // in SWI2 code: SSWI2_QR and SSWI2_QC
                //LOG(debug) << "tipToeFlow: " << tipToeFlow.value() << std::endl;
                t_vol_t pseudoSourceBelowZeta = getPseudoSourceBelowZeta(localZetaID); // in SWI2 code: SSWI2_SD and SSWI2_SR
                //LOG(debug) << "pseudoSourceBelowZeta: " << pseudoSourceBelowZeta.value() << std::endl;
                t_vol_t out = - porosityTerm - sources + tipToeFlow + pseudoSourceBelowZeta;
                //if (out.value() == 0) {
                //    LOG(debug) << "getRHS(localZetaID) is 0 at " << get<large_num,ID>();
                //}
                NANChecker(out.value(), "getRHS(int localZetaID)");
                return out;
            }

            void setHeadChange(t_meter change) noexcept {
                __setHeadChange(change);
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
            quantity<Velocity> getVelocity(map_itter pos) {
                t_vol_t lateral_flow{0 * si::cubic_meter / day};
                auto vert_size = get<t_meter, VerticalSize>();
                t_s_meter_t conductance;
                if (get<int, Layer>() > 0 and get<bool, UseEfolding>()) {
                    conductance = mechanics.calculateEFoldingConductance(createDataTuple<Head>(pos), get<t_meter, EFolding>(), getAt<t_meter, EFolding>(pos));
                } else {
                    conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(pos));
                }

                lateral_flow = conductance * (get<t_meter, Head>() - getAt<t_meter, Head>(pos));
                return lateral_flow / (vert_size * vert_size);
            }

            /**
             * @brief Calculate flow velocity for flow tracking
             * Vx and Vy represent the flow velocity in x and y direction.
             * A negative value represents a flow in the opposite direction.
             * @return Velocity vector (x,y)
             */
            std::pair<double, double> getVelocityVector() {
                quantity<Velocity> Vx{0};
                quantity<Velocity> Vy{0};

                std::forward_list<NeighbourPosition> possible_neighbours = getPossibleNeighbours_horizontal();

                for (const auto &position: possible_neighbours) {
                    auto got = neighbours.find(position);
                    if (got != neighbours.end()) {
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
                         bool hasRefinedNeighbour,
                         bool densityVariable,
                         std::vector<t_dim> delnus,
                         std::vector<t_dim> nusInZones,
                         double effPorosity,
                         double maxTipSlope,
                         double maxToeSlope,
                         double minDepthFactor,
                         double slopeAdjFactor,
                         t_meter vdfLock)
                    : NodeInterface(nodes, lat, lon, area, edgeLengthLeftRight, edgeLengthFrontBack, SpatID, ID, K,
                                    head, aquiferDepth, anisotropy, specificYield, specificStorage, useEfolding,
                                    confined, refID, hasRefinedNeighbour, densityVariable, delnus, nusInZones,
                                    effPorosity, maxTipSlope, maxToeSlope, minDepthFactor, slopeAdjFactor, vdfLock) {}
        private:
            // implementation
            friend class NodeInterface;

            //Learning weight
            t_dim weight = 0.1;

            /**
             * @brief Update heads after one or multiple inner iterations
             * @param delta
             */
            virtual void __setHeadChange(t_meter change) {
                NANChecker(change.value(), "Set Head Change");
                t_meter current_head = get<t_meter, Head>();
                set<t_meter, HeadChange>(change);
                set<t_meter, Head>(current_head + change);
            };

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
                    0.000015, false, true, 0, false, false,{0.0, 0.1}, {0.0, 0.1},
                    0.2, 0.2, 0.2, 0.1, 0.1,0.001 * si::meter) {}

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
                ar & steadyState;
                ar & fields;
                ar & cached;
                ar & eq_flow;
            }

        };

    }
}

#endif //NODE_HPP
