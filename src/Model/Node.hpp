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
 * FRONT (larger ID) * BACK (smaller ID)
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
        //typedef std::hash<underlying_type>::result_type result_type; // todo can this be removed?

        std::size_t operator()(const argument_type &arg) const { // todo correct implementation?
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
            vector<std::string> ZetaPosInNode;
            vector<t_meter> ZetasChange;
            vector<t_meter> ZetasChange_TZero;
            int numOfExternalFlows{0};
            bool nwt{false};
            bool initial_head{true};
            bool simpleDistance{false};
            bool simpleK{false};
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
             * @brief Overload ostream operator "<<" to return various node properties
             */
            friend std::ostream& operator<< (std::ostream& stream, const p_node& pNode) {
                stream << "Node properties:"
                << "\nID [large_num]: "<< pNode->get<large_num, ID>()
                << "\nSpatID [large_num]: " << pNode->get<large_num, SpatID>()
                << "\nLat [double]: " << pNode->get<double, Lat>()
                << "\nLon [double]: " << pNode->get<double, Lon>()
                << "\nLayer [int]: " << pNode->get<int, Layer>()
                << "\nStepModifier [-]: " << pNode->get<t_dim, StepModifier>().value()
                << "\nArea [m²]: " << pNode->get<t_s_meter, Area>().value()
                << "\nVerticalSize [m]: " << pNode->get<t_meter, VerticalSize>().value()
                << "\nElevation [m]: " << pNode->get<t_meter, Elevation>().value()
                << "\nTopElevation [m]: " << pNode->get<t_meter, TopElevation>().value()
                << "\nSlope [-]: " << pNode->get<t_dim, Slope>().value()
                << "\nEFolding [m]: " << pNode->get<t_meter, EFolding>().value()
                << "\nConfinement [bool]: " << pNode->get<bool, Confinement>()
                << "\nK [m/s]: " << pNode->get<t_vel, K>().value()
                << "\nAnisotropy [-]: " << pNode->get<t_dim, Anisotropy>().value()
                << "\nStepSize [d_time]: " << pNode->get<quantity < d_time>, StepSize>().value()
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
                << "\nVolumeOfCell [m³]: " << pNode->get<t_c_meter, VolumeOfCell>().value();

                unordered_map<NeighbourPosition, large_num> neighbourList = pNode->getListOfNeighbours();
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
                          large_num SpatID,
                          large_num ID,
                          t_vel K,
                          int stepModifier,
                          double aquiferDepth,
                          double anisotropy,
                          double specificYield,
                          double specificStorage,
                          bool confined,
                          bool densityVariable);

            virtual ~NodeInterface() = default;

            large_num getID() { return get<large_num, SpatID>(); }

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
             * @param slope_percent
             */
            void setSlope(double slope_percent) {
                set < t_dim, Slope > ((slope_percent / 100) * si::si_dimensionless);
                applyToAllLayers([slope_percent](NodeInterface *nodeInterface) {
                    try { // todo: should slope be added to the nodeInterface? currently only in PhysicalProperties
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

            void setMaxToeSlope(t_dim maxToeSlope) {
                set < t_dim, MaxToeSlope > (maxToeSlope);
            }

            void setMaxTipSlope(t_dim maxTipSlope) {
                set < t_dim, MaxTipSlope > (maxTipSlope);
            }

            void setDelnus(vector<t_dim> delnusVec){ set<vector<t_dim>, Delnus>(delnusVec);
                    // todo make this point to delnusVec (check all other node properties that are defined in config)
            }

            void setNusInZones(vector<t_dim> nusInZones){ set<vector<t_dim>, NusInZones>(nusInZones); }

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
                                       getLengthNeig(got), // length of neighbour node
                                       getLengthSelf(got), // length of this node
                                       //getWidthNeig(got), // width of neighbour node
                                       getWidthSelf(got), // width of this node
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
                            t_vol_t fluxCorrectionTop = getFluxCorrTop();
                            t_vol_t fluxCorrectionDown = getFluxCorrDown();

                            t_vol_t flow = conductance * (get<t_meter, HeadType>() - getAt<t_meter, HeadType>(got) -
                                                          fluxCorrectionTop - fluxCorrectionDown);
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
                } else if (type == FLOODPLAIN_DRAIN) {  // TODO adapt to rectangular nodes (change bottom)
                    externalFlows.insert(std::make_pair(type,
                                                        ExternalFlow(numOfExternalFlows, type,
                                                                        get<t_meter, Elevation>(),
                                                                        get<t_vel, K>() * get<t_meter, VerticalSize>(),
                                                                        bottom))); // todo add edge lenth in different directions
                } else { // RIVER, RIVER_MM, DRAIN, WETLAND, GLOBAL_WETLAND, LAKE, GENERAL_HEAD_BOUNDARY
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
             * @brief Add a zeta surface to the cell (bounded by elevation at top and cell bottom at bottom).
             * @param initialZeta the zeta surface height in meters
             */
            void addInitialZeta(t_meter height){
                NANChecker(height.value(), "height (in addInitialZeta)");

                // in SWI2: lines 660-680
                t_meter swismall = 0.001 * si::meter; // SWISMALL
                t_meter topOfNode = get<t_meter, Elevation>();
                t_meter bottomOfNode = get<t_meter, Elevation>() - get<t_meter, VerticalSize>();
                if (height > topOfNode - swismall) {
                    height = topOfNode;
                } else if (height < bottomOfNode + swismall){
                    height = bottomOfNode;
                }

                Zetas.push_back(height);
                //Zetas_TZero.push_back(height);
                ZetasChange.push_back(0 * si::meter); // todo improve?
                ZetasChange_TZero.push_back(0 * si::meter); // todo improve?
                //todo throw error if height not in correct order
            }

            /**
             * @brief Update zetas after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            virtual void setZeta(int localZetaID, t_meter height) {
                NANChecker(height.value(), "height (in setZeta)");
                if (localZetaID < Zetas.size()) {
                    Zetas[localZetaID] = height;
                }
            }

            /**
             * @brief Update zeta change after one or multiple inner iteration
             * @param localZetaID zeta surface id in this node
             * @param height zeta height
             */
            virtual void setZetaChange(int localZetaID, t_meter height){
                NANChecker(height.value(), "height (in setZetaChange)");
                if (localZetaID < ZetasChange.size()) {
                    ZetasChange[localZetaID] = height - Zetas[localZetaID];
                }
            }

            /**
             * @brief get zeta surface height
             * @return meter
             */
            t_meter getZeta(int localZetaID) noexcept { return Zetas[localZetaID];}

            /**
             * @brief get zeta surface change
             * @return meter
             */
            t_meter getZetaChange(int localZetaID) noexcept { return ZetasChange[localZetaID];}

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
             * - zoneOfSources: an integer between 1 and the number of zones (brackish/saline water)
             * - zoneOfSinks: 0 (fresh water)
             */
            void setZoneOfSinksAndSources(int zoneSinks, int zoneSources, int numOfZones) {
                if (zoneSinks > zoneSources) {
                    throw "Zone number of sinks cannot be larger than zone number of sources";
                }

                if (zoneSinks < zoneSources and zoneSinks > 0) {
                    throw "If zone number of sinks and sources are not equal (for simulating SGD), zone of sinks is 0";
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

            std::string getZetaPosInNode(int localZetaID){ return ZetaPosInNode[localZetaID]; }

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
             * Currently assumes only one recharge as external flow!
             * @param amount The new flow amount
             * @param Should the recharge in the dynamic rivers be locked or updated by this change?
             */
            void updateUniqueFlow(double amount, FlowType flow = RECHARGE, bool lock = true) {
                if (lock and flow == RECHARGE) {
                    if (hasTypeOfExternalFlow(RIVER_MM)) {
                        //get current recharge and lock it bevor setting new recharge
                        //in arid regions recharge might be 0
                        t_vol_t recharge{0 * si::cubic_meter / day};
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

            /**
             * @brief calculate the right hand side for zeta surface equation (b_zetas)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             */
            t_vol_t getZetaRHS(int localZetaID){ // todo: test
                t_vol_t porosityTerm = 0 * (si::cubic_meter / day);
                if (ZetaPosInNode[localZetaID] == "between" or localZetaID == 0) { // if IPLPOS == 0
                    porosityTerm = getEffectivePorosityTerm() * Zetas[localZetaID];
                }
                //LOG(userinfo) << "porosityTerm: " << porosityTerm.value() << std::endl;
                t_vol_t sourcesBelowZeta = getSourcesBelowZeta(localZetaID); // in SWI2 code: part of BRHS; in SWI2 doc: G or known source term below zeta
                //LOG(userinfo) << "sourcesBelowZeta: " << sourcesBelowZeta.value() << std::endl;
                t_vol_t pseudoSourceBelowZeta = getPseudoSourceBelowZeta(localZetaID);
                //LOG(userinfo) << "pseudoSourceBelowZeta: " << pseudoSourceBelowZeta.value() << std::endl;
                t_vol_t tipToeFlow = getTipToeFlow(localZetaID);
                //LOG(userinfo) << "tipToeFlow: " << tipToeFlow.value() << std::endl;
                t_vol_t out = - porosityTerm - sourcesBelowZeta + pseudoSourceBelowZeta + tipToeFlow;
                NANChecker(out.value(), "getZetaRHS");
                return out;
            }

            /**
             * @brief The effective porosity term for this node
             * @return square meter over time
             */
            t_s_meter_t getEffectivePorosityTerm(){ // computed independent of steady or transient flow (see SWI2 doc "Using the SWI2 Package")
                t_s_meter_t out = (get<t_dim, EffectivePorosity>() *
                                   get<t_meter, EdgeLengthLeftRight>() * get<t_meter, EdgeLengthFrontBack>()) /
                                  (day * get<t_dim, StepModifier>()); // todo add get<t_dim, StepModifier>()
                //LOG(debug) << "effective porosity: " << out.value() << std::endl;
                NANChecker(out.value(), "getEffectivePorosityTerm");
                return out;
            }

            /**
             * @brief The matrix entry for the left hand side of the zeta surface equation
             * @param localZetaID zeta surface id in this node
             * @return map <CellID,Conductance>
             */
            std::unordered_map<large_num, t_s_meter_t> getVDFMatrixEntries(int localZetaID) { // Question: Rename to getLeftHandSide_Zeta?
                // todo: test
                size_t numC = 5;
                std::unordered_map<large_num, t_s_meter_t> out;
                out.reserve(numC);
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
                auto delnus = get<vector<t_dim>, Delnus>();

                std::vector<t_s_meter_t> zoneConductances;
                t_s_meter_t zoneConductanceCum;
                t_s_meter_t zetaMovementConductance;

                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    zetaMovementConductance = 0 * (si::square_meter / day);
                    if (got == neighbours.end()) { // no neighbour at position
                    } else { // there is a neighbour at position
                        if (hasGHB()) { // not at boundary nodes
                        } else {
                            if (ZetaPosInNode[localZetaID] == "between" and
                                at(got)->ZetaPosInNode[localZetaID] == "between") {
                                zoneConductances = getZoneConductances(got);
                                zoneConductanceCum = getZoneConductanceCum(localZetaID, zoneConductances);
                                zetaMovementConductance += delnus[localZetaID] * zoneConductanceCum; // in SWI2: SWISOLCC/R
                                //LOG(debug) << "zoneConductanceCum = " << zoneConductanceCum.value() << std::endl;
                            }
                        }

                        NANChecker(zetaMovementConductance.value(), "zetaMovementConductance");
                        // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                        out[nodes->at(got->second)->get<large_num, ID>()] = move(zetaMovementConductance);
                    }
                }

                // To solve for zeta in this node, the conductances to neighbours and porosity term are used
                t_s_meter_t tmp_c = 0 * (si::square_meter / day);
                if (hasGHB()) {// At boundary nodes, tmp_c stays 0
                } else {
                    // subtract the conductances to neighbours (which were calculated above)
                    for (const auto &c: out) { tmp_c = tmp_c - c.second; }
                }
                // subtract effective porosity term
                tmp_c = tmp_c - getEffectivePorosityTerm();
                //LOG(debug) << "effectivePorosityTerm = " << getEffectivePorosityTerm().value() << std::endl;

                // add conductance of this node to out, the key in the unordered map is the ID of this node
                out[get<large_num, ID>()] = tmp_c;
                return out;
            };

            /**
             * @brief Dimensionless density (nus) at the top of this node
             * @return dimensionless
             * @note in SWI2: NUTOP, lines 1230-1249
             */
            t_dim getNusTop(){
                vector<t_dim> nusInZones = get<vector<t_dim>, NusInZones>();
                t_dim out = nusInZones.front();
                for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++){
                    if (ZetaPosInNode[localZetaID] == "top"){
                        auto delnus = get<vector<t_dim>, Delnus>();
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
                vector<t_dim> nusInZones = get<vector<t_dim>, NusInZones>();
                t_dim out = nusInZones.back();
                for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++){
                    if (ZetaPosInNode[localZetaID] == "bottom"){
                        auto delnus = get<vector<t_dim>, Delnus>();
                        out -= delnus[localZetaID];
                    }
                }
                return out;
            }

            /**
             * @brief set the zeta position in the node as at top/bottom/between the first and last zeta surface
             * @param localZetaID zeta surface id in this node
             * @note in SWI2: lines 4358-4383 (line 4362 is ignored since it does not make sense)
             */
            void setZetaPosInNode(int localZetaID){
                std::string tmp;
                /*if (localZetaID == 0) {
                    tmp = "between"; // SWI2 line 4362
                } else */if (Zetas[localZetaID] >= Zetas.front()){
                    tmp = "top";
                } else if (Zetas[localZetaID] <= Zetas.back() or
                           get<t_meter,Head>() < (get<t_meter,Elevation>() - get<t_meter,VerticalSize>())){
                    tmp = "bottom";
                } else { //if (Zetas[localZetaID] < Zetas.front() and Zetas[localZetaID] > Zetas.back()){
                    tmp = "between";
                } // todo: if nodes can be inactive: additional else if

                if (localZetaID >= ZetaPosInNode.size()){
                    ZetaPosInNode.push_back(tmp);
                } else {
                    ZetaPosInNode[localZetaID] = tmp;
                }
                //LOG(debug) << "ZetaPosInNode[localZetaID=" << localZetaID << "]=" << ZetaPosInNode[localZetaID] << std::endl;
            }

            /**
             * @brief The source flow below a zeta surface (for the right hand side in the zeta equation)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note like G in SWI2 doc, but without vertical leakage; in SWI2 code: lines 3523-3569
             * G = RHS (of flow, for constant density) - HCOF_(i,j,k,n)*h^(m)_(i,j,k) + (verticalLeakage_(i,j,k-1,n) - verticalLeakage_(i,j,k,n))
             */
            t_vol_t getSourcesBelowZeta(int localZetaID){ //todo: debug
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                if (hasGHB()) { return out; } // todo: compute BUFF with SSWI2_BDCH
                t_vol_t boundaryFlux = 0.0 * (si::cubic_meter / day);
                int zoneToUse = zoneOfSources; // zone of flow sources in node (for SGD: zoneOfSources > zoneOfSinks)
                //LOG(userinfo) << "zoneToUse: " << zoneToUse << std::endl;

                if (ZetaPosInNode[localZetaID] == "between" or localZetaID == 0) {
                    // Question: adapt functionality to deal with zone numbers over 100 from SWI2 code? lines 3548-3553

                    // if the new groundwater head is above or equal to the node bottom
                    if (get<t_meter, Head>() >= (get<t_meter, Elevation>() - get<t_meter, VerticalSize>())) { // lines 3532-3536
                        // get RHS of flow equation without VDF terms (pseudo source term and flux correction)
                        //t_vol_t RHSConstantDensity = get<t_vol_t, RHSConstantDensity_TZero>(); // in SWI2 code: RHSPRESWI
                        t_vol_t extFlows = -getQ(); //LOG(debug) << "extFlows: " << extFlows.value() << std::endl;
                        t_vol_t dewateredFlow = calculateDewateredFlow(); //LOG(debug) << "dewateredFlow: " << dewateredFlow.value() << std::endl;
                        t_vol_t rivers = calculateNotHeadDependandFlows(); //LOG(debug) << "rivers: " << rivers.value() << std::endl;
                        t_vol_t storageFlow =
                                getStorageCapacity() * (get<t_meter, Head_TZero>() / (day * get<t_dim, StepModifier>())); // todo add get<t_dim, StepModifier>()
                        if (steadyState) {
                            storageFlow = 0 * (si::cubic_meter / day);
                        }

                        t_vol_t RHSConstantDensity = extFlows + dewateredFlow - rivers - storageFlow;
                        // get HCOF of NEW time step
                        t_s_meter_t hcof = mechanics.getHCOF(steadyState,
                                                             get<t_dim, StepModifier>(),
                                                             getStorageCapacity(),
                                                             getP()); // todo add get<t_dim, StepModifier>()
                        // calculate the boundary flux
                        boundaryFlux = RHSConstantDensity - hcof * get<t_meter, Head>();
                        //LOG(userinfo) << "RHSConstantDensity: " << RHSConstantDensity.value() << std::endl;
                        //LOG(userinfo) << "hcof: " << hcof.value() << std::endl;
                        //LOG(userinfo) << "boundaryFlux: " << boundaryFlux.value() << std::endl;
                    }

                    if (zoneOfSources > zoneOfSinks and // if we intend to simulate submarine groundwater discharge
                        boundaryFlux > 0 * (si::cubic_meter / day)) { // and boundary flux is positive, so out of node
                        zoneToUse = zoneOfSinks; // use the zone of sinks (= top of aquifer = fresh water)
                    }
                    //LOG(userinfo) << "zoneToUse: " << zoneToUse << std::endl;
                    // Question: add 3539-3540 with "if (ibound < 0){q=-BUFF}"?

                    if (localZetaID <= zoneToUse and boundaryFlux != 0 * (si::cubic_meter / day)) {
                        out = boundaryFlux; // in SWI2 code, boundary flux is multiplied by "fact = thickb / thick" that depends on IZONENR if IZONENR > 100
                    }
                }

                NANChecker(out.value(), "getSourcesBelowZeta");
                //LOG(userinfo) << "getSourcesBelowZeta[" << localZetaID << "]: " << out.value() << std::endl;
                return out;
            }

            /**
             * @brief The pseudo source for a zeta surface (for the right hand side in the zeta equation)
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note in SWI2 code lines 3574-3635, using SSWI2_SD and SSWI2_SR)
             */
            t_vol_t getPseudoSourceBelowZeta(int localZetaID) { // todo test
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<vector<t_dim>, Delnus>();
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};

                // pseudo source term calculation (in 2 parts)
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) {
                        continue;
                    } else {
                        // calculating zone conductances for pseudo source term calculation
                        std::vector<t_s_meter_t> zoneConductances = getZoneConductances(got);

                        if ((ZetaPosInNode[localZetaID] == "between" and
                            at(got)->ZetaPosInNode[localZetaID] == "between") or localZetaID == 0) { // if IPLPOS == 0
                            //%% head part %%
                            t_s_meter_t zoneCondCum = getZoneConductanceCum(localZetaID,zoneConductances);
                            t_vol_t head_part = -zoneCondCum * (getAt<t_meter, Head>(got) - get<t_meter, Head>());
                            out += head_part;
                            //LOG(userinfo) << "head_part: " << head_part.value() << std::endl;
                            //LOG(userinfo) << "zoneCondCum: " << zoneCondCum.value() << std::endl;
                            //LOG(userinfo) << "getAt<t_meter, Head>(got): " << getAt<t_meter, Head>(got).value() << std::endl;
                            //LOG(userinfo) << "get<t_meter, Head>(): " << get<t_meter, Head>().value() << std::endl;

                            for (int zetaID = 0; zetaID < Zetas.size() - 1; zetaID++) {
                                //%% delnus part %%
                                if (zetaID != localZetaID) {
                                    t_s_meter_t zoneCondCumZetaID = getZoneConductanceCum(zetaID,
                                                                                          zoneConductances);
                                    t_vol_t delnus_part = -delnus[zetaID] *
                                                          (zoneCondCumZetaID * (at(got)->Zetas[zetaID] - Zetas[zetaID]));
                                    out += delnus_part;
                                    //LOG(userinfo) << "Zetas[zetaID]: " << Zetas[zetaID].value() << std::endl;
                                    //LOG(userinfo) << "at(got)->Zetas[zetaID]: " << at(got)->Zetas[zetaID].value() << std::endl;
                                    //LOG(userinfo) << "zoneCondCumZetaID: " << zoneCondCumZetaID.value() << std::endl;
                                    //LOG(userinfo) << "delnus_part: " << delnus_part.value() << std::endl;
                                }
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
            t_vol_t tipToeFlow(map_itter got, int localZetaID){
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                auto delnus = get<vector<t_dim>, Delnus>();

                std::vector<t_s_meter_t> zoneConductances = getZoneConductances(got);
                t_s_meter_t zoneCondCum = getZoneConductanceCum(localZetaID, zoneConductances);
                // %%head part %% for left/back neighbour
                t_vol_t head_part = zoneCondCum * (getAt<t_meter, Head>(got) - get<t_meter, Head>());
                out += head_part;
                //LOG(debug) << "head_part: " << head_part.value() << std::endl;
                //LOG(debug) << "with zoneCondCum = " << zoneCondCum.value() << std::endl;
                //LOG(debug) << "with getAt<t_meter, Head>(got): " << getAt<t_meter, Head>(got).value() << std::endl;
                //LOG(debug) << "get<t_meter, Head>(): " << get<t_meter, Head>().value() << std::endl;

                // %%delnus part %% for left/back neighbour
                for (int zetaID = 0; zetaID < Zetas.size() - 1; zetaID++) {
                    t_s_meter_t zoneCondCumZetaID = getZoneConductanceCum(zetaID,zoneConductances);
                    t_vol_t delnus_part = delnus[zetaID] *
                                          (zoneCondCumZetaID * (at(got)->Zetas[zetaID] - Zetas[zetaID]));
                    out += delnus_part;
                    //LOG(debug) << "delnus_part: " << delnus_part.value() << std::endl;
                    //LOG(debug) << "with zoneCondCumZetaID = " << zoneCondCumZetaID.value() << std::endl;
                    //LOG(debug) << "with Zetas[zetaID = " << zetaID << "] = " << Zetas[zetaID].value() << std::endl;
                    //LOG(debug) << "with at(got)->Zetas[zetaID = " << zetaID << "] = " << at(got)->Zetas[zetaID].value() << std::endl;
                }
                return out;
            }

            /**
             * @brief Specification of boundary condition at tips and toes
             * @param localZetaID zeta surface id in this node
             * @return volume per time
             * @note adapted from SWI2 code lines 3637-3703 (includes usage of SSWI2_QR and SSWI2_QC)
             */
            t_vol_t getTipToeFlow(int localZetaID){ // todo: debug
                t_vol_t out = 0.0 * (si::cubic_meter / day);
                //if (hasGHB()) { return out; } // return 0 at boundary nodes

                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::LEFT, NeighbourPosition::RIGHT,
                         NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::TOP, NeighbourPosition::DOWN};

                // tip and toe flow calculation (in 8 parts)
                if (ZetaPosInNode[localZetaID] == "between" or localZetaID == 0) { // if the Zeta surface at this node is between top and bottom of the node
                    for (const auto &position: possible_neighbours) {
                        map_itter got = neighbours.find(position);
                        if (got == neighbours.end()) { // do nothing if neighbour does not exist
                        } else {
                            if (at(got)->ZetaPosInNode[localZetaID] != "between" and localZetaID != 0) { // if the Zeta surface at neighbouring node is at top or bottom of the node
                                if ((position == NeighbourPosition::LEFT) or
                                    (position == NeighbourPosition::BACK)) {
                                    out -= tipToeFlow(got, localZetaID); // subtract tip/toe flow at LEFT/BACK neighbour

                                } else if ((position == NeighbourPosition::RIGHT) or
                                           (position == NeighbourPosition::FRONT)) {
                                    // add tipToeFlow for RIGHT and FRONT
                                    map_itter thisNode;
                                    if (position == NeighbourPosition::RIGHT){
                                        thisNode = at(got)->neighbours.find(NeighbourPosition::LEFT);
                                    } else if (position == NeighbourPosition::FRONT) {
                                        thisNode = at(got)->neighbours.find(NeighbourPosition::BACK);
                                    }
                                    t_vol_t ttf_right = at(got)->tipToeFlow(thisNode, localZetaID);
                                    LOG(userinfo) << "tipToeFlow(thisNode, localZetaID): " << ttf_right.value() << std::endl;
                                    out += ttf_right; // add tip/toe flow at RIGHT/FRONT neighbour

                                } else if (position == NeighbourPosition::TOP) {
                                    // vertical leakage to TOP neighbour
                                    // todo: debug/test
                                    vector<t_dim> nusInZones = get<vector<t_dim>, NusInZones>();
                                    if (getFluxCorrTop() < 0 * (si::cubic_meter / day) and
                                        at(got)->getNusBot() >= nusInZones[localZetaID] and
                                        getNusBot() >= at(got)->getNusBot()) { // IF ((qztop.LT.0).AND.(NUBOT(i,j,k-1).GE.NUS(iz)).AND.(NUBOT(j,i,k).GE.NUBOT(i,j,k-1))) THEN
                                        out += getFluxCorrTop(); // in SWI2: qztop
                                    }
                                } else if (position == NeighbourPosition::DOWN) {
                                    // vertical leakage to DOWN neighbour
                                    // todo: debug/test
                                    vector<t_dim> nusInZones = get<vector<t_dim>, NusInZones>();
                                    if(getFluxCorrDown() < 0 * (si::cubic_meter / day) and
                                        at(got)->getNusTop() < nusInZones[localZetaID] and
                                        getNusTop() <= at(got)->getNusTop()) { // IF ((qzbot.LT.0).AND.(NUTOP(i,j,k+1).LT.NUS(iz)).AND.(NUTOP(j,i,k).LE.NUTOP(i,j,k+1))) THEN
                                        continue;
                                    } else{
                                        out += getFluxCorrDown(); // in SWI2: qzbot
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
                for (int localZetaID = 0; localZetaID < Zetas.size() - 1; localZetaID++) {
                    //LOG(userinfo) << "Zetas[" << localZetaID << "]:" << Zetas[localZetaID].value() << std::endl;
                    //LOG(userinfo) << "Zetas[" << localZetaID+1 << "]:" << Zetas[localZetaID+1].value() << std::endl;

                    deltaZeta = Zetas[localZetaID] - Zetas[localZetaID + 1];
                    deltaZeta_neig = at(got)->Zetas[localZetaID] - at(got)->Zetas[localZetaID + 1];
                    if (deltaZeta <= (0 * si::meter) or deltaZeta_neig <= (0 * si::meter)){ // adapted from SWI2 code line 1149
                        zoneThickness = 0 * si::meter;
                    } else {
                        zoneThickness = ((getLengthNeig(got) * deltaZeta) / (getLengthNeig(got) + getLengthSelf(got))) +
                                        ((getLengthSelf(got) * deltaZeta_neig) / (getLengthNeig(got) + getLengthSelf(got)));
                        sumOfZoneThicknesses += zoneThickness;
                    }
                    NANChecker(zoneThickness.value(), "zoneThickness");
                    zoneThicknesses.push_back(zoneThickness);
                }

                // calculate the density zone conductances
                for (int localZetaID = 0; localZetaID < Zetas.size() - 1; localZetaID++) {
                    //LOG(userinfo) << "zoneThicknesses[" << localZetaID << "]:" << zoneThicknesses[localZetaID].value();

                    if (sumOfZoneThicknesses == (0 * si::meter)) { // adapted from SWI2 code line 1159
                        zoneConductance = 0 * si::square_meter / day; // zoneConductance is 0 if sum of thicknesses is 0
                    } else {
                        if (localZetaID > 0 and localZetaID < Zetas.size() - 2 and
                            ((ZetaPosInNode[localZetaID] != "between") or // adapted from SWI2 code lines 1212-1222
                            (at(got)->ZetaPosInNode[localZetaID] != "between") or
                            (ZetaPosInNode[localZetaID + 1] != "between") or
                            (at(got)->ZetaPosInNode[localZetaID + 1] != "between"))) {
                                // this section is reached if Zetas.size() >= 4 and numberOfZones >= 3
                                zoneConductance = 0 * si::square_meter / day; // zone is inactive or across layers
                        } else {

                            conductance = mechanics.calculateHarmonicMeanConductance(createDataTuple<Head>(got));
                            //LOG(userinfo) << "conductance:" << conductance.value();
                            zoneConductance = conductance * (zoneThicknesses[localZetaID] / sumOfZoneThicknesses);
                            //LOG(userinfo) << "zoneConductance[" << localZetaID << "] :" << zoneConductance.value();
                        }
                    }
                    //LOG(debug) << "zoneConductance:" << zoneConductance.value() << std::endl;
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

                for (int zoneID = localZetaID; zoneID < Zetas.size() - 1; zoneID++) {
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
                if (hasGHB()) { return out; } // return 0 at boundary nodes

                auto delnus = get<vector<t_dim>, Delnus>();
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};

                // pseudo source term calculation
                for (const auto &position: possible_neighbours) {
                    map_itter got = neighbours.find(position);
                    if (got == neighbours.end()) { // no neighbour at position
                        continue;
                    }
                    std::vector<t_s_meter_t> zoneConductances = getZoneConductances(got);

                    for (int localZetaID = 0; localZetaID < Zetas.size() - 1; localZetaID++){
                        t_s_meter_t zoneConductanceCum = getZoneConductanceCum(localZetaID, zoneConductances);
                        //LOG(userinfo) << "zoneConductanceCum (localZetaID: " << localZetaID << "): " << zoneConductanceCum.value() << std::endl;
                        if (delnus[localZetaID] > 0) {
                            out -= delnus[localZetaID] * (zoneConductanceCum *
                                                          (at(got)->Zetas[localZetaID] - Zetas[localZetaID]));
                        }
                    }
                }
                LOG(debug) << "getPseudoSourceNode: " << out.value() << std::endl;
                return out;
            }

            /**
             * @brief The flux correction term in vertical direction
             * @return volume per time
             * @note in SWI2 documentation: CV*BOUY; in code: QLEXTRA
             */
            t_vol_t getVerticalFluxCorrection(){ // todo debug/test
                vector<t_dim> nusInZones = get<vector<t_dim>, NusInZones>();
                t_meter headdiff = 0 * si::meter;
                t_vol_t out = 0 * (si::cubic_meter / day);

                // find the top neighbor
                map_itter got = neighbours.find(NeighbourPosition::TOP);
                if (got == neighbours.end()) {//No top node
                } else {//Current node has a top node
                    // first part of the flux correction term
                    for (int localZetaID = 0; localZetaID < Zetas.size() - 1; localZetaID++){
                        headdiff -= nusInZones[localZetaID] * (at(got)->Zetas[localZetaID+1] - at(got)->Zetas[localZetaID]); // Question: how to deal with this: in documentation is, BOUY is calculated with the simple sum (would be out +=), MODFLOW code for headdiff is as implemented (like out -=)
                    }
                    // second part of the flux correction term
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                    out = verticalConductance *
                          (headdiff + 0.5 * (at(got)->Zetas.back() - Zetas.front()) *
                           (at(got)->getNusBot() + getNusTop())); // Question: how to deal with this: in documentation, BOUY is calculated with a - between NUBOT and NUTOP, but in the code there is a - in the calculation of QLEXTRA
                    //LOG(userinfo) << "headdiff: " << headdiff.value() << std::endl;
                    //LOG(userinfo) << "verticalConductance: " << verticalConductance.value() << std::endl;
                }
                return out;
            }

            /**
             * @brief Calculates the upward vertical flux correction for variable density flow
             * @return volume per time
             * @note in SWI2 code: qztop
             */
            t_vol_t getFluxCorrTop() { // todo test
                t_vol_t out = 0 * (si::cubic_meter / day);
                map_itter got = neighbours.find(NeighbourPosition::TOP);
                if (got == neighbours.end()) { // no neighbour at position
                } else {
                    t_vol_t fluxFromTopNode = getVerticalFluxCorrection();
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                    out = verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(got)) - fluxFromTopNode;

                }
                return out;
            }

            /**
             * @brief Calculates the downward vertical flux correction for variable density flow
             * @return volume per time
             * @note in SWI2 code: qzbot
             */
            t_vol_t getFluxCorrDown() { // todo test
                t_vol_t out = 0 * (si::cubic_meter / day);
                map_itter got = neighbours.find(NeighbourPosition::DOWN);
                if (got == neighbours.end()) { // no neighbour at position
                } else {
                    t_vol_t fluxFromDownNode = at(got)->getVerticalFluxCorrection();
                    t_s_meter_t verticalConductance = mechanics.calculateVerticalConductance(createDataTuple(got));
                    //LOG(userinfo) << "fluxFromDownNode: " << fluxFromDownNode.value() << std::endl;
                    //LOG(userinfo) << "verticalConductance: " << verticalConductance.value() << std::endl;

                    out = verticalConductance * (get<t_meter, Head>() - getAt<t_meter, Head>(got)) + fluxFromDownNode;
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
                    t_meter head = get<t_meter, Head>();
                    t_meter bottomOfNode = get<t_meter, Elevation>() - get<t_meter, VerticalSize>();
                    t_meter topOfNode = get<t_meter, Elevation>();
                    t_meter updatedZeta;
                    // if groundwater head is BELOW the top of the node
                    if (head < topOfNode) {
                        // if groundwater head is ABOVE the bottom of the node
                        if (head > bottomOfNode) {
                            updatedZeta = head;
                        // if groundwater head is BELOW OR EQUAL to the bottom of the node
                        } else { // head <= bottomOfNode
                            updatedZeta = bottomOfNode;
                        }
                        // update the first zeta surface
                        Zetas.front() = updatedZeta;

                        // update all other zeta surfaces that are ABOVE the updated first zeta surface
                        for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++) {
                            if (Zetas[localZetaID] > updatedZeta) {
                                Zetas[localZetaID] = updatedZeta;
                            }
                        }
                    // if groundwater head is ABOVE OR EQUAL to the top of the node
                    } else { // head >= topOfNode
                        // clip zeta to the top of the node
                        Zetas.front() = topOfNode;
                    }
                }
            }


            /**
             * @brief vertical movement of zeta surfaces through top of this node
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_VERTMOVE
             */
            void verticalZetaMovement() { // todo: debug
                // skip nodes where head is below the bottom of the node
                t_vol_t fluxCorrectionTop; // in SWI2: qztop
                t_s_meter_t verticalConductanceTop;
                t_meter deltaZeta;
                map_itter top = neighbours.find(NeighbourPosition::TOP);
                // skip nodes that do not have a top neighbour
                if (top == neighbours.end()) { // no neighbour at position
                } else {
                    if (get<t_meter, Head>() >= (get<t_meter, Elevation>() - get<t_meter, VerticalSize>()) and
                        getAt<t_meter, Head>(top) >= (getAt<t_meter, Elevation>(top) - getAt<t_meter, VerticalSize>(top))) {

                        // calculate flux through the top
                        fluxCorrectionTop = getFluxCorrTop();
                        LOG(userinfo) << "fluxCorrectionTop: " << fluxCorrectionTop.value() << std::endl;
                        for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++) {
                            // zeta only moves through the top of a node if there is a ZETA surface
                            // - at the top of the current node (in SWI2: IPLPOS_(i,j,k,n) = 1)
                            // - AND at the bottom of the top node (in SWI2: IPLPOS_(i,j,k-1,n) = 2)
                            if (ZetaPosInNode[localZetaID] == "top" and
                                at(top)->ZetaPosInNode[localZetaID] == "bottom") {
                                // if vertical flux through the top of the node is positive...
                                if (fluxCorrectionTop > (0 * si::cubic_meter / day)) {
                                    deltaZeta = (fluxCorrectionTop * (day * get<t_dim, StepModifier>())) /
                                                (get<t_s_meter, Area>() * getAt<t_dim, EffectivePorosity>(top)); // todo add get<t_dim, StepModifier>() again
                                    // ...lift zeta height of the lowest zeta surface in this node
                                    at(top)->Zetas[localZetaID] = at(top)->Zetas.back() + deltaZeta;
                                    // if vertical flux through the top of the node is negative...
                                } else if (fluxCorrectionTop < (0 * si::cubic_meter / day)) {
                                    deltaZeta = (fluxCorrectionTop * (day * get<t_dim, StepModifier>())) /
                                                (get<t_s_meter, Area>() * get<t_dim, EffectivePorosity>()); // todo add get<t_dim, StepModifier>() again
                                    // ...lower zeta height of the highest zeta surface in this node
                                    Zetas[localZetaID] = Zetas.front() + deltaZeta;
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
                if (hasGHB()) { return; } // do nothing at boundary nodes
                std::unordered_map<NeighbourPosition, NeighbourPosition> oppositePositions;
                oppositePositions[NeighbourPosition::BACK] = NeighbourPosition::FRONT;
                oppositePositions[NeighbourPosition::FRONT] = NeighbourPosition::BACK;
                oppositePositions[NeighbourPosition::LEFT] = NeighbourPosition::RIGHT;
                oppositePositions[NeighbourPosition::RIGHT] = NeighbourPosition::LEFT;

                t_meter edgeLength_self;
                t_meter edgeLength_neig;
                t_meter edgeLength_neig_opp; // edge length of opposite neighbour node
                t_meter maxDeltaZeta; // maximum height difference of a zeta surface n between adjacent nodes
                t_meter referenceElevation; // top/bottom of the node for tip/toe tracking
                t_dim slopeAdjustmentFraction = 0.1 * si::si_dimensionless; // slope adjustment fraction (in SWI2: ALPHA) // todo move to config
                t_dim minDepthThreshold = 0.1 * si::si_dimensionless; // minimum depth threshold for zeta surface toes (in SWI2: BETA) // todo move to config
                t_dim effPor_self; // effective porosity of this node
                t_dim effPor_neig; // effective porosity of neighbouring node
                t_dim effPor_neig_opp; // effective porosity of opposite neighbouring node
                t_meter zetaChange_self; // potential zeta height adjustment for this node
                t_meter zetaChange_neig; // potential zeta height adjustment for this node

                std::unordered_map<string, t_dim> maxSlopes;
                maxSlopes["Toe"] = get<t_dim, MaxToeSlope>();
                maxSlopes["Tip"] = get<t_dim, MaxTipSlope>();
                std::forward_list<NeighbourPosition> possible_neighbours =
                        {NeighbourPosition::BACK, NeighbourPosition::FRONT,
                         NeighbourPosition::LEFT, NeighbourPosition::RIGHT};
                for (const auto &maxSlope : maxSlopes) {
                    for (const auto &position : possible_neighbours) {
                        map_itter got = neighbours.find(position);
                        map_itter got_opp = neighbours.find(oppositePositions[position]);
                        if (got == neighbours.end() or got_opp == neighbours.end()) { // no neighbour at position or opposite side
                        } else {
                            for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++) {
                                if (ZetaPosInNode[localZetaID] == "between") {
                                    if ((maxSlope.first == "Toe" and at(got)->ZetaPosInNode[localZetaID] == "bottom") or
                                        (maxSlope.first == "Tip" and at(got)->ZetaPosInNode[localZetaID] == "top")) {

                                        // get length of the edges
                                        edgeLength_self = getLengthSelf(got);
                                        edgeLength_neig = getLengthNeig(got);
                                        // get max delta of zeta between nodes
                                        maxDeltaZeta = 0.5 * (edgeLength_self + edgeLength_neig) * maxSlope.second;
                                        //LOG(debug) << "maxDeltaZeta: " << maxDeltaZeta.value() << std::endl;

                                        // raise/lower zeta surfaces
                                        effPor_self = get<t_dim, EffectivePorosity>();
                                        effPor_neig = getAt<t_dim, EffectivePorosity>(got);
                                        // if tracking tip/toe: raise/lower this zeta surface in this node (ZETA_(i,j,k,n))
                                        zetaChange_self = slopeAdjustmentFraction * maxDeltaZeta *
                                                          ((effPor_neig * edgeLength_neig) /
                                                           ((effPor_self * edgeLength_self) +
                                                            (effPor_neig * edgeLength_neig)));
                                        // if tracking tip/toe: lower/raise this zeta surface in neighbouring node (e.g. ZETA_(i,j+1,k,n))
                                        zetaChange_neig = slopeAdjustmentFraction * maxDeltaZeta *
                                                          ((effPor_self * edgeLength_self) /
                                                           ((effPor_self * edgeLength_self) +
                                                            (effPor_neig * edgeLength_neig)));

                                        if (maxSlope.first == "Toe" and
                                            at(got)->ZetaPosInNode[localZetaID] == "bottom") {
                                            t_meter deltaZeta = abs(Zetas[localZetaID] - at(got)->Zetas.back());
                                            //LOG(debug) << "deltaZeta (toe): " << deltaZeta.value() << std::endl;

                                            if (deltaZeta > maxDeltaZeta) {
                                                Zetas[localZetaID] = Zetas[localZetaID] - zetaChange_self;
                                                //LOG(debug) << "deltaZeta (toe): " << -zetaChange_self.value() << std::endl;
                                                at(got)->Zetas[localZetaID] = at(got)->Zetas.back() + zetaChange_neig;
                                                //LOG(debug) << "deltaZeta (toe): " << zetaChange_neig.value() << std::endl;
                                            }
                                        } else if (maxSlope.first == "Tip" and
                                                   at(got)->ZetaPosInNode[localZetaID] == "top") {
                                            t_meter deltaZeta = abs(at(got)->Zetas.front() - Zetas[localZetaID]);
                                            //LOG(debug) << "deltaZeta (tip): " << deltaZeta.value() << std::endl;
                                            if (deltaZeta > maxDeltaZeta) {
                                                Zetas[localZetaID] = Zetas[localZetaID] + zetaChange_self;
                                                //LOG(debug) << "zetaChange_self (tip): " << zetaChange_self.value() << std::endl;
                                                at(got)->Zetas[localZetaID] = at(got)->Zetas.front() - zetaChange_neig;
                                                //LOG(debug) << "zetaChange_neig (tip): " << -zetaChange_neig.value() << std::endl;
                                            }
                                        }

                                        if ((Zetas[localZetaID] - Zetas.back()) <
                                            (minDepthThreshold * zetaChange_neig)) {
                                            if (at(got_opp)->ZetaPosInNode[localZetaID] == "between") {
                                                // change zeta in other direction neighbour
                                                edgeLength_neig_opp = getLengthNeig(got_opp);
                                                effPor_neig_opp = getAt<t_dim, EffectivePorosity>(got_opp);
                                                at(got_opp)->Zetas[localZetaID] = at(got_opp)->Zetas[localZetaID] +
                                                                                  ((Zetas[localZetaID] - Zetas.back()) *
                                                                                   (edgeLength_self * effPor_self) /
                                                                                   (edgeLength_neig_opp *
                                                                                    effPor_neig_opp));
                                                Zetas[localZetaID] = Zetas.back();
                                            }
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
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ZETACLIP
             */
            void clipInnerZetas() {
                for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++) {
                    //LOG(debug) << "Zetas[localZetaID]: " << Zetas[localZetaID].value() << std::endl;

                    if (ZetaPosInNode[localZetaID] == "between") {
                        if (Zetas[localZetaID] < Zetas.back()) { Zetas[localZetaID] = Zetas.back(); }

                        if (Zetas[localZetaID] > Zetas.front()) { Zetas[localZetaID] = Zetas.front(); }
                    }
                }
            }

            /**
             * @brief Adjust zeta surfaces if they have crossed or their height difference is below a minimum value
             * @param localZetaID zeta surface id in this node
             * @note in SWI2 code: SSWI2_ZETACROSS
             * todo: debug with this info:
             *  n | n2      | n3
             *  1 | -       | -
             *  2 | 1       | 1,2,3
             *  3 | 2,1     | 2,3,4; 1,2,3,4
             *  4 | 3,2,1   | 3,4,5; 2,3,4,5; 1,2,3,4,5
             */
            void correctCrossingZetas(){
                t_meter zetaDifferenceCap = 0.001 * si::meter; // todo move to config
                t_meter zetaAverage;
                t_meter zetaSum;
                int n2_min;
                int n2_max;
                for (int localZetaID = 1; localZetaID < Zetas.size() - 2; localZetaID++) {
                    // if zeta surface n is very close to or lower than the zeta surface that SHOULD be below (n+1)
                    if (Zetas[localZetaID] - Zetas[localZetaID + 1] < zetaDifferenceCap) {
                        // make the zeta height of both surfaces their average
                        zetaAverage = 0.5 * (Zetas[localZetaID] + Zetas[localZetaID + 1]);
                        Zetas[localZetaID] = zetaAverage;
                        Zetas[localZetaID + 1] = zetaAverage;
                        // if there are zeta surfaces above that could possibly be crossing
                        if (localZetaID >= 2) {
                            for (int x = 1; x < localZetaID - 1; x++) {
                                // create a range from n_min (n-x) to n_max (n+1)
                                n2_min = localZetaID - x;
                                n2_max = localZetaID + 1;
                                // if a zeta surface above (n-x) crosses or is very close to zeta surface below (n+1)
                                //  (which is now at the same, averaged, height of zeta surface n)
                                if (Zetas[n2_min] - Zetas[n2_max] < zetaDifferenceCap) {
                                    // calculate the average height from that zeta surface above (n-x) until
                                    //  the zeta surface below (n+1)
                                    for (int n2 = n2_min; n2 <= n2_max; n2++) { zetaSum += Zetas[n2]; }
                                    zetaAverage = zetaSum / ((n2_max - n2_min) *
                                                             si::si_dimensionless); // todo check whether that is true
                                    // set every zeta surface between n-1 and n+1 to the averaged value
                                    for (int n2 = n2_min; n2 <= n2_max; n2++) { Zetas[n2] = zetaAverage; }
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
                if (hasGHB()) { return; } // do nothing at boundary nodes
                t_meter swiLock = 0.001 * si::meter;
                t_dim slopeAdjustmentFraction = 0.1 * si::si_dimensionless; // slope adjustment fraction (in SWI2: ALPHA)
                t_meter maxDeltaZeta;
                t_meter maxAdjustment_self;
                t_meter maxAdjustment_neig;
                t_meter deltaZeta_self;
                t_meter deltaZeta_neig;
                t_dim effectivePorosity_self;
                t_dim effectivePorosity_neig;
                for (int localZetaID = 1; localZetaID < Zetas.size() - 1; localZetaID++) {
                    if (ZetaPosInNode[localZetaID] == "top" or ZetaPosInNode[localZetaID] == "bottom") {
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
                                deltaZeta_self = maxDeltaZeta * (effectivePorosity_neig * getLengthNeig(got)) /
                                                 (effectivePorosity_self * getLengthSelf(got) +
                                                  effectivePorosity_neig * getLengthNeig(got));
                                deltaZeta_neig = maxDeltaZeta * (effectivePorosity_self * getLengthSelf(got)) /
                                                 (effectivePorosity_self * getLengthSelf(got) +
                                                  effectivePorosity_neig * getLengthNeig(got));

                                // if a zeta surface is at the BOTTOM of this node and at the TOP of the neighbour
                                // else if a zeta surface is at the TOP of this node and at the BOTTOM of the neighbour
                                if (ZetaPosInNode[localZetaID] == "bottom" and // IPLPOS_self = 2
                                    at(got)->ZetaPosInNode[localZetaID] == "top") { // IPLPOS_neig = 1
                                    // adjust zeta surface heights
                                    Zetas[localZetaID] = Zetas[localZetaID] + deltaZeta_self;
                                    at(got)->Zetas[localZetaID] = Zetas[localZetaID] - deltaZeta_neig;
                                } else if (ZetaPosInNode[localZetaID] == "top" and // IPLPOS_self = 1
                                           at(got)->ZetaPosInNode[localZetaID] == "bottom") {// IPLPOS_neig = 2
                                    // adjust zeta surface heights
                                    Zetas[localZetaID] = Zetas[localZetaID] - deltaZeta_self;
                                    at(got)->Zetas[localZetaID] = Zetas[localZetaID] + deltaZeta_neig;
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
                if (got->first == NeighbourPosition::LEFT or got->first == NeighbourPosition::RIGHT){
                    return getAt<t_meter, EdgeLengthFrontBack>(got);
                } else if (got->first == NeighbourPosition::FRONT or got->first == NeighbourPosition::BACK){
                    return getAt<t_meter, EdgeLengthLeftRight>(got);
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            /**
             * @brief The length of this node, in direction to the neighbouring node
             * @param got neighbour node
             * @return meter
             */
            t_meter getLengthSelf(map_itter got){
                if (got->first == NeighbourPosition::LEFT or got->first == NeighbourPosition::RIGHT){
                    return get<t_meter, EdgeLengthFrontBack>();
                } else if (got->first == NeighbourPosition::FRONT or got->first == NeighbourPosition::BACK){
                    return get<t_meter, EdgeLengthLeftRight>();
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
                }
            }

            /**
             * @brief The width of this node, perpendicular to the direction to the neighbouring node
             * @param got neighbour node
             * @return meter
             */
            t_meter getWidthSelf(map_itter got) {
                if (got->first == NeighbourPosition::LEFT or got->first == NeighbourPosition::RIGHT){
                    return get<t_meter, EdgeLengthLeftRight>();
                } else if (got->first == NeighbourPosition::FRONT or got->first == NeighbourPosition::BACK){
                    return get<t_meter, EdgeLengthFrontBack>();
                } else {
                    throw "Horizontal neighbour (LEFT/RIGHT or FRONT/BACK) required as input";
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
                        if (flow.second.flowIsHeadDependant(head)) {
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
                        if (not flow.second.flowIsHeadDependant(head)) {
                            out += flow.second.getP(eq_head, head, recharge, slope, eqFlow) * get<t_dim, StepModifier>() *
                                   flow.second.getBottom();
                        }
                    }
                }
                return out;
            }

            /**
             * @brief The Jacobian entry for the cell (NWT approach)
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

                        // add conductance to out, the key in the unordered map is the ID of the neighbouring node
                        // (used to solve for the head at the neighbouring node)
                        out[nodes->at(got->second)->get<large_num, ID>()] = move(conduct);
                    }
                }

                // To solve for the head at this node, the conductances to neighbours and HCOF are used
                t_s_meter_t tmp_c = 0; // Question: multiply with units?
                // subtract the conductances to neighbours (which were calculated above)
                for (const auto &c : out) { tmp_c = tmp_c - c.second; }
                // add HCOF
                t_s_meter_t hcof = mechanics.getHCOF(steadyState,
                                                     get<t_dim, StepModifier>(),
                                                     getStorageCapacity(),
                                                     getP());
                tmp_c = tmp_c + hcof;
                // check for nan
                NANChecker(tmp_c.value(), "HCOF");  //if(tmp_c.value() == 0){LOG(numerics) << "HCOF term is 0";}

                // add resulting conductance to solve for the head at this node to out
                out[get<large_num, ID>()] = tmp_c;
                return out;
            };

            /**
             * @brief The right hand side for constant density groundwater conditions
             * @return volume per time
             */
            t_vol_t getRHSConstantDensity(){
                t_vol_t extFlows = -getQ();
                //LOG(userinfo) << "extFlows: " << extFlows.value() << std::endl;
                t_vol_t dewateredFlow = calculateDewateredFlow();
                //LOG(userinfo) << "dewateredFlow: " << dewateredFlow.value() << std::endl;
                t_vol_t rivers = calculateNotHeadDependandFlows();
                //LOG(userinfo) << "rivers: " << rivers.value() << std::endl;
                t_vol_t storageFlow =
                        getStorageCapacity() * (get<t_meter, Head_TZero>() / (day * get<t_dim, StepModifier>()));
                //LOG(userinfo) << "storageFlow: " << storageFlow.value() << std::endl;
                if (steadyState) {
                    storageFlow = 0 * (si::cubic_meter / day);
                }

                t_vol_t out = extFlows + dewateredFlow - rivers - storageFlow; //LOG(debug) << "RHS constant density: " << out.value() << std::endl;
                NANChecker(out.value(), "RHS constant density");
                return out;
            }

            /**
             * @brief The right hand side of the flow equation
             * @return volume per time
             */
            t_vol_t getRHS() {
                t_vol_t out = getRHSConstantDensity();
                //LOG(userinfo) << "getRHSConstantDensity: " << out.value() << std::endl;

                if (get<bool, DensityVariable>()) {
                    // save constant density RHS (without variable density terms) for calculation of zeta movement
                    set<t_vol_t, RHSConstantDensity_TZero>(out);

                    // calculate variable density terms Pseudo-Source Flow
                    t_vol_t pseudoSourceNode = getPseudoSourceNode();
                    //LOG(userinfo) << "pseudoSourceNode: " << pseudoSourceNode.value() << std::endl;

                    // calculate Vertical Flux Correction (for top and down neighbour)
                    t_vol_t fluxCorrectionTop = getFluxCorrTop();
                    //LOG(userinfo) << "fluxCorrectionTop: " << fluxCorrectionTop.value() << std::endl;
                    t_vol_t fluxCorrectionDown = getFluxCorrDown();
                    //LOG(userinfo) << "fluxCorrectionDown: " << fluxCorrectionDown.value() << std::endl;

                    out += pseudoSourceNode + fluxCorrectionTop + fluxCorrectionDown;
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

            void setHead(t_meter delta) noexcept { // Question: change to "addDeltaToHead"?
                __setHead(delta);
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
             * @brief Calculate the lateral flow velocity
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
             * A negative value represents a flow in the opposite direction.
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
                         large_num SpatID,
                         large_num ID,
                         t_vel K,
                         int stepmodifier,
                         double aquiferDepth,
                         double anisotropy,
                         double specificYield,
                         double specificStorage,
                         bool confined,
                         bool densityVariable)
                    : NodeInterface(nodes, lat, lon, area, edgeLengthLeftRight, edgeLengthFrontBack, SpatID, ID, K,
                                    stepmodifier, aquiferDepth, anisotropy, specificYield, specificStorage, confined, densityVariable) {}
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
                    0.3 * (si::meter / day), 1, 100, 10, 0.15, 0.000015, true, true) {}

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
