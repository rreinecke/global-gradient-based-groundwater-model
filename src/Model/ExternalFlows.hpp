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

#ifndef GLOBAL_FLOW_EXTERNALFLOWS_HPP
#define GLOBAL_FLOW_EXTERNALFLOWS_HPP

#include "Units.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {
        /**
         * @enum External Flow Types
         * Flow types as defined in MODFLOW:
         *  - RECHARGE
         *  - RIVER
         *  - DRAIN
         *  - GENERAL_HEAD_BOUNDARY
         *
         * Additional flows:
         *
         * RECHARGE:
         *
         * RIVER_MM:
         *  dynamic river conductance as defined by Miguez-Macho 2007
         *
         * FLOODPLAIN_DRAIN:
         *  as defined in Inge de Graaf. 2014
         *
         *
         * LAKE, WETLAND
         *  similar to modflow river definition

         */
        enum FlowType : int {
            RECHARGE = 1,
            FAST_SURFACE_RUNOFF,    // 2
            NET_ABSTRACTION,        // 3
            EVAPOTRANSPIRATION,     // 4
            RIVER,                  // 5
            RIVER_MM,               // 6
            DRAIN,                  // 7
            FLOODPLAIN_DRAIN,       // 8
            WETLAND,                // 9
            GLOBAL_WETLAND,         // 10
            LAKE,                   // 11
            GLOBAL_LAKE,            // 12
            GENERAL_HEAD_BOUNDARY   // 13
        };

        struct FlowTypeHash {
            template<typename T>
            std::size_t operator()(T t) const {
                return static_cast<std::size_t>(t);
            }
        };

        /**
         * @class ExternalFlow
         *
         * TODO add flow equation here
         */
        class ExternalFlow {
        public:
            /**
             * @brief Constructor for RIVER, RIVER_MM, DRAIN, WETLAND, GLOBAL_WETLAND, LAKE, GENERAL_HEAD_BOUNDARY
             * @param id
             * @param type
             * @param flowHead
             * @param cond
             * @param bottomElev
             */
            ExternalFlow(int id,
                         FlowType type,
                         t_meter flowHead,
                         t_s_meter_t cond,
                         t_meter bottomElev)
                    : ID(id), type(type), flowHead(flowHead), conductance(cond), bottomElev(bottomElev) {}

            /**
             * @brief Constructor for RECHARGE, FAST_SURFACE_RUNOFF and NET_ABSTRACTION
             * @param id
             * @param flow
             * @param type
             */
            ExternalFlow(int id, t_vol_t flow, FlowType type)
                    : ID(id), type(type), flowHead(0), conductance(0), bottomElev(0), special_flow(flow) {}

            /**
             * @brief Constructor for Evapotranspiration
             * @param id
             * @param flowHead
             * @param bottomElev
             * @param evapotrans
             * @return
             */
            ExternalFlow(int id, t_meter flowHead, t_meter bottomElev, t_vol_t evapotrans)
                    : ID(id), type(EVAPOTRANSPIRATION), flowHead(0), conductance(0), bottomElev(0),
                      special_flow(evapotrans) {}

            /**
             * Check if flow can be calculated on the right hand side
             * @param head The current hydraulic head
             * @return Bool
             */
            bool isFlowHeadDependent(t_meter gw_head) const noexcept {
                return (gw_head > bottomElev);
            }

            /**
             * The head dependent part of the external flow equation:
             * This is the total conductance of all head-dependent external source terms in a cell
             * @param head The current hydraulic head
             * @param eq_head The equilibrium head
             * @param recharge The current recharge
             * @param eqFlow
             * @return
             */
            t_s_meter_t getP(t_meter head,
                             t_meter eq_head,
                             t_vol_t recharge,
                             t_vol_t eqFlow) const noexcept;

            /**
             * The head independent part of the external flow equation:
             * This is the total specified external source term
             * @param head
             * @param eq_head
             * @param recharge
             * @param eqFlow
             * @return
             */
            t_vol_t getQ(t_meter head,
                         t_meter eq_head,
                         t_vol_t recharge,
                         t_vol_t eqFlow) const noexcept;


            FlowType getType() const noexcept { return type; }

            t_meter getBottomElev() const noexcept { return bottomElev; }

            t_vol_t getRecharge() const noexcept { return special_flow; }

            t_meter getFlowHead() const noexcept { return flowHead; }

            t_s_meter_t getDyn(t_vol_t current_recharge,
                               t_meter eq_head,
                               t_meter head,
                               t_vol_t eq_flow) const noexcept {
                t_s_meter_t out = calcERC(current_recharge, eq_head, head, eq_flow);
                return out;
            }

            t_meter getRiverDiff(t_meter eqHead) const noexcept;

            t_s_meter_t getConductance() const noexcept { return conductance; }

            t_s_meter_t getInitConductance() const noexcept { return initConductance; }

            double getRiverDepthSteadyState() {return RiverDepthSteadyState; }

            //void setInitConductance(double initCond) { initConductance = initCond * boost::units::quantity<MeterSquaredPerTime>(); }
            void setInitConductance(double initCond) { initConductance = initCond * (si::square_meter / day); }

            void setRiverDepthSteadyState (double RiverDepth) {RiverDepthSteadyState = RiverDepth;}

            int getID() const noexcept { return ID; }

            void setMult(double mult) {
                this->mult = mult;
                return;
            }

            void setLock() { lock_recharge = true; }

            bool getLock() { return lock_recharge; }

            void setLockRecharge(t_vol_t re) { locked_recharge = re; }

            t_vol_t getLockRecharge() { return locked_recharge; }

            void setLockConduct(t_s_meter_t c) { locked_conductance = c; }

            t_s_meter_t getLockConduct() { return locked_conductance; }

            void getERC(t_vol_t current_recharge,
                    t_meter eq_head,
                    t_meter current_head,
                    t_vol_t eq_flow) { locked_conductance = calcERC(current_recharge,eq_head,current_head,eq_flow); };

        private:
            const int ID;
            const FlowType type;
            const t_meter flowHead;
            const t_s_meter_t conductance; //for special_flow same as q
            const t_vol_t special_flow;
            const t_meter bottomElev;
            t_dim mult{1 * si::si_dimensionless}; //Multiplier only used for SA
            t_s_meter_t initConductance = 0. * (si::square_meter / day);
            double RiverDepthSteadyState = -99.;
            t_vol_t locked_recharge;
            t_s_meter_t locked_conductance;
            bool lock_recharge{false};

            t_vol_t
            calculateFloodplainDrainage(t_meter head) const noexcept;

            /**
            * Calculate ERC (must be repeated every time step)
            * ERC = (GW_Recharge + eq_flow) / (eq_head - river_elevation)
            */
            t_s_meter_t
            calcERC(t_vol_t current_recharge,
                    t_meter eq_head,
                    t_meter current_head,
                    t_vol_t eq_flow) const noexcept;
        };

    }
}//ns
#endif //EXTERNALFLOWS_HPP