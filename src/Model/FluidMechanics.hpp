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
#ifndef GLOBAL_FLOW_FLUID_MECHANICS_HPP
#define GLOBAL_FLOW_FLUID_MECHANICS_HPP

#include "Units.hpp"
#include <unordered_map>

namespace GlobalFlow {
namespace Model {

/**
 * @class FluidMechanics
 * Provides helper functions for conductance calulcations
 */
class FluidMechanics {
    public:
        FluidMechanics() {}

        /**
         * Used to calculate if a cell is dry
         */
        t_meter
        calcDeltaV(t_meter head, t_meter elevation, t_meter depth) noexcept;

        /**
         * @brief Calculates the horizontal flow between two nodes
         * @param k_neig K-value of neighbouring node
         * @param edgeLength_neig The edge size of the neighbouring node (assumes rectengular cells)
         * @return A weighted conductance value for the flow between two nodes
         *
         * Calculates the harmonic mean conductance between two nodes.
         * C = 2*EdgeLenght_1 * (TR_1 * TR_2) / (TR_1 * EdgeLenght_1 + TR_2 * EdgeLenght_2)
         */
        t_s_meter_t calculateHarmonicMeanConductance
                (t_vel k_neig,
                 t_vel k_self,
                 t_meter edgeLength_neig,
                 t_meter edgeLength_self,
                 t_meter head_neig,
                 t_meter head_self,
                 t_meter ele_neig,
                 t_meter ele_self,
                 t_meter deltaV_neig,
                 t_meter deltaV_self,
                 bool confined) noexcept;

        /**
         * Simple smoother function to buffer iteration steps in NWT approach
         * @param elevation
         * @param verticalSize
         * @param head
         * @return
         */
        double smoothFunction__NWT(t_meter elevation, t_meter verticalSize, t_meter head);

        /**
         * Get the coeffiecients for storage and P components
         * @param steadyState
         * @param stepModifier
         * @param storageCapacity
         * @param P
         * @return
         */
        t_s_meter_t getHCOF(bool steadyState, quantity<Dimensionless> stepModifier,
                t_s_meter storageCapacity, t_s_meter_t P) noexcept;

        /**
         * Calculates the vertical flow between two nodes
         * @param k_vert_neig
         * @param k_vert_self
         * @param verticalSize_self
         * @param head_self
         * @param elevation_self
         * @param area_self
         * @param elevation_neig
         * @param depth_neig
         * @param head_neig
         * @param confined
         * @return
         */
        t_s_meter_t calculateVerticalConductance(t_vel k_vert_neig,
                                     t_vel k_vert_self,
                                     t_meter verticalSize_self,
                                     t_meter head_self,
                                     t_meter elevation_self,
                                     t_s_meter area_self,
                                     t_meter elevation_neig,
                                     t_meter depth_neig,
                                     t_meter head_neig,
                                     bool confined) noexcept;

        /**
         * Calculate derivates for NWT approach
         * @param elevation
         * @param verticalSize
         * @param head
         * @return
         */
        double getDerivate__NWT(t_meter elevation,
                         t_meter verticalSize,
                         t_meter head);

};

}
}//ns
#endif
