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
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {

        using FlowInputHor = std::tuple<t_vel, t_vel, t_meter, t_meter, t_meter, t_meter, t_meter, t_meter, t_meter, t_meter, t_meter, bool>;
        using FlowInputVert = std::tuple<t_vel, t_vel, t_meter, t_meter, t_meter, t_meter, t_meter, t_meter, t_s_meter, bool>;

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
            t_meter calcDeltaV(t_meter head, t_meter elevation, t_meter verticalSize) noexcept;

            t_s_meter_t calculateEFoldingConductance(FlowInputHor flow, t_meter folding_self, t_meter folding_neig);

            /**
             * @brief Calculates the horizontal flow between two nodes
             * @param flow a tuple of inputs about the aquifer
             * @return A weighted conductance value for the flow between two nodes
             * Calculates the harmonic mean conductance between two nodes.
             * $C = 2 \times EdgeLenght_1 \times \frac{ (TR_1 \times TR_2)}{(TR_1 \times EdgeLenght_1 + TR_2 \times EdgeLenght_2)}$
             */
            t_s_meter_t calculateHarmonicMeanConductance(FlowInputHor flow) noexcept;

            /**
             * Get the coefficients for storage and P components
             * @param steadyState
             * @param stepModifier
             * @param storageCapacity
             * @param P
             * @return HCOF
             */
            t_s_meter_t getHCOF(bool steadyState, quantity<Dimensionless> stepModifier,
                                t_s_meter storageCapacity, t_s_meter_t P) noexcept;

            /**
             * Calculates the vertical flow between two nodes
             * @param flow a tuple of inputs about the aquifer
             * @return the vertical conductance
             */
            t_s_meter_t calculateVerticalConductance(FlowInputVert flow) noexcept;

            /**
             * Criv = Krb/e*L*W
             * Krb/e = 2*Kh/W*(gamma/(1-G*gamma/Daq)
             * gamma = 1/(2*(1+1/pi*ln(2/(1-sqrt(e^(-pi*WpN)))))) where WpN is a normalized wetted perimeter which is W/Daq
             *
             * Krb - hydraulic conductivity of the river bed (we will just use that of the aquifer in that cell)
             * e - thickness of the river bed (this is totally unknown, but this falls out in the approximation that we will use)
             * L - the length of the river in the cell (intersect length...not sure you have this or not?)
             * W - width of the river in the cell
             * Kh - horizontal hydraulic conductivity of the aquifer (can use the same as Krb above as a start)
             * G - is the side length of the finite difference cell
             * Daq - is the average depth (thickness) of the aquifer in the cell (layer thickness)
             *
             */
           t_s_meter_t
           estimateConductance(t_vel K, t_meter length, t_meter width, t_meter Daq, t_meter G, t_meter depth) {
                assert(K.value() != 0 and length.value() != 0 and width.value() != 0 and Daq.value() != 0 and G.value() != 0 && "Inputs can't be 0!");
                const t_dim two = 2 * d;
                const t_dim four = 4 * d;
                t_s_meter_t criv;
                t_dim WpN = width / Daq;
                t_dim dpN = depth / Daq;
                t_dim ani = 0.1 * d;
                t_dim rho = sqrt(ani);
                t_dim Xi = (1 - sqrt(dpN)) * (1 - rho);
                assert(Xi.value() > 0 && "Depth of river and thickness of cell are not fit for this equation");
                t_dim Redc = 1 - 0.333 * Xi - 0.294 * pow(Xi, 2);
                t_meter Std_fd = two * Daq / rho;
                t_meter Delta_std = Std_fd - (two * Daq);
                t_meter B = width / two;
                //t_meter delta = G - (2 / rho * Daq + B);
                //t_meter delta = G - ((8*d) * Daq/rho + two*B);
                t_meter delta = G/four - B - two* Daq/rho;

                t_dim a1 = 1 * d;
                t_dim a2 = -1 * d;
                if (WpN < 1 and dpN < 0.2) {
                    a1 = 0.89 * d;
                    a2 = -2.43 * d;
                } else if (WpN < 1 and dpN < 0.5) {
                    a1 = 0.538 * d;
                    a2 = -0.387 * d;
                } else if (WpN < 3 and dpN < 0.2) {
                    a1 = 0.819 * d;
                    a2 = -1.34 * d;
                } else if (WpN < 3 and dpN < 0.2) {
                    a1 = 0.672 * d;
                    a2 = -0.542 * d;
                } else if (WpN < 3 and dpN < 0.09) {
                    a1 = 0.667 * d;
                    a2 = -0.33 * d;
                } else{
                    //default case
                    a1 = 0.89 * d;
                    a2 = -2.43 * d;
		}
                t_dim gamma_flat = 1 / (2 * (1 + 1 / pi<double> * log(2 / (1 - sqrt(exp(-pi<double> * WpN))))));
                t_dim gamma_flat_rect = gamma_flat * (1 + a1 * dpN + a2 * pow(dpN, 2));
                t_dim gamma_rect_iso = gamma_flat_rect / (1 + gamma_flat_rect * Delta_std / Daq);

                t_dim gamma = Redc * gamma_rect_iso / (1 + Redc * gamma_rect_iso * delta / Daq);

                criv = two * length * K * gamma;
		if(criv.value() <=0){
			std::cout << gamma.value() << "\n";
			std::cout << K.value() << "\n";
			std::cout << length.value() << "\n";
			std::cout << width.value() << "\n";
			std::cout << Daq.value() << "\n";
			std::cout << G.value() << "\n";
			std::cout << depth.value() << "\n";
            		std::exit(-1);
		}
        //        assert(criv.value() <= 0 &&
        //               "Conductance of riverbed is negative or zero - wrong inputs or eqation not fit for this scale");
		return criv;
            } 
        };

    }
}//ns
#endif
