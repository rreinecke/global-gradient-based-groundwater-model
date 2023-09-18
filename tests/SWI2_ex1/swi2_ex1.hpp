/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef GLOBAL_FLOW__MAIN_HPP
#define GLOBAL_FLOW__MAIN_HPP

#include "../../src/GW_Interface.hpp"
#include "../../src/Simulation/Stepper.hpp"
#include "../../src/Simulation/Simulation.hpp"
#include "../../src/Simulation/Options.hpp"
#include "../../src/DataProcessing/DataOutput/OutputManager.hpp"
#include "set"
#include <fstream>
#include "../../src/Misc/colors.hpp"
#include "SWI2_ex1_DataReader.hpp"

namespace GlobalFlow {

    class NotCoupled;

    class StandaloneRunner : GW_Interface<NotCoupled> {
        public:

            StandaloneRunner();

            void
            loadSettings();

            void
            setupSimulation();

            void
            writeNodeInfosToCSV();

            void
            simulate();

            void
            getResults();

            void
            writeData();

        private:
            Solver::Equation *_eq;
            Simulation::Options op;
            Simulation::Simulation sim;
            DataProcessing::SWI2_ex1_DataReader *reader;
    };

}//ns
#endif //GROUNDWATER_SOLVER_MAIN_HPP

