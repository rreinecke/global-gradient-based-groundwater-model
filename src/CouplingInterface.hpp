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

#ifndef GLOBAL_FLOW_COUPLINGINTERFACE_HPP
#define GLOBAL_FLOW_COUPLINGINTERFACE_HPP

#include "Misc/Helpers.hpp"
#include "Model/Node.hpp"
#include "DataProcessing/DataReader.hpp"

namespace GlobalFlow {
    using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<Model::NodeInterface>>>;

    /**
     * This interface needs to be implemented to access datafields in the simulation class of G³M-f
     * How data is transferred from the coupled model to G³M-f is the responsibility of the implemented class
     */
    template<class Container>
    class CouplingInterface {
    public:
        virtual ~CouplingInterface() {}

        CouplingInterface(NodeVector nodeVector, DataReader *reader) : nodes(nodeVector), reader(reader) {}

        virtual void updateRechargeInGW(Container data, short month, int numberOfGridCells) {}

        virtual void updateNetAbstractionInGW(Container data, short month, int numberOfGridCells) {}

        virtual void updateRiversInGW(Container data, short month, int numberOfGridCells) {}

        virtual void updateGlobalWetlandsInGW(Container data, short month, int numberOfGridCells) {}

        virtual void updateLocalWetlandsInGW(Container data, short month, int numberOfGridCells) {}

        virtual void updateGlobalLakesInGW(Container data, short month, int numberOfGridCells) {}

        virtual void updateLakesInGW(Container data, short month, int numberOfGridCells) {}

        virtual void writeRiver(Container data, short month, int numberOfGridCells) {}

        virtual void writeGwFlow(Container data, int daysNotSimulated, short month, int numberOfGridCells) {}

        virtual void writeGlobalWetlands(Container data, short month, int numberOfGridCells) {}

        virtual void writeLocalWetlands(Container data, short month, int numberOfGridCells) {}

        virtual void writeGlobalLakes(Container data, short month, int numberOfGridCells) {}

        virtual void writeLakes(Container data, short month, int numberOfGridCells) {}

        virtual void getStorageData(Container data, short month, int numberOfGridCells) {}

        virtual void writeGwStorage(Container data, short month, int numberOfGridCells) {}

        virtual void writeFlowsIntoGW(Container gwr_loclak, Container gwr_locwet, Container gwr_glolac, Container gwr_res, Container gwr_glowet, int daysNotSimulated, short month, int numberOfGridCells){}

        virtual void deleteSWB(int cell, int numberOfGridCells, std::array<bool, 3> SWBType){}

        virtual void deleteGloLake(int cell, int numberOfGridCells){}

        virtual void activateDeactivateGloLakes(int cell, int numberOfGridCells, bool activate){}

        virtual void updateSwbConduct(int cell, int numberOfGridCells, std::array<double, 4> changeSwb){}

        virtual void saveSteadyStateConduct(int cell){}

        virtual void saveSteadyStateRiverDepth(int cell){}

        virtual void NoFiveMinuteCells(int cell, int numberOfGridCells, int *temp){}

        virtual void initMapping(int numberOfGridCells){}

    protected:
        NodeVector nodes;
        DataReader *reader;
    };
}
#endif
