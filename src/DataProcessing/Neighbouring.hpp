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

#ifndef GLOBAL_FLOW_NEIGHBOURING_HPP
#define GLOBAL_FLOW_NEIGHBOURING_HPP

#include "DataReader.hpp"
#include "../Simulation/Options.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow
{
namespace DataProcessing
{
template<class T>
using Matrix = std::vector<std::vector<T>>;
using large_num = unsigned long int;
using n_array = std::array<int,4>;

/**
 * Modflow like grid file
 * @param nodes
 * @param grid
 * @param layers
 */
void buildByGrid(NodeVector nodes, Matrix<int> grid, int nodesPerLayer, int layers);
/**
* Builds a map of neighbouring nodes based spatial Id's and resolution
* Missing neighbours or empty spaces lead to adding of a General Head Boundary Flow addition
*/
void buildBySpatID(NodeVector nodes, std::unordered_map<int, std::vector<int>> spatIDtoNodeIDs, int resolution,
                   int numberOfLayers, double oceanCoduct, Simulation::Options::BoundaryCondition boundaryCondition);

void copyNeighboursToBottomLayers(NodeVector nodes, int layers);

n_array getNeighbourByNodeID(int nodeID, int res);


/**
* Builds a map of neighbouring nodes based on x and y coordinates
* Missing neighbours or empty spaces lead to adding of a General Head Boundary Flow addition
*/
int buildNeighbourMap(NodeVector nodes, int numberOfTOPNodes, int layers, double ghbConduct, Simulation::Options::BoundaryCondition boundaryCondition);

/**
* Extend the current nodes with an arbitrary number of layers
* @param layers: Number of total layers
* Thus 1 is the default setting, 0 is undefined
*/
void
buildBottomLayers(NodeVector nodes,
                  int layers,
                  std::vector<bool> conf,
                  std::vector<int> aquifer_depth,
                  std::vector<double> conductances,
                  std::vector<double> anisotropies);
}
}

#endif
