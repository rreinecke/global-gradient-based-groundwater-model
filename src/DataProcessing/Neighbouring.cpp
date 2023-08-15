#include "Neighbouring.hpp"

namespace GlobalFlow {
    namespace DataProcessing {

/**
 * For regular grids (modflow like)
 * @param nodes
 * @param grid
 * @param nodesPerLayer
 * @param layers
 */
        void buildByGrid(NodeVector nodes, Matrix<int> grid, large_num nodesPerLayer, int layers) {
            //id->row,col
            int rows = grid[0].size();
            int cols = grid.size();
            LOG(debug) << "cols: " << cols << ", rows: " << rows << ", nodes per layer: " << nodesPerLayer
                       << ", layers: " << layers << std::endl;

            auto check = [grid](int i, int j) {
                try {
                    int p = grid.at(i).at(j);
                    return true;
                } catch (std::out_of_range e) {
                    return false;
                }
            };

            for (int layer = 0; layer < layers; ++layer) {
                int l_mult = layer * nodesPerLayer;
                //id->row,col
                for (int i = 0; i < cols; ++i) {
                    for (int j = 0; j < rows; ++j) {
                        int id = grid[i][j];
                        if (check(i + 1, j)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i + 1][j] + l_mult, Model::RIGHT);
                        }
                        if (check(i - 1, j)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i - 1][j] + l_mult, Model::LEFT);
                        }
                        if (check(i, j - 1)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i][j - 1] + l_mult, Model::FRONT);
                        }
                        if (check(i, j + 1)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i][j + 1] + l_mult, Model::BACK);
                        }
                    }
                }
            }
        }



/**
 * @brief function to add boundaries to nodes with no neighbours
 *
 */
void addBoundary(NodeVector const& nodes,
                 double boundaryConduct,
                 Simulation::Options::BoundaryCondition boundaryCondition,
                 large_num nodeID,
                 int layer, bool isGlobal) {

    if (layer > 0) {
        return;
    }

    if (isGlobal){
        switch (boundaryCondition) {
            case Simulation::Options::GENERAL_HEAD_NEIGHBOUR: {
                auto head = nodes->at(nodeID)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
                nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, head, boundaryConduct, head);
            }
                break;
            case Simulation::Options::GENERAL_HEAD_BOUNDARY: {
                nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, 0 * Model::si::meter,
                                                   boundaryConduct, 0 * Model::si::meter);
            }
            default:
                break;
        }
    }
}

/**
 * @brief
 * @param nodes
 * @param id_mapping
 * @param resolution
 * @param layers
 * @param boundaryConduct
 * @param boundaryCondition
 */
void buildBySpatID(NodeVector nodes,
                   std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<int, large_num>>> spatIDtoNodeIDs,
                   double resolution,
                   int lonRange, int latRange, bool isGlobal,
                   int layers, large_num numberOfNodesPerLayer,
                   double boundaryConduct,
                   Simulation::Options::BoundaryCondition boundaryCondition) {
    auto lu = setNeighbourPositions();

    std::unordered_map<int, large_num> nodeIDs_neig;
    large_num spatID;
    int spatID_neig;
    int refID;
    large_num nodeID;

    for (large_num nodeIDTopLayer = 0; nodeIDTopLayer < numberOfNodesPerLayer; ++nodeIDTopLayer) {
        spatID = nodes->at(nodeIDTopLayer)->getSpatID();
        refID = nodes->at(nodeIDTopLayer)->getRefID();
        for (int layer = 0; layer < layers; ++layer) {
            nodeID = nodeIDTopLayer + (numberOfNodesPerLayer * layer);
            for (int j = 0; j < lu.size(); ++j) {
                if (refID == 0) { // #### set neighbour(s) of unrefined node
                    spatID_neig = getNeighbourSpatID((int) spatID, j, resolution, lonRange, latRange, isGlobal);
                    if (spatIDtoNodeIDs.contains(spatID_neig)) {
                        nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig).at(layer);
                        nodes->at(nodeID)->setNeighbours(nodeIDs_neig, lu[j]);
                    } else {
                        addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer, isGlobal);
                    }
                } else { // ####  set neighbour of refined node ####
                    setNeigOfRefinedNode(nodes, spatID, j, resolution, lonRange, latRange, isGlobal, refID, nodeID,
                                         layer,spatIDtoNodeIDs, boundaryConduct, boundaryCondition);
                }
            }
        }
    }

}


void setNeigOfRefinedNode(NodeVector nodes, large_num spatID, int j, double resolution,
                      int lonRange, int latRange, bool isGlobal, int refID, large_num nodeID, int layer,
                      std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<int, large_num>>> spatIDtoNodeIDs,
                      double boundaryConduct, Simulation::Options::BoundaryCondition boundaryCondition) {

    auto neighbourPositions = setNeighbourPositions();

    std::unordered_map<int, std::vector<int>> mapOutside; // mapping of neighbour position to refIDs
    mapOutside[0] = {1, 2, 3, 4}; // neighbour is at FRONT -> if refID of this node is 1 or 2: neighbour outside refined node
    mapOutside[1] = {3, 4, 1, 2}; // BACK
    mapOutside[2] = {2, 4, 1, 3}; // RIGHT
    mapOutside[3] = {1, 3, 2, 4}; // LEFT

    auto neighbourOutsideRefinedNode = [spatID, j, resolution, lonRange, latRange, isGlobal, layer, spatIDtoNodeIDs,
                                       nodes, boundaryConduct, boundaryCondition, nodeID]
                                               (int ref_id, Model::NeighbourPosition neighbourPosition) {
        int spatID_neig = getNeighbourSpatID((int) spatID, j, resolution, lonRange, latRange, isGlobal);
        if (spatIDtoNodeIDs.contains(spatID_neig)) {
            std::unordered_map<int, large_num> nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig).at(layer); // layer = 0
            if (nodeIDs_neig.size() == 1) {
                nodes->at(nodeID)->setNeighbour(nodeIDs_neig.at(0), neighbourPosition); // refID = 0
            } else {
                nodes->at(nodeID)->setNeighbour(nodeIDs_neig.at(ref_id), neighbourPosition);
            }
        } else {
                addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, 0, isGlobal); // layer = 0
        }
    };

    // set neighbour outside refined node (need to find spatID and find out whether that node is refined or not)
    if (refID == mapOutside.at(j)[0]) {
        neighbourOutsideRefinedNode(mapOutside.at(j)[2], neighbourPositions[j]);
    } else if (refID == mapOutside.at(j)[1]) {
        neighbourOutsideRefinedNode(mapOutside.at(j)[3], neighbourPositions[j]);
    }

    // set neighbour within refined node (same spatID, thus we can just set the neighbour)
    std::unordered_map<int, std::vector<int>> mapInside; // mapping of neighbour position to refIDs
    mapInside[0] = {3, 4, 1, 2}; // neighbour is at FRONT -> if refID of this node is 3 or 4: neighbour inside refined node
    mapInside[1] = {1, 2, 3, 4}; // BACK
    mapInside[2] = {1, 3, 2, 4}; // RIGHT
    mapInside[3] = {2, 4, 1, 3}; // LEFT

    if (refID == mapInside.at(j)[0]) {
        large_num nodeID_ref_neig = spatIDtoNodeIDs.at(spatID).at(layer).at(mapInside.at(j)[2]); // layer = 0
        nodes->at(nodeID)->setNeighbour(nodeID_ref_neig, neighbourPositions[j]);
    } else if (refID == mapInside.at(j)[1]){
        large_num nodeID_ref_neig = spatIDtoNodeIDs.at(spatID).at(layer).at(mapInside.at(j)[3]); // layer = 0
        nodes->at(nodeID)->setNeighbour(nodeID_ref_neig, neighbourPositions[j]);
    }
}

std::unordered_map<int, Model::NeighbourPosition> setNeighbourPositions() {
    // unrefined positions
    std::unordered_map<int, Model::NeighbourPosition> neighbourPositions;
    neighbourPositions[0] = Model::FRONT;
    neighbourPositions[1] = Model::BACK;
    neighbourPositions[2] = Model::RIGHT;
    neighbourPositions[3] = Model::LEFT;
    return neighbourPositions;
}

/**
 * @brief calculates all neighbours based on a spatial ID
 * Assumes no landmask and calculates all possible neighbours
 * @param id
 * @param res
 * @return
 */
        int getNeighbourSpatID(int spatID, int j, double res, int lonRange, int latRange, bool isGlobal) {
            int rowLength{(int) std::round(lonRange / res)};
            int colLength{(int) std::round(latRange / res)};
            if (rowLength >= 2) {
                assert(rowLength % 2 == 0 && "resolution is impossible");
            }
            switch (j) {
                case 0: //FRONT
                    if (spatID >= rowLength) { return spatID - rowLength; } else { return -1; }
                    break;
                case 1: //BACK
                    if (spatID <= colLength * rowLength - rowLength) { return spatID + rowLength; } else { return -1; }
                    break;
                case 2: //RIGHT
                    if ((spatID + 1) % rowLength == 0) {
                        if (isGlobal) { return spatID - rowLength + 1; } else { return -1; }
                    } else { return spatID + 1; }
                    break;
                case 3: //LEFT
                    if (spatID % rowLength == 0) {
                        if (isGlobal) { return spatID + rowLength - 1; } else { return -1; }
                    } else { return spatID - 1; }
                    break;
                default:
                    return -1;
            }
        }

/**
 * @brief Copies neighbours of top layer to all bottom layers
 * @param nodes
 * @param numberOfLayers
 */
        void copyNeighboursToBottomLayers(NodeVector nodes, int numberOfLayers) {
            assert(numberOfLayers && "0 layers does not make sense");
            if (numberOfLayers == 1) { return; }
            large_num nodesPerLayer = nodes->size() / numberOfLayers;

            // iterate through nodeIDs of top layer
            for (large_num from_nodeID = 0; from_nodeID < nodesPerLayer; ++from_nodeID) {
                // get neighbours of node in top layer
                auto neighbours = nodes->at(from_nodeID)->getListOfNeighbours();

                // iterate through bottom layers
                for (int to_layer = 1; to_layer < numberOfLayers; ++to_layer) {
                    // calculate nodeIDs in bottom layer
                    large_num to_nodeID = from_nodeID + (nodesPerLayer * to_layer);

                    // iterate through neighbours of top layer
                    for (const auto &n: neighbours) {
                        if (n.first == Model::DOWN or n.first == Model::TOP) {
                            continue; // ignore TOP and DOWN neighbours
                        }
                        // copy neighbour
                        large_num newNeighbourNodeID = n.second + (nodesPerLayer * to_layer);
                        nodes->at(to_nodeID)->setNeighbour(newNeighbourNodeID, n.first);
                    }
                }
            }
        }

/**
 * Adds additional layers and connects nodes
 * @param nodes
 * @param layers
 * @param confined
 * @param aquifer_thickness
 */
        void buildBottomLayers(NodeVector nodes,
                               int numberOfLayers,
                               std::vector<bool> confined,
                               std::vector<int> aquifer_thickness,
                               std::vector<double> conductances,
                               std::vector<double> anisotropies) {
            assert(numberOfLayers && "AsModel::signing 0 layers does not make any sense");
            if (numberOfLayers == 1) {
                return;
            }

            size_t nodesPerLayer = nodes->size();
            nodes->reserve(numberOfLayers * nodesPerLayer);

            LOG(debug) << "Building additional layers with node count: " << nodesPerLayer << " for " << numberOfLayers
                       << " layers";

            size_t id = nodesPerLayer;
            large_num spatID;
            double lat, lon;
            Model::quantity<Model::SquareMeter> area;
            Model::quantity<Model::Meter> edgeLengthLeftRight;
            Model::quantity<Model::Meter> edgeLengthFrontBack;
            Model::quantity<Model::Velocity> K;
            Model::quantity<Model::Meter> head;
            double aquiferDepth;
            double anisotropy;
            double specificYield;
            double specificStorage;
            bool useEfolding;
            int refID;
            bool densityVariable;
            std::vector<Model::quantity<Model::Dimensionless>> delnus;
            std::vector<Model::quantity<Model::Dimensionless>> nusInZones;
            double effPorosity;
            double maxTipSlope;
            double maxToeSlope;
            double minDepthFactor;
            double slopeAdjFactor;
            Model::quantity<Model::Meter> vdfLock;

            for (int layer = 0; layer < numberOfLayers - 1; ++layer) {
                //1) Add a Model::similar node in z direction for each layer
                //TODO Parallell?
                for (int i = 0; i < nodesPerLayer; ++i) {
                    //for each node in top layer
                    spatID = nodes->at(i)->getProperties().get<large_num, Model::SpatID>();
                    lat = nodes->at(i)->getProperties().get<double, Model::Lat>();
                    lon = nodes->at(i)->getProperties().get<double, Model::Lon>();
                    area = nodes->at(i)->getProperties().get<Model::quantity<Model::SquareMeter>, Model::Area>();
                    edgeLengthLeftRight = nodes->at(
                            i)->getProperties().get<Model::quantity<Model::Meter>, Model::EdgeLengthLeftRight>();
                    edgeLengthFrontBack = nodes->at(
                            i)->getProperties().get<Model::quantity<Model::Meter>, Model::EdgeLengthFrontBack>();
                    K = conductances[layer + 1] * Model::si::meter / Model::day;
                    head = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::Head>();
                    aquiferDepth = aquifer_thickness[layer + 1];
                    anisotropy = anisotropies[layer + 1];
                    specificYield = nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                                    Model::SpecificYield>().value();
                    specificStorage = nodes->at(i)->getProperties().get<Model::quantity<Model::perUnit>,
                                    Model::SpecificStorage>().value();
                    useEfolding = nodes->at(i)->getProperties().get<bool, Model::UseEfolding>();
                    refID = nodes->at(i)->getProperties().get<int, Model::RefID>();
                    densityVariable = nodes->at(i)->getProperties().get<bool, Model::DensityVariable>();
                    delnus = nodes->at(i)->getProperties().
                            get<std::vector<Model::quantity<Model::Dimensionless>>, Model::Delnus>();
                    nusInZones = nodes->at(i)->getProperties().
                            get<std::vector<Model::quantity<Model::Dimensionless>>, Model::NusInZones>();
                    effPorosity = nodes->at(i)->getProperties().
                            get<Model::quantity<Model::Dimensionless>, Model::EffectivePorosity>();
                    maxTipSlope = nodes->at(i)->getProperties().
                            get<Model::quantity<Model::Dimensionless>, Model::MaxTipSlope>();
                    maxToeSlope = nodes->at(i)->getProperties().
                            get<Model::quantity<Model::Dimensionless>, Model::MaxToeSlope>();
                    minDepthFactor = nodes->at(i)->getProperties().
                            get<Model::quantity<Model::Dimensionless>, Model::MinDepthFactor>();
                    slopeAdjFactor = nodes->at(i)->getProperties().
                            get<Model::quantity<Model::Dimensionless>, Model::SlopeAdjFactor>();
                    vdfLock = nodes->at(i)->getProperties().
                            get<Model::quantity<Model::Meter>, Model::VDFLock>();

                    if (nodes->at(i)->isStaticNode()) {
                        //is taken care of by neighbouring algorithm
                        continue;
                    } else {
                        if (id > nodesPerLayer * numberOfLayers) {
                            LOG(critical) << "This is not possible!";
                            exit(9);
                        }
                        nodes->emplace_back(new Model::StandardNode(nodes, lat, lon, area, edgeLengthLeftRight,
                                                                    edgeLengthFrontBack,
                                                                    spatID,
                                                                    id,
                                                                    K,
                                                                    head,
                                                                    aquiferDepth,
                                                                    anisotropy,
                                                                    specificYield,
                                                                    specificStorage,
                                                                    useEfolding,
                                                                    confined[layer + 1],
                                                                    refID,
                                                                    densityVariable,
                                                                    delnus,
                                                                    nusInZones,
                                                                    effPorosity,
                                                                    maxTipSlope,
                                                                    maxToeSlope,
                                                                    minDepthFactor,
                                                                    slopeAdjFactor,
                                                                    vdfLock));
                        nodes->at(id)->getProperties().set<int, Model::Layer>(layer + 1);
                        nodes->at(id)->getProperties().set<Model::quantity<Model::Meter>, Model::Elevation>(
                                nodes->at(id)->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>()
                                - (aquiferDepth * Model::si::meter));
                    }
                    //2) Neighbouring for top and bottom
                    nodes->at(id)->setNeighbour(i + (layer * nodesPerLayer), Model::TOP);
                    nodes->at(i + (layer * nodesPerLayer))->setNeighbour(id, Model::DOWN);

                    id++;
                    if (id > (nodesPerLayer * numberOfLayers) - 1) {
                        break;
                    }
                }
            }
            LOG(debug) << "Last nodeID was " << id << " with max ID (with non static nodes) "
                       << nodesPerLayer * numberOfLayers;
        };

    }//ns
}