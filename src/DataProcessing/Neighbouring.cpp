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
                 int layer) {

    if (layer > 0) {
        return;
    }

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
                   int numberOfLayers,
                   double boundaryConduct,
                   Simulation::Options::BoundaryCondition boundaryCondition) {
    large_num nodes_per_layer = nodes->size() / numberOfLayers;

    auto lu = setNeighbourPositions();

    large_num nodeID;
    std::unordered_map<int, large_num> nodeIDs_neig;
    large_num spatID;
    int spatID_neig;
    int refID;

    for (int layer = 0; layer < numberOfLayers; ++layer) {
        for (int i = 0; i < nodes_per_layer; ++i) {
            nodeID = i + (nodes_per_layer * layer);
            spatID = nodes->at(nodeID)->getSpatID();
            refID = nodes->at(nodeID)->getRefID();
            for (int j = 0; j < lu.size(); ++j) {
                if (refID == 0) { // #### set neighbour(s) of unrefined node
                    spatID_neig = getNeighbourSpatID((int) spatID, j, resolution, lonRange, latRange, isGlobal);
                    if (spatIDtoNodeIDs.contains(spatID_neig)) {
                        nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig).at(layer);
                        if (nodeIDs_neig.empty()) {
                        } else if (nodeIDs_neig.size() == 1) {
                            nodes->at(nodeID)->setNeighbour(nodeIDs_neig[0], lu[j]);
                        } else {
                            nodes->at(nodeID)->setNeighbours(nodeIDs_neig, lu[j]);
                        }
                    } else {
                        addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer);
                    }
                } else { // ####  set neighbour of refined node ####
                    setNeigOfRefinedNode(nodes, spatID, j, resolution, lonRange, latRange, isGlobal, layer, refID, nodeID,
                                     spatIDtoNodeIDs, boundaryConduct, boundaryCondition);
                }
            }
        }
    }

}


void setNeigOfRefinedNode(NodeVector nodes, large_num spatID, int j, double resolution,
                      int lonRange, int latRange, bool isGlobal, int layer, int refID, large_num nodeID,
                      std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<int, large_num>>> spatIDtoNodeIDs,
                      double boundaryConduct, Simulation::Options::BoundaryCondition boundaryCondition) {

    auto neighbourPositions = setNeighbourPositions();

    std::unordered_map<int, std::vector<int>> mapping;
    mapping[0] = {1, 2, 3, 4}; // FRONT
    mapping[1] = {3, 4, 1, 2}; // BACK
    mapping[2] = {2, 4, 1, 3}; // RIGHT
    mapping[3] = {1, 3, 2, 4}; // LEFT

    auto neighbourOutsideRefinedNode = [spatID, j, resolution, lonRange, latRange, isGlobal, spatIDtoNodeIDs, layer,
                                       nodes, boundaryConduct, boundaryCondition, nodeID]
                                               (int index, Model::NeighbourPosition neighbourPosition) {
        int spatID_neig = getNeighbourSpatID((int) spatID, j, resolution, lonRange, latRange, isGlobal);
        if (spatIDtoNodeIDs.contains(spatID_neig)) {
            std::unordered_map<int, large_num> nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig).at(layer);
            if (nodeIDs_neig.size() == 1) {
                nodes->at(nodeID)->setNeighbour(nodeIDs_neig.at(0), neighbourPosition);
            } else {
                nodes->at(nodeID)->setNeighbour(nodeIDs_neig.at(index), neighbourPosition);
            }
        } else { addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer); }
    };

    if (refID == mapping.at(j)[0]) {
        neighbourOutsideRefinedNode(mapping.at(j)[2], neighbourPositions[j]);
    } else if (refID == mapping.at(j)[1]) {
        neighbourOutsideRefinedNode(mapping.at(j)[3], neighbourPositions[j]);
    }
    if (refID == mapping.at(j)[2] or refID == mapping.at(j)[3]) { // set neighbour within refined node (same spatID)
        large_num nodeID_ref_neig = spatIDtoNodeIDs.at(spatID).at(layer).at(refID);
        nodes->at(nodeID)->setNeighbour(nodeID_ref_neig, neighbourPositions[j]);
    }
}

std::unordered_map<int, Model::NeighbourPosition> setNeighbourPositions() {
    // unrefined positions
    std::unordered_map<int, Model::NeighbourPosition> neighbourPositions;
    neighbourPositions[0] = Model::FRONT; // formerly NORTH
    neighbourPositions[1] = Model::BACK; // formerly SOUTH
    neighbourPositions[2] = Model::RIGHT; // formerly EAST
    neighbourPositions[3] = Model::LEFT; // formerly WEST
    /*if (refined) {
        // refined positions
        neighbourPositions[4] = Model::FRONTLEFT;
        neighbourPositions[5] = Model::FRONTRIGHT;
        neighbourPositions[6] = Model::BACKLEFT;
        neighbourPositions[7] = Model::BACKRIGHT;
        neighbourPositions[8] = Model::RIGHTFRONT;
        neighbourPositions[9] = Model::RIGHTBACK;
        neighbourPositions[10] = Model::LEFTFRONT;
        neighbourPositions[11] = Model::LEFTBACK;
    }*/
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
            assert(rowLength % 2 == 0 && "resolution is impossible");
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
 * @brief
 * @param from is position in vector of top layer node
 * @param to is position in vector of node that receive neighbouring information
 */
        void copyNeighbours(large_num from, large_num to, NodeVector nodes, large_num layer_shift) {
            auto neighbours = nodes->at(from)->getListOfNeighbours();
            for (const auto &n: neighbours) {
                if (n.first == Model::DOWN or n.first == Model::TOP) {
                    continue;
                }
                nodes->at(to)->setNeighbour(n.second + layer_shift, n.first);
            }
        }

/**
 * @brief Copies cardinal points of top layer to all bottom layers
 * @param nodes
 * @param layers
 */
        void copyNeighboursToBottomLayers(NodeVector nodes, int numberOfLayers) {
            assert(numberOfLayers && "0 layers does not make sense");
            if (numberOfLayers == 1) {
                return;
            }
            size_t nodesPerLayer = nodes->size() / numberOfLayers;

            for (large_num from_nodeID = 0; from_nodeID < nodesPerLayer; ++from_nodeID) {
                for (int to_layer = 1; to_layer < numberOfLayers; ++to_layer) {
                    large_num to_nodeID = from_nodeID + (nodesPerLayer * to_layer);

                    copyNeighbours(from_nodeID, to_nodeID, nodes, nodesPerLayer * to_layer);
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

            for (int j = 0; j < numberOfLayers - 1; ++j) {
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
                    K = conductances[j + 1] * Model::si::meter / Model::day;
                    head = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::Head>();
                    aquiferDepth = aquifer_thickness[j + 1];
                    anisotropy = anisotropies[j + 1];
                    specificYield =
                            nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                                    Model::SpecificYield>().value();
                    specificStorage =
                            nodes->at(i)->getProperties().get<Model::quantity<Model::perUnit>, Model::SpecificStorage>
                                    ().value();
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
                                                                    confined[j + 1],
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
                        nodes->at(id)->getProperties().set<int, Model::Layer>(j + 1);
                        nodes->at(id)->getProperties().set<Model::quantity<Model::Meter>, Model::Elevation>(
                                nodes->at(id)->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>()
                                - (aquiferDepth * Model::si::meter));
                    }
                    //2) Neighbouring for top and bottom

                    if (j > 0) {
                        //Layer above is not top layer
                        nodes->at(id)->setNeighbour(i + (j * nodesPerLayer), Model::TOP);
                        nodes->at(i + (j * nodesPerLayer))->setNeighbour(id, Model::DOWN);
                    } else {
                        //Layer above is top layer
                        nodes->at(id)->setNeighbour(i, Model::TOP);
                        nodes->at(i)->setNeighbour(id, Model::DOWN);
                    }

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