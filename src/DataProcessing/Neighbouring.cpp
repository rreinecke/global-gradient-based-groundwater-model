#include "Neighbouring.hpp"

namespace GlobalFlow {
    namespace DataProcessing {

/**
 * @brief function to add boundaries to nodes with no neighbours
 * @return total sum of boundaries added
 *
 */
int addBoundary(NodeVector const& nodes,
                 double boundaryConduct,
                 Simulation::Options::BoundaryCondition boundaryCondition,
                 large_num nodeID,
                 int layer, bool isGlobal, int sumBoundaries) {

    if (layer > 0 or nodes->at(nodeID)->hasGHB()) {
        return sumBoundaries;
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
    return sumBoundaries + 1;
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
                   std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<large_num, large_num>>> spatIDtoNodeIDs,
                   double resolution, large_num xRange, large_num yRange, bool isGlobal, int layers,
                   large_num numberOfNodesPerLayer, double boundaryConduct,
                   Simulation::Options::BoundaryCondition boundaryCondition) {
    // unrefined positions
    std::unordered_map<int, Model::NeighbourPosition> neigPositions;
    neigPositions[0] = Model::FRONT;
    neigPositions[1] = Model::BACK;
    neigPositions[2] = Model::RIGHT;
    neigPositions[3] = Model::LEFT;

    std::unordered_map<large_num, large_num> nodeIDs_neig;
    large_num spatID;
    large_num spatID_neig;
    large_num refID;
    large_num nodeID;
    Model::NeighbourPosition neigPos;
    int sumBoundaries{0};

    for (large_num nodeIDTopLayer = 0; nodeIDTopLayer < numberOfNodesPerLayer; ++nodeIDTopLayer) {
        spatID = nodes->at(nodeIDTopLayer)->getSpatID();
        refID = nodes->at(nodeIDTopLayer)->getRefID();
        for (int layer = 0; layer < layers; ++layer) {
            nodeID = nodeIDTopLayer + (numberOfNodesPerLayer * layer);
            large_num refinedInto = spatIDtoNodeIDs.at(spatID).at(layer).size();
            nodes->at(nodeID)->getProperties().set<large_num, Model::RefinedInto>(refinedInto);
            //LOG(debug) << "nodeID:" << nodeID << ", refinedInto: " << refinedInto;
            for (int neigPosID = 0; neigPosID < neigPositions.size(); ++neigPosID) {
                neigPos = neigPositions[neigPosID];
                if (refID == 0) { // #### set neighbour(s) of unrefined node
                    spatID_neig = getNeighbourSpatID((int) spatID, neigPos, resolution, xRange, yRange, isGlobal);
                    if (spatIDtoNodeIDs.contains(spatID_neig)) {
                        nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig).at(layer);
                        nodes->at(nodeID)->setNeighbours(nodeIDs_neig, neigPositions[neigPosID]);
                    } else {
                        sumBoundaries = addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer, isGlobal,
                                                    sumBoundaries);
                    }
                } else { // ####  set neighbour of refined node ####
                    sumBoundaries = setNeigOfRefinedNode(nodes, spatID, neigPos, resolution, xRange, yRange, isGlobal,
                                                         refID, nodeID, layer,spatIDtoNodeIDs, boundaryConduct,
                                                         boundaryCondition, sumBoundaries);
                }
            }
        }
    }
    LOG(debug) << "    Added " << sumBoundaries << " boundaries";
}


int setNeigOfRefinedNode(NodeVector nodes, large_num spatID, Model::NeighbourPosition neigPos, double resolution,
                          large_num xRange, large_num yRange, bool isGlobal, large_num refID, large_num nodeID, int layer,
                          std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<large_num, large_num>>> spatIDtoNodeIDs,
                      double boundaryConduct, Simulation::Options::BoundaryCondition boundaryCondition,
                      int sumBoundaries) {
    large_num refID_neig{};
    std::unordered_map<large_num, large_num> nodeIDs = spatIDtoNodeIDs.at(spatID).at(layer);
    // define map to set neighbour OUTSIDE refined node: neigPos >maps to> refID >maps to> refID_neig
    auto mapOutside = defineMapOutside(nodeIDs.size());

    // set neighbour OUTSIDE refined node (need to find spatID and find out whether that node is refined or not)
    try {
        refID_neig = mapOutside.at(neigPos).at(refID);
        setNeighbourOutsideRefinedNode(nodes, spatID, neigPos, resolution, xRange, yRange, isGlobal, nodeID, layer,
                                       spatIDtoNodeIDs, boundaryConduct, boundaryCondition, refID_neig, neigPos,
                                       sumBoundaries);
    } catch (const std::out_of_range &ex) {}

    // define map to set neighbour WITHIN refined node: neigPos >maps to> refID >maps to> refID_neig
    auto mapInside = defineMapInside(nodeIDs.size());

    // set neighbour WITHIN refined node (same spatID, thus we can just set the neighbour)
    try {
        refID_neig = mapInside.at(neigPos).at(refID);
        large_num nodeID_ref_neig = spatIDtoNodeIDs.at(spatID).at(layer).at(refID_neig); // layer = 0
        nodes->at(nodeID)->setNeighbour(nodeID_ref_neig, neigPos);
    } catch (const std::out_of_range &ex) {}
    return sumBoundaries;
}

/**
 * @brief Define map to set neighbour OUTSIDE refined node: neigPos >maps to> refID >maps to> refID_neig
 * @param refinedInto Number of refined nodes the original grid's node size (at spatID) is refined to
 * @return
 */
std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, large_num>>
defineMapOutside(large_num refinedInto) {
    std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, large_num>> mapOutside;
    if (refinedInto == 4) {
        /*
         *                  FRONT (neig)
         *                   3   4
         *                  /\  /\
         *                  ||  ||
         *             2 <=  1   2  => 1
         * RIGHT (neig)     (this)     RIGHT (neig)
         *             4 <=  3   4  => 3
         *                  ||  ||
         *                  \/  \/
         *                   1   2
         *                  BACK (neig)
         *
         */
        mapOutside = { { Model::FRONT, { {1, 3}, {2, 4} } },
                       { Model::BACK,  { {3, 1}, {4, 2} } },
                       { Model::RIGHT, { {2, 1}, {4, 3} } },
                       { Model::LEFT,  { {1, 2}, {3, 4} } } };
    } else if (refinedInto == 9) {
        /*
         *                   7  8  9  FRONT (neig)
         *                   /\ /\ /\
         *                   || || ||
         *              3 <=  1  2  3  => 1
         * LEFT (neig)  6 <=  4  5  6  => 2  RIGHT (neig)
         *                     (this)
         *              9 <=  7  8  9  => 3
         *                   || || ||
         *                   \/ \/ \/
         *                   1  2  3  BACK (neig)
         *
         */
        mapOutside = { { Model::FRONT, { {1, 7}, {2, 8}, {3, 9}} },
                       { Model::BACK,  { {7, 1}, {8, 2}, {9, 3} } },
                       { Model::RIGHT, { {3, 1}, {6, 4}, {9, 7}} },
                       { Model::LEFT,  { {1, 3}, {4, 6}, {7, 9}} } };
    } else if (refinedInto == 16) {
        /*
         *                   13 14 15 16 FRONT (neig)
         *                   /\ /\ /\ /\
         *                   || || || ||
         *              4 <=  1  2  3  4 =>  1
         * LEFT (neig)  8 <=  5  6  7  8 =>  5  RIGHT (neig)
         *                     (this)
         *             12 <=  9 10 11 12 =>  9
         *             16 <= 13 14 15 16 => 13
         *                   || || || ||
         *                   \/ \/ \/ \/
         *                    1  2  3  4 BACK (neig)
         *
         */
        mapOutside = { { Model::FRONT, { {1, 13}, {2, 14}, {3, 15}, {4, 16}  } },
                       { Model::BACK,  { {13, 1}, {14, 2}, {15, 3}, {16, 4}  } },
                       { Model::RIGHT, { {4, 1},  {8, 5},  {12, 9}, {16, 13} } },
                       { Model::LEFT,  { {1, 4},  {5, 8},  {9, 12}, {13, 16} } } };
    }
    return mapOutside;
}


/**
 * @brief Define map to set neighbour OUTSIDE refined node: neigPos >maps to> refID >maps to> refID_neig
 * @param refinedInto Number of refined nodes the original grid's node size (at spatID) is refined to
 * @return
 */
std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, large_num>>
defineMapInside(large_num refinedInto) {
    std::unordered_map<Model::NeighbourPosition, std::unordered_map<large_num, large_num>> mapInside;
    if (refinedInto == 4) {
        /*
         *   1 - 2
         *   |   |
         *   3 - 4
         */
        mapInside = { {Model::FRONT, { {3, 1}, {4, 2} } },
                      {Model::BACK,  { {1, 3}, {2, 4} } },
                      {Model::RIGHT, { {1, 2}, {3, 4} } },
                      {Model::LEFT,  { {2, 1}, {4, 3} } } };
    } else if (refinedInto == 9) {
        /*
         *   1 - 2 - 3
         *   |   |   |
         *   4 - 5 - 6
         *   |   |   |
         *   7 - 8 - 9
         */
        mapInside = { {Model::FRONT, { {4, 1}, {5, 2}, {6, 3}, {7, 4}, {8, 5}, {9, 6} } },
                      {Model::BACK,  { {1, 4}, {2, 5}, {3, 6}, {4, 7}, {5, 8}, {6, 9} } },
                      {Model::RIGHT, { {1, 2}, {2, 3}, {4, 5}, {5, 6}, {7, 8}, {8, 9} } },
                      {Model::LEFT,  { {2, 1}, {3, 2}, {5, 4}, {6, 5}, {8, 7}, {9, 8} } } };
    } else if (refinedInto == 16) {
        /*
         *   1 -  2 -  3 -  4
         *   |    |    |    |
         *   5 -  6 -  7 -  8
         *   |    |    |    |
         *   9 - 10 - 11 - 12
         *   |    |    |    |
         *  13 - 14 - 15 - 16
         */
        mapInside = { {Model::FRONT, { {5, 1},  {6, 2},   {7, 3},   {8, 4},   {9, 5},   {10, 6},
                                       {11, 7}, {12, 8},  {13, 9},  {14, 10}, {15, 11}, {16, 12}} },
                      {Model::BACK,  { {1, 5},  {2, 6},   {3, 7},   {4, 8},   {5, 9},   {6, 10},
                                       {7, 11}, {8, 12},  {9, 13},  {10, 14}, {11, 15}, {12, 16}} },
                      {Model::RIGHT, { {1, 2},  {2, 3},   {3, 4},   {5, 6},   {6, 7},   {7, 8},
                                       {9, 10}, {10, 11}, {11, 12}, {13, 14}, {14, 15}, {15, 16} } },
                      {Model::LEFT,  { {2, 1},  {3, 2},   {4, 3},   {6, 5},   {7, 6},   {8, 7},
                                       {10, 9}, {11, 10}, {12, 11}, {14, 13}, {15, 14}, {16, 15}} } };
    }
    return mapInside;
}

int setNeighbourOutsideRefinedNode(NodeVector nodes, large_num spatID, Model::NeighbourPosition neigPos, double resolution,
                                    large_num xRange, large_num yRange, bool isGlobal, large_num nodeID, int layer,
                                    std::unordered_map<large_num, std::unordered_map<int, std::unordered_map<large_num, large_num>>> spatIDtoNodeIDs,
                                    double boundaryConduct, Simulation::Options::BoundaryCondition boundaryCondition,
                                    large_num ref_id_neig, Model::NeighbourPosition neighbourPosition,
                                    int sumBoundaries) {
    int spatID_neig = getNeighbourSpatID((int) spatID, neigPos, resolution, xRange, yRange, isGlobal);
    if (spatIDtoNodeIDs.contains(spatID_neig)) {
        std::unordered_map<large_num, large_num> nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig).at(layer); // layer = 0
        if (nodeIDs_neig.size() == 1) {
            nodes->at(nodeID)->setNeighbour(nodeIDs_neig.at(0), neighbourPosition); // refID = 0
        } else {
            nodes->at(nodeID)->setNeighbour(nodeIDs_neig.at(ref_id_neig), neighbourPosition);
        }
    } else {
        sumBoundaries = addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer, isGlobal, sumBoundaries); // layer = 0
    }
    return sumBoundaries;
};



/**
 * @brief calculates all neighbours based on a spatial ID
 * Assumes no landmask and calculates all possible neighbours
 * @param id
 * @param res
 * @return
 */
int getNeighbourSpatID(int spatID, Model::NeighbourPosition neigPos, double res, large_num xRange, large_num yRange,
                       bool isGlobal) {
    int rowLength{(int) std::round((double) xRange / res)};
    int colLength{(int) std::round((double) yRange / res)};
    /*if (rowLength >= 2) {
        assert(rowLength % 2 == 0 && "resolution is impossible");
    }*/
    switch (neigPos) {
        case Model::FRONT:
            if (spatID >= rowLength) { return spatID - rowLength; } else { return -1; }
            break;
        case Model::BACK:
            if (spatID <= colLength * rowLength - rowLength) { return spatID + rowLength; } else { return -1; }
            break;
        case Model::RIGHT:
            if ((spatID + 1) % rowLength == 0) {
                if (isGlobal) { return spatID - rowLength + 1; } else { return -1; }
            } else { return spatID + 1; }
            break;
        case Model::LEFT:
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

    //LOG(debug) << "Building additional layers with node count: " << nodesPerLayer << " for " << numberOfLayers
    //           << " layers";

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
    large_num refID;
    large_num maxRefinement;
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
            K = conductances[layer + 1] * Model::si::meter / Model::day; // or: nodes->at(i)->getK__pure();
            head = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::Head>();
            aquiferDepth = aquifer_thickness[layer + 1];
            anisotropy = anisotropies[layer + 1]; // or nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>, Model::Anisotropy>().value();
            specificYield = nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                            Model::SpecificYield>().value();
            specificStorage = nodes->at(i)->getProperties().get<Model::quantity<Model::perUnit>,
                            Model::SpecificStorage>().value();
            useEfolding = nodes->at(i)->getProperties().get<bool, Model::UseEfolding>();
            refID = nodes->at(i)->getProperties().get<large_num, Model::RefID>();
            maxRefinement = nodes->at(i)->getProperties().get<large_num, Model::MaxRefinement>();
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
                                                            maxRefinement,
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
    //LOG(debug) << "Last nodeID was " << id << " with max ID (with non static nodes) "
    //           << nodesPerLayer * numberOfLayers;
};

}//ns
}