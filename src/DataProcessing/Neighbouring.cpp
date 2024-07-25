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
                sumBoundaries++;
            }
                break;
            case Simulation::Options::GENERAL_HEAD_BOUNDARY: {
                nodes->at(nodeID)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, 0 * Model::si::meter,
                                                   boundaryConduct, 0 * Model::si::meter);
                sumBoundaries++;
            }
            default:
                break;
        }
    }
    return sumBoundaries;
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
                   std::unordered_map<large_num, std::unordered_map<int, large_num>> spatIDtoNodeID,
                   double resolution, large_num xRange, large_num yRange, bool isGlobal, int layers,
                   large_num numberOfNodesPerLayer, double boundaryConduct,
                   Simulation::Options::BoundaryCondition boundaryCondition) {
    // unrefined positions
    std::unordered_map<int, Model::NeighbourPosition> neigPositions;
    neigPositions[0] = Model::FRONT;
    neigPositions[1] = Model::BACK;
    neigPositions[2] = Model::RIGHT;
    neigPositions[3] = Model::LEFT;

    large_num nodeID_neig;
    large_num spatID;
    large_num spatID_neig;
    large_num nodeID;
    Model::NeighbourPosition neigPos;
    int sumBoundaries{0};

    for (large_num nodeIDTopLayer = 0; nodeIDTopLayer < numberOfNodesPerLayer; ++nodeIDTopLayer) {
        spatID = nodes->at(nodeIDTopLayer)->getSpatID();
        for (int layer = 0; layer < layers; ++layer) {
            nodeID = nodeIDTopLayer + (numberOfNodesPerLayer * layer);
            for (int neigPosID = 0; neigPosID < neigPositions.size(); ++neigPosID) {
                neigPos = neigPositions[neigPosID];
                // #### set neighbour of node
                spatID_neig = getNeighbourSpatID((int) spatID, neigPos, resolution, xRange, yRange, isGlobal);
                if (spatIDtoNodeID.contains(spatID_neig)) {
                    nodeID_neig = spatIDtoNodeID.at(spatID_neig).at(layer);
                    nodes->at(nodeID)->setNeighbour(nodeID_neig, neigPositions[neigPosID]);
                } else {
                    sumBoundaries = addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer, isGlobal,
                                                    sumBoundaries);
                }
            }
        }
    }
    LOG(debug) << "    Added " << sumBoundaries << " boundaries";
}


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

    LOG(debug) << "Building additional layers for " << numberOfLayers << " layers";

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
    bool isSteadyState;
    bool isDensityVariable;
    std::vector<Model::quantity<Model::Dimensionless>> delnus;
    std::vector<Model::quantity<Model::Dimensionless>> nusInZones;
    double effPorosity;
    double maxTipSlope;
    double maxToeSlope;
    double minDepthFactor;
    double slopeAdjFactor;
    Model::quantity<Model::Meter> vdfLock;
    int sourceZoneGHB;
    int sourceZoneRecharge;

    for (int layer = 0; layer < numberOfLayers - 1; ++layer) {
        //1) Add a Model::similar node in z direction for each layer
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
            isSteadyState = nodes->at(i)->getProperties().get<bool, Model::IsSteadyState>();
            isDensityVariable = nodes->at(i)->getProperties().get<bool, Model::IsDensityVariable>();
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
            vdfLock = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::VDFLock>();
            sourceZoneGHB = nodes->at(i)->getProperties().get<int, Model::SourceZoneGHB>();
            sourceZoneRecharge = nodes->at(i)->getProperties().get<int, Model::SourceZoneRecharge>();

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
                                                            isSteadyState,
                                                            isDensityVariable,
                                                            delnus,
                                                            nusInZones,
                                                            effPorosity,
                                                            maxTipSlope,
                                                            maxToeSlope,
                                                            minDepthFactor,
                                                            slopeAdjFactor,
                                                            vdfLock,
                                                            sourceZoneGHB,
                                                            sourceZoneRecharge));
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
    LOG(debug) << "Done: building additional layers";
};

}//ns
}