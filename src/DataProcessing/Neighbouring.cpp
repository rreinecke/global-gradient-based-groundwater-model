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
void buildByGrid(NodeVector nodes, Matrix<int> grid, int nodesPerLayer, int layers) {
    //id->row,col
    int rows = grid[0].size();
    int cols = grid.size();
    LOG(debug) << "cols: " << cols << ", rows: " << rows << ", nodes per layer: " << nodesPerLayer << ", layers: " << layers << std::endl;

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
 * @brief
 * @param nodes
 * @param id_mapping
 * @param resolution
 * @param layers
 * @param boundaryConduct
 * @param boundaryCondition
 */
void buildByLatLon(NodeVector nodes, int resolution, int numberOfLayers, double boundaryConduct,
                   Simulation::Options::BoundaryCondition boundaryCondition) {
    large_num nodes_per_layer = nodes->size() / numberOfLayers;

    large_num nodeID;
    large_num nodeID_neig;
    std::vector<large_num> nodeIDs_neig;
    large_num spatID;
    large_num spatID_neig;
    for (int layer = 0; layer < numberOfLayers; ++layer) {
        for (int i = 0; i < nodes_per_layer; ++i) {
            nodeID = i + (nodes_per_layer * layer);
            addNeighbourOrBoundary(nodeID, boundaryConduct, boundaryCondition, layer, resolution);
        }
    }

}

void addNeighbourOrBoundary(NodeVector nodes, large_num nodeID, double boundaryConduct, Simulation::Options::BoundaryCondition boundaryCondition, int layer, int resolution){

    std::unordered_map<int, Model::NeighbourPosition> lu;
    // unrefined positions
    lu[0] = Model::FRONT; // formerly NORTH
    lu[1] = Model::BACK; // formerly SOUTH
    lu[2] = Model::RIGHT; // formerly EAST
    lu[3] = Model::LEFT; // formerly WEST
    // refined positions
    lu[4] = Model::FRONTLEFT;
    lu[5] = Model::FRONTRIGHT;
    lu[6] = Model::BACKLEFT;
    lu[7] = Model::BACKRIGHT;
    lu[8] = Model::RIGHTFRONT;
    lu[9] = Model::RIGHTBACK;
    lu[10] = Model::LEFTFRONT;
    lu[11] = Model::LEFTBACK;

    int rowLength{360 * 60 * 60 / resolution}; // todo change so that simple models work too

    for (int j = 0; j < 4; ++j) {
        spatID = nodes->at(nodeID)->getSpatID();
        spatID_neig = getNeighbourLatLon(spatID, j, resolution);
        if (spatID_neig > 0) {
            if (spatIDtoNodeIDs.contains(spatID_neig) and
                not spatIDtoNodeIDs.at(spatID_neig).empty()) { //Neighbour id exists in landmask
                nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig);
                if (nodeIDs_neig.empty()){
                    continue;
                }
                nodeID_neig = nodeIDs_neig[layer];
                nodes->at(nodeID)->setNeighbour(nodeID_neig, lu[j]);
            } else {// Neighbour is not in landmask
                // Calling function defined above to add boundary
                //LOG(userinfo) << "adding boundary at nodeID" << nodeID;
                addBoundary(nodes, boundaryConduct, boundaryCondition, nodeID, layer, lu[j]);
            }
        }
    }
}

/***
 *
 */
void // defining function to add boundaries to nodes with no neighbours (called later in code)
addBoundary(NodeVector nodes, double boundaryConduct, Simulation::Options::BoundaryCondition boundaryCondition,
            large_num pos, int layer, Model::NeighbourPosition positionOfBoundary) {

    if (layer > 0) {
        return;
    }

    switch (boundaryCondition) {
        case Simulation::Options::GENERAL_HEAD_NEIGHBOUR: {
            auto head = nodes->at(pos)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
            nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, head, boundaryConduct, head);
        }
            break;
        case Simulation::Options::GENERAL_HEAD_BOUNDARY: {
            nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, 0 * Model::si::meter,
                                            boundaryConduct, 0 * Model::si::meter);
        }
        default:
            break;
    }
};

/**
 * @brief calculates neighbour based on lat and lon
 * Assumes no landmask and calculates all possible neighbours
 * @param id
 * @param res
 * @return
 */
    large_num getNeighbourLatLon(double lat, double lon, int j, int rowLength) {

        assert(row_l % 2 == 0 && "resolution is impossible");
        switch(j) {
            case 0:
                if (spatID > row_l) {
                    //NORTH
                    return spatID - row_l;
                }
            case 1:
                if (spatID < (row_l / 2) * row_l - row_l) {
                    //SOUTH
                    return spatID + row_l;
                }
            case 2:
                if (spatID % row_l == 0) {
                    //EAST
                    return spatID - row_l + 1;
                } else { return spatID + 1; }
            case 3:
                if ((spatID - 1) % row_l == 0) {
                    //WEST
                    return spatID + row_l - 1;
                } else { return spatID - 1; }
            default:
                return -1;
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
void buildBySpatID(NodeVector nodes, const std::unordered_map<large_num, std::vector<large_num>> spatIDtoNodeIDs,
                   int resolution, int numberOfLayers, double boundaryConduct,
                   Simulation::Options::BoundaryCondition boundaryCondition) {
        large_num nodes_per_layer = nodes->size() / numberOfLayers;

    // defining function to add boundaries to nodes with no neighbours (called later in code)
    auto addBoundary = [nodes, boundaryConduct, boundaryCondition](
            large_num pos, int layer, Model::NeighbourPosition positionOfBoundary) {

        if (layer > 0) {
            return;
        }

        switch (boundaryCondition) {
            case Simulation::Options::GENERAL_HEAD_NEIGHBOUR: {
                auto head = nodes->at(pos)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, head, boundaryConduct, head);
            }
                break;
            case Simulation::Options::GENERAL_HEAD_BOUNDARY: {
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY, 0 * Model::si::meter,
                                                boundaryConduct, 0 * Model::si::meter);
            }
            default:
                break;
        }
    };

    std::unordered_map<int, Model::NeighbourPosition> lu;
    lu[0] = Model::FRONT; // formerly NORTH
    lu[1] = Model::BACK; // formerly SOUTH
    lu[2] = Model::RIGHT; // formerly EAST
    lu[3] = Model::LEFT; // formerly WEST

    large_num nodeID;
    large_num nodeID_neig;
    std::vector<large_num> nodeIDs_neig;
    large_num spatID;
    large_num spatID_neig;
    for (int layer = 0; layer < numberOfLayers; ++layer) {
        for (int i = 0; i < nodes_per_layer; ++i) {
            nodeID = i + (nodes_per_layer * layer);
            for (int j = 0; j < 4; ++j) {
                spatID = nodes->at(nodeID)->getSpatID();
                spatID_neig = getNeighbourSpatID(spatID, j, resolution);
                if (spatID_neig > 0) {
                    if (spatIDtoNodeIDs.contains(spatID_neig) and
                        not spatIDtoNodeIDs.at(spatID_neig).empty()) { //Neighbour id exists in landmask
                        nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig);
                        if (nodeIDs_neig.empty()){
                            continue;
                        }
                        nodeID_neig = nodeIDs_neig[layer];
                        nodes->at(nodeID)->setNeighbour(nodeID_neig, lu[j]);
                    } else {// Neighbour is not in landmask
                        // Calling function defined above to add boundary
                        //LOG(userinfo) << "adding boundary at nodeID" << nodeID;
                        addBoundary(nodeID, layer, lu[j]);
                    }
                }
            }
        }
    }

}

/**
 * @brief calculates all neighbours based on a spatial ID
 * Assumes no landmask and calculates all possible neighbours
 * @param id
 * @param res
 * @return
 */
large_num getNeighbourSpatID(large_num spatID, int j, int res) {
    int row_l{360 * 60 * 60 / res};
    assert(row_l % 2 == 0 && "resolution is impossible");
    switch(j) {
        case 0:
            if (spatID > row_l) {
                //NORTH
                return spatID - row_l;
            }
        case 1:
            if (spatID < (row_l / 2) * row_l - row_l) {
                //SOUTH
                return spatID + row_l;
            }
        case 2:
            if (spatID % row_l == 0) {
                //EAST
                return spatID - row_l + 1;
            } else { return spatID + 1; }
        case 3:
            if ((spatID - 1) % row_l == 0) {
                //WEST
                return spatID + row_l - 1;
            } else { return spatID - 1; }
        default:
            return -1;
    }
}


/**
 * @brief
 * @param from is position in vector of top layer node
 * @param to is position in vector of node that receive neighbouring information
 */
void copyNeighbour(size_t from, size_t to, NodeVector nodes, size_t layer_shift){
    auto neighbours = nodes->at(from)->getListOfNeighbours();
    for(const auto &n : neighbours){
        if(n.first == Model::DOWN or n.first == Model::TOP){
            continue;
        }
        nodes->at(to)->setNeighbour(n.second + layer_shift ,n.first);
    }
}

/**
 * @brief Copies cardinal points of top layer to all bottom layers
 * @param nodes
 * @param layers
 */
void copyNeighboursToBottomLayers(NodeVector nodes, int numberOfLayers){
    assert(numberOfLayers && "0 layers does not make sense");
    if (numberOfLayers == 1) {
        return;
    }
    size_t nodesPerLayer = nodes->size() / numberOfLayers;
    for (int i = 0; i < nodesPerLayer; ++i) {
        for (int j = 0; j < numberOfLayers - 1; ++j) {
            copyNeighbour(i ,i + (nodesPerLayer * j), nodes, (int) nodesPerLayer  * j);
        }
    }
}

/**
 * Adds additional layers and connects nodes
 * @param nodes
 * @param layers
 * @param conf
 * @param aquifer_thickness
 */
void buildBottomLayers(NodeVector nodes,
                       int numberOfLayers,
                       std::vector<bool> conf,
                       std::vector<int> aquifer_thickness,
                       std::vector<double> conductances,
                       std::vector<double> anisotropies) {
    assert(numberOfLayers && "AsModel::signing 0 layers does not make any sense");
    if (numberOfLayers == 1) {
        return;
    }

    size_t nodesPerLayer = nodes->size();
    nodes->reserve(numberOfLayers * nodesPerLayer);

    LOG(debug) << "Building additional layers with node count: " << nodesPerLayer << " for " << numberOfLayers << " layers";

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
    int refinementLevel;
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
            edgeLengthLeftRight = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::EdgeLengthLeftRight>();
            edgeLengthFrontBack = nodes->at(i)->getProperties().get<Model::quantity<Model::Meter>, Model::EdgeLengthFrontBack>();
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
            refinementLevel = nodes->at(i)->getProperties().get<int, Model::RefinementLevel>();
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
                                                            conf[j + 1],
                                                            refinementLevel,
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
    LOG(debug) << "Last nodeID was " << id << " with max ID (with non static nodes) " << nodesPerLayer * numberOfLayers;
};

}
}//ns
