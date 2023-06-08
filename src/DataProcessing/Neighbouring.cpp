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
void buildBySpatID(NodeVector nodes, const std::unordered_map<large_num, std::vector<large_num>> spatIDtoNodeIDs,
                   large_num resolution, int numberOfLayers, double boundaryConduct,
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
                Model::quantity<Model::Meter> head =
                        nodes->at(pos)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
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
    int spatID_neig;
    for (int layer = 0; layer < numberOfLayers; ++layer) {
        for (int i = 0; i < nodes_per_layer; ++i) {
            nodeID = i + (nodes_per_layer * layer);
            for (int j = 0; j < 4; ++j) {
                spatID = nodes->at(nodeID)->getSpatID();
                spatID_neig = getNeighbourSpatID(spatID, j, resolution);
                if (spatID_neig > 0) { //
                    if (spatIDtoNodeIDs.contains((unsigned long) spatID_neig) and
                        not spatIDtoNodeIDs.at((unsigned long) spatID_neig).empty()) { //Neighbour id exists in landmask
                        nodeIDs_neig = spatIDtoNodeIDs.at(spatID_neig);
                        nodeID_neig = spatIDtoNodeIDs.at(spatID_neig)[layer];
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
/*
n_array getNeighboursBySpatID(large_num spatID, int res) {
    n_array neighbours{-1, -1, -1, -1};
    large_num row_l{360 * 60 * 60 / res};
    assert(row_l % 2 == 0 && "resolution is impossible");

    if (spatID > row_l) {
        //NORTH
        neighbours[0] = spatID - row_l;
    }
    if (spatID < (row_l / 2) * row_l - row_l) {
        //SOUTH
        neighbours[1] = spatID + row_l;
    }
    if (spatID % row_l == 0) {
        //EAST
        neighbours[2] = spatID - row_l + 1;
    } else { neighbours[2] = spatID + 1; }
    if ((spatID - 1) % row_l == 0) {
        //WEST
        neighbours[3] = spatID + row_l - 1;
    } else { neighbours[3] = spatID - 1; }
    return neighbours;
}*/

/**
 * @brief calculates all neighbours based on a spatial ID
 * Assumes no landmask and calculates all possible neighbours
 * @param id
 * @param res
 * @return
 */
int getNeighbourSpatID(int spatID, int j, int res) {
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
 * @brief Connects neighbouring nodes
 * @param nodes
 * @param numberOfTOPNodes
 * @param layers
 * @param ghbConduct
 * @param boundaryCondition
 * @return Number of new nodes
 */
int buildNeighbourMap(NodeVector nodes, int numberOfTOPNodes, int layers, double ghbConduct,
                      Simulation::Options::BoundaryCondition boundaryCondition) {
    //Key is x-poModel::sition of node, value node ID
    std::unordered_map<double, int> previousRow;
    std::unordered_map<double, int> currentRow;
    previousRow.reserve(1000);
    currentRow.reserve(1000);
    previousRow[0] = 0;

    auto id = nodes->size() - 1;
    int numOfStaticHeads{0};

    //auto rd = [](double num){return (int) num * 10000;};

    auto setNeighbouring = [nodes](large_num positionOfNewNode,
                                   large_num positionOfExistingNode,
                                   Model::NeighbourPosition pos1,
                                   Model::NeighbourPosition pos2) {
        nodes->at(positionOfNewNode)->setNeighbour(positionOfExistingNode, pos1);
        nodes->at(positionOfExistingNode)->setNeighbour(positionOfNewNode, pos2);
    };

    auto getLat = [nodes](large_num pos) {
        return nodes->at(pos)->getProperties().get<double, Model::Lat>();
    };

    auto getLon = [nodes](large_num pos) {
        return nodes->at(pos)->getProperties().get<double, Model::Lon>();
    };

    auto addBoundary = [nodes, ghbConduct, boundaryCondition, id, &numOfStaticHeads, setNeighbouring](
            large_num pos, int layer, Model::NeighbourPosition positionOfBoundary) {
        if (layer > 0) {
            return;
        }

        switch (boundaryCondition) {
            case Simulation::Options::GENERAL_HEAD_NEIGHBOUR: {
                Model::quantity<Model::Meter> head =
                        nodes->at(pos)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                head,
                                                ghbConduct,
                                                head);
            }
                break;
            case Simulation::Options::GENERAL_HEAD_BOUNDARY: {
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                0 * Model::si::meter,
                                                ghbConduct,
                                                0 * Model::si::meter);
            }
                break;
            case Simulation::Options::STATIC_HEAD_SEA_LEVEL: {
                LOG(debug) << "Model::Using static head boundary";
                //Add a constant head boundary
                auto staticID = id + numOfStaticHeads;
                Model::quantity<Model::SquareMeter> area = 1 * Model::si::square_meter;
                Model::quantity<Model::Meter> edgeLengthLeftRight = 1 * Model::si::meter;
                Model::quantity<Model::Meter> edgeLengthFrontBack = 1 * Model::si::meter;
                nodes->emplace_back(new Model::StaticHeadNode(nodes, staticID, area, edgeLengthLeftRight,
                                                              edgeLengthFrontBack));

                switch (positionOfBoundary) {
                    case Model::LEFT:
                        setNeighbouring(pos, staticID, Model::LEFT, Model::RIGHT);
                        break;
                    case Model::RIGHT:
                        setNeighbouring(pos, staticID, Model::RIGHT, Model::LEFT);
                        break;
                    case Model::BACK:
                        setNeighbouring(pos, staticID, Model::BACK, Model::FRONT);
                        break;
                    case Model::FRONT:
                        setNeighbouring(pos, staticID, Model::FRONT, Model::BACK);
                        break;
                    default :
                        break;
                }
                numOfStaticHeads++;
            }
                break;
            default:
                break;
        }
    };

    const double epsLon{0.01}; //allow a minimal deviation
    const double epsLat{0.125}; //Allow 5' + 1/2 5', 5 arcmin = 0.08333 in decimal degree

    //FIXME Inefficient
    //currently does the same work for all layers all over again
    for (int j = 0; j < layers; ++j) {
        previousRow.clear();
        currentRow.clear();

        //First node
        //-> add left GHB
        addBoundary(0 + (j * numberOfTOPNodes), j, Model::LEFT);

        currentRow[getLat(0 + (j * numberOfTOPNodes))] = 0 + (j * numberOfTOPNodes);
        //FrontN = TopRow
        //BackN = BottomRow
        for (int i = 1 + (j * numberOfTOPNodes); i < numberOfTOPNodes + (j * numberOfTOPNodes); ++i) {

            //Still in same row?
            //
            if (std::abs(getLon(i) - getLon(i - 1)) <= epsLon) {
                if (std::abs(getLat(i) - getLat(i - 1)) < epsLat) {
                    //Previous node is direct neighbour in x-direction
                    setNeighbouring(i, i - 1, Model::LEFT, Model::RIGHT);
                } else {
                    //Still in same row
                    //But x diff > Model::epsilon thus add GHB cell
                    addBoundary(i, j, Model::RIGHT);
                    addBoundary(i + 1, j, Model::LEFT);
                }
            } else {
                //New row

                //AsModel::sign a GHB to last in row to the right
                addBoundary(i - 1, j, Model::RIGHT);

                //If there are nodes left with no back node asModel::signed -> asModel::sign an GHB
                for (const auto &node : previousRow) {
                    //Nodes which were not asModel::signed asModel::sign GHB
                    addBoundary(node.second, j, Model::BACK);
                }

                previousRow = currentRow;
                currentRow.clear();

                //First left is always GHB
                addBoundary(i, j, Model::LEFT);
            }

            if (previousRow.empty()) {
                // Should only be at first row
                //. add top GHB
                addBoundary(i, j, Model::FRONT);
            } else {
                //Not first Row
                //Set front and back
                auto val = getLat(i);
                std::unordered_map<double, int>::const_iterator item = std::find_if(previousRow.begin(),
                                                                                    previousRow.end(), [val]
                                                                                            (std::pair<const double, int> item) {
                            return std::abs(val - item.first) < 0.001;
                        });
                if (item != previousRow.end()) {
                    //Node at same x existed
                    nodes->at(i)->setNeighbour(nodes->at(item->second)->getProperties().get<large_num, Model::ID>(),
                                               Model::FRONT);
                    nodes->at(item->second)->setNeighbour(nodes->at(i)->getProperties().get<large_num, Model::ID>(),
                                                          Model::BACK);
                    //Delete already used member
                    previousRow.erase(item->first);
                } else {
                    //AsModel::sign GHB to others
                    addBoundary(i, j, Model::FRONT);
                }
            }

            //Store current row
            currentRow[getLat(i)] = i;
        }
        //That was the last row add GHB nodes to all at BACK
        for (const auto &item : currentRow) {
            addBoundary(item.second, j, Model::BACK);
        }

        //AsModel::sign a GHB to last in row to the right
        addBoundary(id + j * id, j, Model::RIGHT);
    }
    return numOfStaticHeads;
};

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
    int stepMod;
    Model::quantity<Model::SquareMeter> area;
    Model::quantity<Model::Meter> edgeLengthLeftRight;
    Model::quantity<Model::Meter> edgeLengthFrontBack;
    Model::quantity<Model::Velocity> K;
    Model::quantity<Model::Meter> head;
    double aquiferDepth;
    double anisotropy;
    double specificYield;
    double specificStorage;
    bool densityVariable;
    std::vector<Model::quantity<Model::Meter>> zetas;
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
            stepMod = nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                    Model::StepModifier>();
            aquiferDepth = aquifer_thickness[j + 1];
            anisotropy = anisotropies[j + 1];
            specificYield =
                    nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                            Model::SpecificYield>().value();
            specificStorage =
                    nodes->at(i)->getProperties().get<Model::quantity<Model::perUnit>, Model::SpecificStorage>
                            ().value();
            densityVariable = nodes->at(i)->getProperties().get<bool, Model::DensityVariable>();
            delnus = nodes->at(i)->getProperties().
                    get<std::vector<Model::quantity<Model::Dimensionless>>, Model::Delnus>();
            nusInZones = nodes->at(i)->getProperties().
                    get<Model::vector<Model::quantity<Model::Dimensionless>>, Model::NusInZones>();
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
