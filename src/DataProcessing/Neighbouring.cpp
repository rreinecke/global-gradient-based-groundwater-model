#include "Neighbouring.hpp"

namespace GlobalFlow {
    namespace DataProcessing {

/**
 * For regular grids (modflow like)
 * @param nodes
 * @param grid
 * @param layers
 * @param oceanCoduct
 * @param staticHeadBoundary
 */
        void buildByGrid(NodeVector nodes, Matrix<int> grid, int layers, double oceanCoduct, bool staticHeadBoundary) {
            //id->row,col
            int rows = grid[0].size();
            int cols = grid.size();

            auto check = [grid](int i, int j) {
                try {
                    int p = grid.at(i).at(j);
                    return true;
                } catch (std::out_of_range e) {
                    return false;
                }
            };

            for (int layer = 0; layer < layers; ++layer) {

                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        int id = grid[i][j];
                        int l_mult = layer * (rows * cols);
                        if (check(i + 1, j)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i + 1][j] + l_mult, Model::EAST);
                        }
                        if (check(i - 1, j)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i - 1][j] + l_mult, Model::WEST);
                        }
                        if (check(i, j - 1)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i][j - 1] + l_mult, Model::NORTH);
                        }
                        if (check(i, j + 1)) {
                            nodes->at(id + l_mult)->setNeighbour(grid[i][j + 1] + l_mult, Model::SOUTH);
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
 * @param oceanCoduct
 * @param boundaryCondition
 */
        void buildBySpatID(NodeVector nodes, const std::unordered_map<int, int> id_mapping, int resolution, int layers,
                           double oceanCoduct, Simulation::Options::BoundaryCondition boundaryCondition) {
            int nodes_per_layer = nodes->size() / layers;

            auto addBoundary = [nodes, oceanCoduct, boundaryCondition](
                           large_num pos, int layer,
                           Model::NeighbourPosition
                           positionOfBoundary) {
                if (layer > 0) {
                    return;
                }

                switch (boundaryCondition) {
                    case Simulation::Options::CONSTANT_HEAD_NEIGHBOUR: {
                        Model::quantity<Model::Meter> head =
                                       nodes->at(pos)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
                        nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                        head,
                                                        oceanCoduct,
                                                        head);
                    }
                        break;
                    case Simulation::Options::CONSTANT_HEAD_SEA_LEVEL: {
                        nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                        0 * Model::si::meter,
                                                        oceanCoduct,
                                                        0 * Model::si::meter);
                    }
                    default:
                        break;
                }
            };


            std::unordered_map<int, Model::NeighbourPosition> lu;
            lu[0] = Model::NORTH;
            lu[1] = Model::SOUTH;
            lu[2] = Model::EAST;
            lu[3] = Model::WEST;

            for (int i = 0; i < nodes_per_layer; ++i) {
                n_array nei = getNeighbourBySpatialID((int) nodes->at(i)->getSpatialID(), resolution);
                for (int j = 0; j < 4; ++j) {
                    if (id_mapping.count(nei[j]) > 0) {
                        //Neighbour id exists in landmask
                        nodes->at(i)->setNeighbour(id_mapping.at(nei[j]), lu[j]);
                    } else {
                        addBoundary(i, 0, lu[j]);
                    }
                }
            }

        }

/**
 * @brief calculates all neighbours based on a spatial ID
 * Assumes no landmask and caclulates all possible neighbours
 * @param id
 * @param res
 * @return
 */
        n_array getNeighbourBySpatialID(int id, int res) {
            n_array neighbours{-1, -1, -1, -1};
            int row_l{360 * 60 * 60 / res};
            assert(row_l % 2 == 0 && "resolution is impossible");

            if (id > row_l) {
                //NORTH
                neighbours[0] = id - row_l;
            }
            if (id < (row_l / 2) * row_l - row_l) {
                //SOUTH
                neighbours[1] = id + row_l;
            }
            if (id % row_l == 0) {
                //EAST
                neighbours[2] = id - row_l + 1;
            } else { neighbours[2] = id + 1; }
            if ((id - 1) % row_l == 0) {
                //WEST
                neighbours[3] = id + row_l - 1;
            } else { neighbours[3] = id - 1; }
            return neighbours;
        }


/**
 * @brief Connects neighbouring nodes
 * @param nodes
 * @param numberOfTOPNodes
 * @param layers
 * @param oceanCoduct
 * @param boundaryCondition
 * @return Number of new nodes
 */
        int buildNeighbourMap(NodeVector nodes, int numberOfTOPNodes, int layers, double oceanCoduct,
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

            auto addBoundary = [nodes, oceanCoduct, boundaryCondition, id, &numOfStaticHeads, setNeighbouring](
                           large_num pos, int layer,
                           Model::NeighbourPosition
                           positionOfBoundary) {
                if (layer > 0) {
                    return;
                }

        switch (boundaryCondition) {
            case Simulation::Options::CONSTANT_HEAD_NEIGHBOUR: {
                Model::quantity<Model::Meter> head =
                        nodes->at(pos)->getProperties().get<Model::quantity<Model::Meter>, Model::EQHead>();
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                head,
                                                oceanCoduct,
                                                head);
            }
                break;
            case Simulation::Options::CONSTANT_HEAD_SEA_LEVEL: {
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                0 * Model::si::meter,
                                                oceanCoduct,
                                                0 * Model::si::meter);
            }
                break;
            case Simulation::Options::STATIC_HEAD_SEA_LEVEL: {
                        LOG(debug) << "UModel::sing static head boundary";
                        //Add a constant head boundary
                        auto staticID = id + numOfStaticHeads;
                        Model::quantity<Model::SquareMeter> area = 1 * Model::si::square_meter;
                        nodes->emplace_back(new Model::StaticHeadNode(nodes, staticID, area));

                        switch (positionOfBoundary) {
                            case Model::WEST:
                                setNeighbouring(pos, staticID, Model::WEST, Model::EAST);
                                break;
                            case Model::EAST:
                                setNeighbouring(pos, staticID, Model::EAST, Model::WEST);
                                break;
                            case Model::SOUTH:
                                setNeighbouring(pos, staticID, Model::SOUTH, Model::NORTH);
                                break;
                            case Model::NORTH:
                                setNeighbouring(pos, staticID, Model::NORTH, Model::SOUTH);
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

            //FIXME Inifficient
            //currently does the same work for all layers all over again
            for (int j = 0; j < layers; ++j) {
                previousRow.clear();
                currentRow.clear();

                //First node
                //-> add left Ocean
                addBoundary(0 + (j * numberOfTOPNodes), j, Model::WEST);

                currentRow[getLat(0 + (j * numberOfTOPNodes))] = 0 + (j * numberOfTOPNodes);
                //FrontN = TopRow
                //BackN = BottomRow
                for (int i = 1 + (j * numberOfTOPNodes); i < numberOfTOPNodes + (j * numberOfTOPNodes); ++i) {

                    //Still in same row?
                    //
                    if (std::abs(getLon(i) - getLon(i - 1)) <= epsLon) {
                        if (std::abs(getLat(i) - getLat(i - 1)) < epsLat) {
                            //Previous node is direct neighbour in x-direction
                            setNeighbouring(i, i - 1, Model::WEST, Model::EAST);
                        } else {
                            //Still in same row
                            //But x diff > epModel::silon thus add ocean cell
                            addBoundary(i, j, Model::EAST);
                            addBoundary(i + 1, j, Model::WEST);
                        }
                    } else {
                        //New row

                        //AsModel::sign an ocean to last in row to the right
                        addBoundary(i - 1, j, Model::EAST);

                        //If there are nodes left with no back node asModel::signed -> asModel::sign an ocean
                        for (const auto &node : previousRow) {
                            //Nodes which were not asModel::signed asModel::sign ocean
                            addBoundary(node.second, j, Model::SOUTH);
                        }

                        previousRow = currentRow;
                        currentRow.clear();

                        //First left is always ocean
                        addBoundary(i, j, Model::WEST);
                    }

                    if (previousRow.empty()) {
                        // Should only be at first row
                        //. add top Ocean
                        addBoundary(i, j, Model::NORTH);
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
                            nodes->at(i)->setNeighbour(
                                           nodes->at(item->second)->getProperties().get<large_num, Model::ID>(),
                                           Model::NORTH);
                            nodes->at(item->second)->setNeighbour(
                                           nodes->at(i)->getProperties().get<large_num, Model::ID>(),
                                           Model::SOUTH);
                            //Delete already used member
                            previousRow.erase(item->first);
                        } else {
                            //AsModel::sign ocean to others
                            addBoundary(i, j, Model::NORTH);
                        }
                    }

                    //Store current row
                    currentRow[getLat(i)] = i;
                }
                //That was the last row add ocean nodes to all at SOUTH
                for (const auto &item : currentRow) {
                    addBoundary(item.second, j, Model::SOUTH);
                }

                //AsModel::sign a ocean to last in row to the right
                addBoundary(id + j * id, j, Model::EAST);
            }
            return numOfStaticHeads;
        };

        /**
         * @brief
         * @param from is position in vector of top layer node
         * @param to is position in vector of node that recieve neighbouring information
         */
        void copyNeighbour(size_t from, size_t to, NodeVector nodes, int layer_shift){
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
        void copyNeighboursToBottomLayers(NodeVector nodes, int layers){
            assert(layers && "AsModel::signing 0 layers does not make any sense");
            if (layers == 1) {
                return;
            }
            size_t top_layer_size = nodes->size() / layers;
            for (int i = 0; i < top_layer_size; ++i) {
                for (int j = 0; j < layers - 1; ++j) {
                    copyNeighbour(i ,i + (top_layer_size * j), nodes, top_layer_size  * j);
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
        void
        buildBottomLayers(NodeVector nodes, int layers, std::vector<bool> conf, std::vector<int> aquifer_thickness) {
            assert(layers && "AsModel::signing 0 layers does not make any sense");
            if (layers == 1) {
                return;
            }

            size_t layersize = nodes->size();
            nodes->reserve(layers * layersize);

            LOG(debug) << "Building additional layers with node count: " << layersize << " for " << layers << " layers";

            size_t id = layersize;
            large_num arcID;
            double lat, lon;
            int stepMod;
            Model::quantity<Model::SquareMeter> area;
            Model::quantity<Model::Velocity> K;
            double aquiferDepth;
            double anisotropy;
            double specificYield;
            double specificStorage;

            for (int j = 0; j < layers - 1; ++j) {
                //1) Add a Model::similar node in z direction for each layer
                //TODO Parallell?
                for (int i = 0; i < layersize; ++i) {
                    //for each node in top layer

                    arcID = nodes->at(i)->getProperties().get<large_num, Model::SpatID>();
                    lat = nodes->at(i)->getProperties().get<double, Model::Lat>();
                    lon = nodes->at(i)->getProperties().get<double, Model::Lon>();
                    area = nodes->at(i)->getProperties().get<Model::quantity<Model::SquareMeter>, Model::Area>();
                    K = nodes->at(i)->getK__pure();
                    stepMod = nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                                   Model::StepModifier>();
                    aquiferDepth = aquifer_thickness[j + 1];
                    anisotropy = nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                                   Model::Anisotropy>().value();
                    specificYield =
                                   nodes->at(i)->getProperties().get<Model::quantity<Model::Dimensionless>,
                                                  Model::SpecificYield>().value();
                    specificStorage =
                                   nodes->at(i)->getProperties().get<Model::quantity<Model::perUnit>, Model::SpecificStorage>
                                                  ().value();

                    if (nodes->at(i)->isStaticNode()) {
                        //is taken care of by neighbouring algorithm
                        continue;
                    } else {
                        if (id > layersize * layers) {
                            LOG(critical) << "This is not possible!";
                            exit(9);
                        }
                        nodes->emplace_back(new Model::StandardNode(nodes, lat, lon, area, arcID, id, K, stepMod,
                                                                    aquiferDepth,
                                                                    anisotropy,
                                                                    specificYield,
                                                                    specificStorage, conf[j + 1]));
                        nodes->at(id)->getProperties().set<int, Model::Layer>(j + 1);
                        nodes->at(id)->getProperties().set<Model::quantity<Model::Meter>, Model::Elevation>(
                                       nodes->at(id)->getProperties().get<Model::quantity<Model::Meter>, Model::Elevation>()
                                       -
                                       (aquiferDepth *
                                        Model::si::meter));
                    }
                    //2) Neighbouring for top and bottom

                    if (j > 0) {
                        //Layer above is not top layer
                        nodes->at(id)->setNeighbour(i + (j * layersize), Model::TOP);
                        nodes->at(i + (j * layersize))->setNeighbour(id, Model::DOWN);
                    } else {
                        //Layer above is top layer
                        nodes->at(id)->setNeighbour(i, Model::TOP);
                        nodes->at(i)->setNeighbour(id, Model::DOWN);
                    }
                    id++;
                    if (id > (layersize * layers) - 1) {
                        break;
                    }
                }
            }
            LOG(debug) << "Last nodeID was " << id << " with max ID (with non static nodes) " << layersize * layers;
        };

    }
}//ns
