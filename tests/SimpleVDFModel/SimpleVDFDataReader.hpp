#ifndef TESTING_SIMPLEVDFDATAREADER_HPP
#define TESTING_SIMPLEVDFDATAREADER_HPP

#include "../../src/DataProcessing/DataReader.hpp"

namespace GlobalFlow {
namespace DataProcessing {


class SimpleVDFDataReader : public DataReader {
    public:
        SimpleVDFDataReader(int step) { stepMod = step; }

        virtual void readData(Simulation::Options op) {
            LOG(userinfo) << "Building the initial model layer";
            std::vector<std::vector<int>> grid;
            grid = readGrid(nodes,
                            buildDir(op.getNodesDir()),
                            op.getNumberOfNodes(),
                            op.getNumberOfRows(),
                            op.getNumberOfCols(),
                            op.getInitialK(),
                            op.getAquiferDepth()[0],
                            op.getAnisotropy(),
                            op.getSpecificYield(),
                            op.getSpecificStorage(),
                            op.getEdgeLengthLeftRight(),
                            op.getEdgeLengthFrontBack(),
                            op.isConfined(0));

            LOG(userinfo) << "Reading hydraulic parameters";
            readConduct(buildDir(op.getLithology()));
            readElevation(buildDir(op.getElevation()));

            LOG(userinfo) << "Reading the groundwater recharge";
            readGWRecharge(buildDir(op.getRecharge()));

            LOG(userinfo) << "Reading the boundary condition";
            readHeadBoundary(buildDir(op.getKGHBDir()));

            LOG(userinfo) << "Reading variable density information";
            readVariableDensity(op.isDensityVariable(), op.isDensityStratified(), op.getMaxDensity(), op.getMinDensity(),
                                op.getNumberOfDensityZones());

            LOG(userinfo) << "Initializing head";
            readInitialHeads((buildDir(op.getInitialHeadsDir())));

            LOG(userinfo) << "Connecting the model cells";
            DataProcessing::buildByGrid(nodes, grid, op.getNumberOfLayers(), op.getGHBConduct(),
                                        op.getBoundaryCondition());
        }

    private:
        template<class T>
        using Matrix = std::vector<std::vector<T>>;

        Matrix<int>
        readGrid(NodeVector nodes, std::string path, int numberOfNodes, int numberOfRows, int numberOfCols,
                 double defaultK,
                 double aquiferDepth,
                 double anisotropy,
                 double specificYield,
                 double specificStorage,
                 double edgeLengthLeftRight,
                 double edgeLengthFrontBack,
                 bool confined) {
            Matrix<int> out = Matrix<int>(numberOfCols, std::vector<int>(numberOfRows));

            io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "global_ID", "X", "Y", "cell_area", "row", "col");

            double x{0};
            double y{0};
            double area{0};
            int globid{0};
            int i{0};
            int row{0};
            int col{0};
            lookupglobIDtoID.reserve(numberOfNodes);

            while (in.read_row(globid, x, y, area, row, col)) {
                out[row][col] = i;
                nodes->emplace_back(new Model::StandardNode(nodes,
                                                            x,
                                                            y,
                                                            area * Model::si::square_meter,
                                                            edgeLengthLeftRight * Model::si::meter,
                                                            edgeLengthFrontBack * Model::si::meter,
                                                            (unsigned long) globid,
                                                            i,
                                                            defaultK * (Model::si::meter / Model::day),
                                                            stepMod,
                                                            aquiferDepth,
                                                            anisotropy,
                                                            specificYield,
                                                            specificStorage, confined));
                lookupglobIDtoID[globid] = i;
                i++;
            }

            return out;
        }

        void readConduct(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                if (data > 10) {
                    LOG(debug) << "Very high conductance value at global_ID "<<pos<<". Possible Data Error";
                }
                nodes->at(pos)->setK(data * (Model::si::meter / Model::day));
            });
        };

        void
        readElevation(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->setElevation(data * Model::si::meter);
            });
        };

        void readHeadBoundary(std::string path) {
            io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "global_ID", "elevation", "conduct");
            int arcid{0};
            double elevation{0};
            double conduct{0};

            while (in.read_row(arcid, elevation, conduct)) {
                int pos = 0;
                try {
                    pos = lookupglobIDtoID.at(arcid);
                }
                catch (const std::out_of_range &ex) {
                    //if Node does not exist ignore entry
                    continue;
                }
                nodes->at(pos)->addExternalFlow(Model::GENERAL_HEAD_BOUNDARY,
                                                elevation * Model::si::meter,
                                                conduct,
                                                elevation * Model::si::meter);
            }
        }

        void readGWRecharge(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->addExternalFlow(Model::RECHARGE,
                                                0 * Model::si::meter,
                                                data * nodes->at(pos)->getProperties().get < Model::quantity <
                                                Model::SquareMeter > ,
                                                Model::Area > ().value(),
                                                0 * Model::si::meter);
            });
        }

        void readVariableDensity(bool densityVariable, bool densityStratified, vector<double> densities,
                                 int densityZones) {

            if (densityVariable){
                vector<double> delnusZone;
                vector<double> epsZone;
                vector<double> nusZeta;
                vector<double> nusZone;
                // TODO make variables Model::si::si_dimensionless
                double densityFresh = 1000;

                for (int n = 1; n <= densityZones+1; n++) {
                    nusZeta[n] = ( densities[n] - densityFresh ) / densityFresh;
                }

                for (int n = 1; n <= densityZones; n++) {
                    if (densityStratified) {
                        nusZone[n] = nusZeta[n+1];
                        epsZone[n] = 0;
                    } else { // if continuous
                        nusZone[n] = 0.5 * ( nusZeta[n] + nusZeta[n+1] );
                        epsZone[n] = ( nusZeta[n+1] - nusZeta[n] ) / 6;
                    }

                    if (n == 1) {
                        delnusZone[n] = nusZone[n];
                    } else {
                        delnusZone[n] = nusZone[n] - nusZone[n-1];
                    }
                }
            }
        }

        void readInitialHeads(std::string path) {
            readTwoColumns(path, [this](double data, int pos) {
                nodes->at(pos)->setHead_direct(data);
            });
        }



};
}
}
#endif //TESTING_SIMPLEVDFDATAREADER_HPP

