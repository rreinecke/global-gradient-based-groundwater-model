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

#ifndef GLOBAL_FLOW_OUTPUTFACT_HPP
#define GLOBAL_FLOW_OUTPUTFACT_HPP

#include "Converter.hpp"
#include "../../Logging/Logging.hpp"

#include <ctime>
#include <iomanip>
//#include <netcdf>

namespace GlobalFlow {
    namespace DataProcessing {
        namespace DataOutput {

            /**
             * @enum Output Type
             * @brief What kind of output format
             */
            enum class OutputType {
                CSV,
                NET_CDF,
                GFS_JSON,
                NON_VALID
            };

            /**
             * @var outputMapping
             * Mapping of string representation to internal field
             */
            const std::unordered_map<std::string, OutputType> outputMapping{
                    {"csv",      OutputType::CSV},
                    {"netcdf",   OutputType::NET_CDF},
                    {"gfs-json", OutputType::GFS_JSON}
            };

            using path = std::string;
            using pos_v = std::vector<std::pair<double, double>>;
            using a_vector = std::vector<large_num>;

            using string_vector= std::pair<std::string, std::string>;
            using vector_2d = std::vector<std::vector<double>>;
            using d_pair = std::pair<double, double>;
            using pair_vector = std::vector<std::pair<double, double>>;
            using vector_d = std::vector<double>;

            struct VectorPack {
                VectorPack(vector_d u, vector_d v) : u(u), v(v) {};

                VectorPack() : u(vector_d()), v(vector_d()) {};
                vector_d u;
                vector_d v;

                void unpack(pair_vector d) {
                    u.reserve(d.size());
                    v.reserve(d.size());
                    for (int j = 0; j < d.size(); ++j) {
                        u[j] = d[j].first;
                        v[j] = d[j].second;
                    }
                }
            };

            /**
             * Able to convert a container to a comma seperated string
             * @tparam Container_T
             * @param begin
             * @param end
             * @return
             */
            template<typename Container_T>
            string str(Container_T begin, Container_T end) {
                stringstream ss;
                bool first = true;
                for (; begin != end; begin++) {
                    double val = *begin;
                    if (!first)
                        ss << ",";
                    if (std::isnan(val)) {
                        ss << "null";
                    } else {
                        ss << std::scientific << setprecision(17) << val;
                    }
                    first = false;
                }
                return ss.str();
            }

            template<typename P>
            std::vector<double> linspace(P start_in, P end_in, int num_in) {

                std::vector<double> linspaced;

                double start = static_cast<double>(start_in);
                double end = static_cast<double>(end_in);
                double num = static_cast<double>(num_in);

                if (num == 0) { return linspaced; }
                if (num == 1) {
                    linspaced.push_back(start);
                    return linspaced;
                }

                double delta = (end - start) / (num - 1);

                for (int i = 0; i < num - 1; ++i) {
                    linspaced.push_back(start + delta * i);
                }
                linspaced.push_back(end);
                return linspaced;
            }

            template<typename InputIterator>
            inline double accumulate(InputIterator first, InputIterator last, int b_max) {
                double init{0.0};
                double boundary{0};
                for (; first != last; ++first) {
                    if (boundary > b_max) { return NAN; }
                    if (std::isnan(*first)) {
                        boundary++;
                        continue;
                    }
                    init = init + *first;
                }
                return init;
            };

            void push(vector_d &u, vector_d &v, string x) {
                u.push_back(::atof(x.c_str()));
                v.push_back(::atof(x.c_str()));
            }

            void push(vector_d &u, vector_d &v, double x) {
                u.push_back(x);
            }

            void push(vector_d &u, vector_d &v, d_pair x) {
                u.push_back(x.first);
                v.push_back(x.second);
            }

            template<typename P>
            void push_nan(vector_d &u, vector_d &v) {/*Genric implementation not valid*/}

            template<>
            void push_nan<d_pair>(vector_d &u, vector_d &v) {
                u.push_back(NAN);
                v.push_back(NAN);
            }

            template<>
            void push_nan<double>(vector_d &u, vector_d &v) {
                u.push_back(NAN);
            }

            /**
             * @brief Downscales a vector to a certain degree
             * @note Currently only 0.5' to 1°
             * Averages the cells
             * @param The continous vector to downscale
             * @return A downscaled represantation of the vector
             */
            vector_d down_scale_v(vector_d &v) {
                vector_d d;
                const int steps{6};
                const int total{36};
                const int boundary{25};
                const int ni{4320};

                int l{1};
                for (int j = 0; j < v.size();) {
                    double val{0};
                    for (int k = 0; k < steps; ++k) {
                        val = val +
                              accumulate(v.begin() + (j + (ni * k)), v.begin() + (j + (steps - 1) + (ni * k)),
                                         boundary);
                    }
                    val = val / total;
                    d.push_back(val);
                    j = j + steps;
                    if (j % ni == 0) {
                        j = (steps * ni) * l;
                        l++;
                    }
                }
                return d;
            }


            /**
             * Convert to double vectors to a vector of double pairs
             * @param v A vector of doubles
             * @param u A vector of doubles
             * @return A vector of pairs of doubles
             */
            pair_vector make_pair_vector(vector_d &v, vector_d &u) {
                assert(v.size() == u.size() && "Size of vectors don't match!");
                pair_vector out;
                for (int j = 0; j < v.size(); ++j) {
                    out.push_back(make_pair(v[j], u[j]));
                }
                return out;
            }

            /**
             * Converts a grid-based vector to a continously spaced vector representation
             * @note Currently only support for 5' resolution
             * @tparam DataType The type of data inside the vector
             * @param data A vector of doubles or double pairs
             * @param p A vector of positions, aligned with the data
             * @param down_scale Should the return vector be downscaled?
             * @return continously spaced vector representation
             */
            template<typename T>
            VectorPack scale(std::vector<T> data, pos_v p, bool down_scale) {
                assert(data.size() == p.size() && "Size of positions and vector don't match!");
                std::deque<T> d_data(data.begin(), data.end());
                std::deque<d_pair> d_pos(p.begin(), p.end());

                double smallest_approach{10};
                double average_dist{0};
                int num = p.size();

                /*
                 * Very simple epsilon function - !do not use to compare floats in general
                 */
                std::function<bool(double, double)> almost_near = [&smallest_approach, &average_dist](const double a,
                                                                                                      const double b) {
                    constexpr double EPSILON_FIVE_MIN{0.06}; //0.0421 = 1/2 5'
                    double dist = std::abs(a - b);
                    bool sm = dist < EPSILON_FIVE_MIN;
                    if (sm) {
                        smallest_approach = dist < smallest_approach ? dist : smallest_approach;
                        average_dist += dist;
                    }
                    return sm;
                };

                std::function<bool(double, double, d_pair)> check = [almost_near](const double x, const double y,
                                                                                  d_pair pos) {
                    return almost_near(pos.first, x) and almost_near(pos.second, y);
                };

                std::function<double(double)> log_mod = [](double x) {
                    return std::signbit(x) ? -std::log(std::abs(x) + 1) : std::log(std::abs(x) + 1);
                };

                vector_d u;
                vector_d v; //only needed if we deal with 2d data like velocities

                /**
                 * vector contains "U-comp,V-comp"
                 * We need to split and fill with null for empty frame
                 * Currently only support for 5' resolution
                 * 360° = 21600' /5 = 4320 W-E
                 * 180° = 10800' /5 = 2160 N-S
                 * Start: 0.0W , 90.0N
                 *
                 * If downscaling enabled:
                 * 1° = 12 * 5'
                 * 12 x 12 -> 144 5' in a 1° cell
                 * 360 W-E
                 * 180 N-S
                 *
                 * 0,5° = 6 * 5'
                 * 6 x 6 = 36 5' in 0.5°
                 * 720 W-E
                 * 360 N-S
                 */
                int ni{4320}; // grid points in W-E
                int nj{2160}; // grid points in N-S

                vector_d pos_y = linspace(90, -90, nj);
                vector_d pos_x = linspace(-180, 180, ni);

                for (int j = 0; j < nj; ++j) {
                    //Moving from N to S
                    for (int k = 0; k < ni; ++k) {
                        //Moving from W-E
                        d_pair pos = d_pos.front();
                        pos.first += -0.042;
                        pos.first += -0.042;
                        if (check(pos_x[k], pos_y[j], pos)) {
                            if (d_data.empty()) { break; }
                            auto d = d_data.front();
                            d_data.pop_front();
                            d_pos.pop_front();
                            push(u, v, d);
                        } else { push_nan<T>(u, v); }
                    }
                }
                LOG(debug) << "Positions lefts unprocessed" << d_pos.size();
                LOG(debug) << "Data left unprocessed" << d_data.size();
                LOG(debug) << "Smallest distance: " << smallest_approach;
                average_dist = average_dist / num;
                LOG(debug) << "Average dist: " << average_dist;

                if (down_scale) {
                    if (u.size() == v.size()) {
                        v = down_scale_v(v);
                    }
                    u = down_scale_v(u);
                }
                if (u.size() == v.size()) {
                    return VectorPack(u, v);
                }
                return VectorPack(u, v);
            }

            /**
             *
             * @param data
             * @param p
             * @return
             */
            std::string scaleDataToString(std::vector<double> data, pos_v p, bool down_scale) {
                VectorPack d = scale(data, p, down_scale);
                return str(d.u.begin(), d.u.end());
            }

            /**
             *
             * @param data
             * @param p
             * @return
             */
            string_vector scaleDataToString(std::vector<d_pair> data, pos_v p, bool down_scale) {
                VectorPack d = scale(data, p, down_scale);
                return std::make_pair(str(d.u.begin(), d.u.end()), str(d.v.begin(), d.v.end()));
            }


            /**
             * @class OutputInterface
             * Writes data to a file
             */
            template<typename T>
            class OutputInterface {
            public:
                virtual ~OutputInterface() {};

                /**
                 * Needs to be implemented
                 * @param filePath
                 * @param printID Bool
                 * @param printXY Bool
                 * @param data Data vector
                 * @param p Position vector
                 */
                virtual void
                write(path filePath, bool printID, bool printXY, std::vector<T> data, pos_v p, a_vector ids) = 0;

            };

            /**
             * @class CSVOutput
             * Writes data to a CSV file
             */
            template<typename T>
            class CSVOutput : public OutputInterface<T> {
            public:
                void
                write(path filePath, bool printID, bool printXY, std::vector<std::pair<double, double>> data, pos_v p,
                      a_vector ids) {}

                void write(path filePath, bool printID, bool printXY, std::vector<bool> data, pos_v p, a_vector ids) {}

                void write(path filePath, bool printID, bool printXY, std::vector<double> data, pos_v p, a_vector ids) {
                    std::vector<std::string> d;
                    d.reserve(data.size());
                    for (int i = 0; i < data.size(); ++i) {
                        std::ostringstream out;
                        out << std::scientific << std::setprecision(17) << data[i];
                        d.emplace_back(out.str());
                    }
                    write(filePath, printID, printXY, d, p, ids);
                }

                void write(path filePath, bool printID, bool printXY, std::vector<std::string> data, pos_v p, a_vector ids) {
                    std::ofstream ofs;
                    ofs.open(filePath + ".csv", std::ofstream::out | std::ofstream::trunc);
                    if (printID) {
                        ofs << "nodeID,";
                    }
                    if (printXY) {
                        ofs << "lat,lon,";
                    }
                    ofs << "data";
                    ofs << std::endl;

                    for (int i = 0; i < data.size(); ++i) {
                        if (printID) { ofs << ids[i] << ","; }
                        if (printXY) { ofs << p[i].first << "," << p[i].second << ","; }

                        ofs << data[i] << std::endl;
                    }
                    ofs.close();
                }
            };

            /**
             * @class NETCDFOutput
             * Writes data to a NETCDF file
             */
            template<typename T>
            class NETCDFOutput : public OutputInterface<T> {
            public:

                void
                write(path filePath, bool printID, bool printXY, std::vector<std::pair<double, double>> data, pos_v p,
                      a_vector ids) {
                    LOG(userinfo) << "not implemented yet";
                }

                void write(path filePath, bool printID, bool printXY, std::vector<bool> data, pos_v p, a_vector ids) {
                    LOG(userinfo) << "not implemented yet";
                }

                void
                write(path filePath, bool printID, bool printXY, std::vector<std::string> data, pos_v p, a_vector ids) {
                    LOG(userinfo) << "not implemented yet";
                }

                void
                write(path filePath, bool printID, bool printXY, std::vector<double> geo_data, pos_v p, a_vector ids) {
                /**

		    try {
	                        netCDF::NcFile dataFile(filePath + ".nc", netCDF::NcFile::replace);

                        //Get continous data
                        auto data = scale(geo_data, p, false);

                        VectorPack positions = makePositions();

                        const int x_num{static_cast<int>(positions.u.size())};
                        const int y_num{static_cast<int>(positions.v.size())};

                        // Create netCDF dimensions
                        netCDF::NcDim xDim = dataFile.addDim("latitude", x_num);
                        netCDF::NcDim yDim = dataFile.addDim("longitude", y_num);
                        netCDF::NcVar latVar = dataFile.addVar("latitude", netCDF::ncFloat, xDim);
                        netCDF::NcVar lonVar = dataFile.addVar("longitude", netCDF::ncFloat, yDim);

                        const double *u = positions.u.data();
                        const double *v = positions.v.data();
                        latVar.putVar(u);
                        lonVar.putVar(v);
                        latVar.putAtt("units", "degrees_north");
                        lonVar.putAtt("units", "degrees_east");

                        std::vector<netCDF::NcDim> dims;
                        dims.push_back(xDim);
                        dims.push_back(yDim);
                        netCDF::NcVar nc_data = dataFile.addVar("data", netCDF::ncDouble, dims);
                        nc_data.putAtt("units", "m");

                        // Write the data to the file. Although netCDF supports
                        // reading and writing subsets of data, in this case we write all
                        // the data in one operation.
                        std::vector<std::vector<double>> a = get2darray(data.u);
                        double **data_p = vectorToArray(a);
                        nc_data.putVar(data_p);
                        dataFile.sync();
                        dataFile.close();
                        for (int i = 0; i < a.size(); ++i) { delete[] data_p[i]; }
                        delete[] data_p;
                    } catch (netCDF::exceptions::NcException &e) {
                        LOG(critical) << "NetCDF exception: " << e.what();
                    }
		**/
                }

            private:

                const int ni{4320}; // grid points in W-E
                const int nj{2160}; // grid points in N-S

                VectorPack makePositions() {
                    VectorPack out;
                    out.u = linspace(90, -90, nj);
                    out.v = linspace(-180, 180, ni);
                    return out;
                }

                std::vector<std::vector<double>> get2darray(std::vector<double> data) {
                    std::vector<std::vector<double>> out(nj, std::vector<double>(ni));
                    for (int j = 0; j < nj; ++j) {
                        //Moving from N to S
                        for (int k = 0; k < ni; ++k) {
                            //Moving from W-E
                            out[j][k] = data[(j * ni) + k];
                        }
                    }
                    return out;
                }

                double **vectorToArray(vector <vector<double>> &vals) {
                    int N = vals.size();
                    int M = vals[0].size();
                    double **whereto = new double *[N];
                    for (unsigned i = 0; (i < N); i++) {
                        whereto[i] = new double[M];
                        for (unsigned j = 0; (j < M); j++) {
                            whereto[i][j] = vals[i][j];
                        }
                    }
                    return whereto;
                }

            };

            /**
            * JSON-ified GRIB files. Example:
            *
            *     [
            *       {
            *         "header": {
            *           "refTime": "2013-11-30T18:00:00.000Z",
            *           "parameterCategory": 2,
            *           "parameterNumber": 2,
            *           "surface1Type": 100,
            *           "surface1Value": 100000.0,
            *           "forecastTime": 6,
            *           "scanMode": 0,
            *           "nx": 360,
            *           "ny": 181,
            *           "lo1": 0,
            *           "la1": 90,
            *           "lo2": 359,
            *           "la2": -90,
            *           "dx": 1,
            *           "dy": 1
            *         },
            *         "data": [3.42, 3.31, 3.19, 3.08, 2.96, 2.84, 2.72, 2.6, 2.47, ...]
            *       }
            *     ]
            *
            * TODO rebuild me with boost ptree
            */


            /**
             * @class Write data to a pseudo GFS json-like file
             * @note provides additonal methods for downscaling
             */
            template<typename T>
            class GFS_JSONOutput : public OutputInterface<T> {
            private:

                string_vector buildData(std::vector<double> data, pos_v p, std::ofstream &ofs) {
                    std::string v = scaleDataToString(data, p, true);
                    ofs << "[\n" << buildHeader() << v << "]\n}]";
                    return string_vector();
                }

                /**
                 * Scan mode 0 assumed. Longitude increases from λ0, and latitude decreases from φ0.
                 * http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-4.shtml
                 * @return a string of an array containing data in steps
                 *
                 * First data of U-component then data of V-component
                 * TODO support non vector data
                 */
                string_vector buildData(std::vector<d_pair> data, pos_v p, std::ofstream &ofs) {
                    assert(data.size() == p.size() && "Size of positions and vector don't match!");
                    std::deque<d_pair> d_data(data.begin(), data.end());
                    std::deque<d_pair> d_pos(p.begin(), p.end());

                    string_vector v = scaleDataToString(data, p, true);
                    ofs << "[\n" << buildHeader() << v.first << "]\n}," << buildHeader() << v.second << "]\n}]";
                    return v;
                }

                string_vector buildData(std::vector<bool> data, pos_v p, std::ofstream &ofs) { return string_vector(); }

                string_vector
                buildData(std::vector<std::string> data, pos_v p, std::ofstream &ofs) { return string_vector(); }

                std::string buildHeader() {
                    std::time_t result = std::time(nullptr);
                    std::string t = std::asctime(std::localtime(&result));
                    t.pop_back();
                    const std::string time =
                            "\"refTime\":\"" + t + "\",\n" + "\"forecastTime\": 0,";
                    //nx, ny  number of grid points W-E and N-S
                    //const std::string params = "\"lo1\": -180, \"la1\": 90, \"dx\":0.0833333, \"dy\":0.0833333, \"nx\": 4320, \"ny\":2160";
                    //const std::string params = "\"lo1\": -180, \"la1\": 90, \"dx\":1, \"dy\":1, \"nx\": 360, \"ny\":180";
                    const std::string params = "\"lo1\": -180, \"la1\": 90, \"dx\":0.5, \"dy\":0.5, \"nx\": 720, \"ny\":360";
                    const std::string header = "{\"header\": {\n" + time + params + "}, \"data\":[";
                    return header;
                }

            public:
                void write(path filePath, bool printID, bool printXY, std::vector<T> data, pos_v p, a_vector ids) {
                    std::ofstream ofs;
                    ofs.open(filePath + ".json", std::ofstream::out | std::ofstream::trunc);
                    buildData(data, p, ofs);
                    ofs.close();
                }
            };

            /**
             * @class Factory for yielding correct data output file
             */
            template<typename T>
            class OutputFactory {
            private:
                static OutputInterface<T> *create(OutputType type) {
                    switch (type) {
                        case OutputType::CSV:
                            return (OutputInterface<T> *) new CSVOutput<T>();
                        case OutputType::NET_CDF :
                            return (OutputInterface<T> *) new NETCDFOutput<T>();
                        case OutputType::GFS_JSON :
                            return (OutputInterface<T> *) new GFS_JSONOutput<T>();
                        case OutputType::NON_VALID:
                            throw std::bad_function_call();
                    }
                }

            public:
                static OutputInterface<T> *getOutput(OutputType type) {
                    return create(type);
                };
            };

        }
    }
}
#endif
