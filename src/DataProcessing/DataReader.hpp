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

#ifndef GLOBAL_FLOW_DATAREADER_HPP
#define GLOBAL_FLOW_DATAREADER_HPP

#include "../../lib/csv.h"
#include "../Model/Node.hpp"
#include "../Logging/Logging.hpp"
#include "../Simulation/Options.hpp"
#include <boost/filesystem.hpp>
#include <fstream>
#include <vector>


namespace GlobalFlow {
using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<Model::NodeInterface>>>;
namespace fs = boost::filesystem;
using large_num = unsigned long int;

/**
 * @interface Interface that needs to be implemented for reading in required data for the model
 */
class DataReader {
    protected:
        NodeVector nodes;
        int stepMod;
        //<GlobalID, ID>
        std::unordered_map<int, int> lookupglobIDtoID;
        //<ARCID(0.5°), vector<GlobalID(5')>>
        std::unordered_map<int, std::vector<int>> lookupZeroPointFivetoFiveMinute;
    public:
        virtual ~DataReader() {}

        void initNodes(NodeVector nodes) { this->nodes = nodes; }

        /**
         * @attention Needs to be implemented!
         * @note Is called by simulation at startup
         * @param op Options object
         */
        virtual void readData(Simulation::Options op) = 0;

        /**
         * Gereric method for looping through files inside a directory and apply a gernic function
         * @param path
         * @param files
         * @param fun a function
         */
        template<class Fun>
        void loopFiles(std::string path, std::vector<std::string> files, Fun fun) {
            for (std::string file : files) {
                std::string real_path = path + '/' + file;
                fun(real_path);
            }
        }

        inline int check(int globid) {
            int i{0};
            try {
                i = lookupglobIDtoID.at(globid);
            }
            catch (const std::out_of_range &ex) {
                return -1;
            }
            return i;
        }

        /**
         * Read data from a two-column csv file and apply function to data
         * @param path
         * @param processData A processing function e.g. upscaling of data
         */
        template<class ProcessDataFunction>
        void readTwoColumns(std::string path, ProcessDataFunction processData) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "global_ID", "data");
            int globid = 0;
            double data = 0;
            int pos = 0;
            while (in.read_row(globid, data)) {
                pos = check(globid);
                if (pos == -1) {
                    continue;
                }
                if (nodes->at(pos)->getProperties().get<large_num, Model::ArcID>() != globid) {
                    throw "Error in reading globID";
                }
                processData(data, pos);
            }
        }

        /**
         * Creates a mapping of 0.5° ArcIDs to a list of contained 5' GlobIDs
         * @param path
         */
        void readZeroPointFiveToFiveMin(std::string path) {
            io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>> in(path);
            in.read_header(io::ignore_no_column, "GLOBALID", "ARC_ID");
            int globID = 0;
            int arcID = 0;

            while (in.read_row(globID, arcID)) {
                lookupZeroPointFivetoFiveMinute[arcID].push_back(std::move(globID));
            }
        }


        const std::unordered_map<int, std::vector<int>> &
        getArcIDMapping() {
            return lookupZeroPointFivetoFiveMinute;
        };

        const std::unordered_map<int, int> &
        getGlobIDMapping() {
            return lookupglobIDtoID;
        };

        fs::path data_dir{"data"};

        std::string buildDir(std::string path) {
            return (data_dir / fs::path(path)).string();
        };
};
}
#endif //DATAREADER_HPP
