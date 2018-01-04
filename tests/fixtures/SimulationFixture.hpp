//
// Created by robert on 21.01.16.
//

#ifndef STANDALONE_DATAREADERFIXTURE_HPP
#define STANDALONE_DATAREADERFIXTURE_HPP

#include <gtest/gtest.h>
#include <fstream>
#include "../../src/Node.hpp"
#include "../../src/Equation.hpp"
#include <unordered_set>

using namespace std;

class SimulationFixture : public ::testing::Test {
public:
    shared_ptr<vector<std::unique_ptr<NodeInterface>>> nodes;
    size_t numberOfNodes = 300;
    size_t numberOfNodes2 = 1200;

    void SetUp(){
        shared_ptr<vector<std::unique_ptr<NodeInterface>>> ptr(new vector<std::unique_ptr<NodeInterface>>);
        ptr.get()->reserve(numberOfNodes);
        nodes.swap(ptr);
    };

    void reset(){
        shared_ptr<vector<std::unique_ptr<NodeInterface>>> ptr(new vector<std::unique_ptr<NodeInterface>>);
        ptr.get()->reserve(numberOfNodes2);
        nodes.swap(ptr);
    }

    void TearDown(){

    };

    int Test__readNodesFromFile(){
        return readNodesFromFile("tests/files/grid_oneLayer.txt",1,15,20);
    };

    int Test__readNodesLayerFromFile(){
        reset();
        return readNodesFromFile("tests/files/grid_multipleLayers.txt",4,15,20);
    }

    int Test__runSimpleWithrivers(){
        reset();
        readNodesFromFile2("tests/files/simple.txt",1,10,11);
        Options op = Options();
        op.load("tests/files/test_config_simple.json");

        Equation *eq = new Equation(110, nodes, op);

        eq->toogleSteadyState();
        eq->solve();

        quantity<CubicMeter> out = 0;
        quantity<CubicMeter> in = 0;
        quantity<CubicMeter> storage = 0;
        quantity<CubicMeter> totalConstantHeadFlow = 0;
        double wellFlow = 0;
        for (int j = 0; j < nodes.get()->size(); ++j) {
            if(nodes.get()->at(j)->getArea() != 1 * si::square_meter){
                out = out + nodes.get()->at(j).get()->getOUT();
                in = in + nodes.get()->at(j).get()->getIN();
                storage = storage + nodes->at(j)->getTotalStorageFlow();
                totalConstantHeadFlow = totalConstantHeadFlow + nodes->at(j)->getConstantHeadFlow();
                wellFlow = wellFlow + nodes->at(j)->getExternlFlowVolumeByName(RECHARGE);
            }
        }

        cerr << "Error: "  << eq->getError() << "\n";
        cerr << "Iterrations: " << eq->getItter() << "\n";
        cerr << "IN " << in.value() << "\n";
        cerr << "OUT " << out.value() << "\n";
        cerr << "TotalStorageFlow" << storage.value() << "\n";
        cerr << "ConstantHeadFlow" << totalConstantHeadFlow.value() << "\n";
        cerr << "WellFlow" << wellFlow << "\n";

        int row = 0;
        for(const auto& node: *nodes){
            if(node.get()){
                cerr << std::setprecision (6) <<  node.get()->getHead().value() << "  ";
            }
            else{
                cerr << "0 ";
            }
            if(row == 9){
                cerr << "\n";
                row = 0;
            } else
                row++;
        }

        return 0;
    }

    int Test__runSimpleWithrivers3D(){
        reset();
        readNodesFromFile3("tests/files/simple_3D.txt",2,10,11);
        Options op = Options();
        op.load("tests/files/test_config_simple.json");

        Equation *eq = new Equation(220, nodes, op);

        eq->toogleSteadyState();
        eq->solve();

        quantity<CubicMeter> out = 0;
        quantity<CubicMeter> in = 0;
        quantity<CubicMeter> storage = 0;
        quantity<CubicMeter> totalConstantHeadFlow = 0;
        double wellFlow = 0;
        for (int j = 0; j < nodes.get()->size(); ++j) {
            if(nodes.get()->at(j)->getArea() != 1 * si::square_meter){
                out = out + nodes.get()->at(j).get()->getOUT();
                in = in + nodes.get()->at(j).get()->getIN();
                storage = storage + nodes->at(j)->getTotalStorageFlow();
                totalConstantHeadFlow = totalConstantHeadFlow + nodes->at(j)->getConstantHeadFlow();
                wellFlow = wellFlow + nodes->at(j)->getExternlFlowVolumeByName(RECHARGE);
            }
        }

        cerr << "Error: "  << eq->getError() << "\n";
        cerr << "Iterrations: " << eq->getItter() << "\n";
        cerr << "IN " << in.value() << "\n";
        cerr << "OUT " << out.value() << "\n";
        cerr << "TotalStorageFlow" << storage.value() << "\n";
        cerr << "ConstantHeadFlow" << totalConstantHeadFlow.value() << "\n";
        cerr << "WellFlow" << wellFlow << "\n";

        int row = 0;
        for(const auto& node: *nodes){
            if(node.get()){
                cerr << std::setprecision (6) <<  node.get()->getHead().value() << "  ";
            }
            else{
                cerr << "0 ";
            }
            if(row == 9){
                cerr << "\n";
                row = 0;
            } else
                row++;
        }

        return 0;
    }

    int Test__runEquationStep(){

        reset();
        readNodesFromFile("tests/files/grid_multipleLayers.txt",4,15,20);

        //1) make equation
        Options op = Options();
        op.load("tests/files/test_config.json");
        //210 * 4
        Equation *eq = new Equation(1200, nodes, op);
        //2) run step

        double row = 0;
        /**
        for(const auto& node: *nodes){
            if(node.get()){
                cerr << std::setprecision (6) <<  node.get()->getElevation().value() << "  ";
            }
            else{
                cerr << "0 ";
            }
            if(row == 14){
                cerr << "\n";
                row = 0;
            } else
                row++;
        }
         */

        //Stress 1
        eq->toogleSteadyState();
        eq->solve();


        //Stress 2
        //eq->toogleSteadyState();
        //for (int k = 0; k < 60; ++k) {
        //    eq->solve();
        //}

        //Stress 3
        //for (int k = 0; k < 60; ++k) {
        //    eq->solve();
        //}
        //3) check results
        quantity<CubicMeter> out = 0;
        quantity<CubicMeter> in = 0;
        quantity<CubicMeter> storage = 0;
        quantity<CubicMeter> totalConstantHeadFlow = 0;
        double wellFlow = 0;
        for (int j = 0; j < nodes.get()->size(); ++j) {
            if(nodes.get()->at(j)->getArea() != 1 * si::square_meter){
                out = out + nodes.get()->at(j).get()->getOUT();
                in = in + nodes.get()->at(j).get()->getIN();
                storage = storage + nodes->at(j)->getTotalStorageFlow();
                totalConstantHeadFlow = totalConstantHeadFlow + nodes->at(j)->getConstantHeadFlow();
                wellFlow = wellFlow + nodes->at(j)->getExternlFlowVolumeByName(RECHARGE);
            }
        }

        cerr << "Error: "  << eq->getError() << "\n";
        cerr << "Iterrations: " << eq->getItter() << "\n";
        cerr << "IN " << in.value() << "\n";
        cerr << "OUT " << out.value() << "\n";
        cerr << "TotalStorageFlow" << storage.value() << "\n";
        cerr << "ConstantHeadFlow" << totalConstantHeadFlow.value() << "\n";
        cerr << "WellFlow" << wellFlow << "\n";

        row = 0;
        for(const auto& node: *nodes){
            if(node.get()){
                cerr << std::setprecision (6) <<  node.get()->getHead().value() << "  ";
            }
            else{
                cerr << "0 ";
            }
            if(row == 14){
                cerr << "\n";
                row = 0;
            } else
                row++;
        }
        //cerr << "------------------------------------------------------------------------";

        //for(const auto& node: *nodes){
        //    cerr << node.get() << "\n";
        //}


        return 0;
    }

private:
    //Small case
    double initalHead_s = 0;
    double elevation_s = 0;
    double depth_s = 100;
    double hydroConRows_s = 1.550000e-04;
    double vertCon_s = 1.550000E-05;


    //Big case
    double initialHead = 100;
    double specificSorage = 1.400000e-6;
    double yield = 0.3000000;
    double hydConRows = 4.0000;
    double hydConVertial = 0.400000;

    //Only Confined
    double hydConRows2 = 0.3;
    double hydConVertial2 = 4.00000;
    double specificStorage2 = 0.4;

    double hydConRows3 = 1.4e-6;
    double hydConVertial3 = 0.3;
    double specificStorage3 = 1.0e-2;

    double hydroCondRows4 = 1.0e-2;
    double hydConVetical4 = 1.4e-6;
    double specficStorage4 = 0.3;


    //double hydConRows3 = 1.0e-02;
    //double hydConVertial3 = 1.0e-02;

    inline int calcPosition(int row, int column){
        int row_lenght = 15;
        return ((row - 1) * row_lenght) + (column -1);
    }

    inline int calcPosition2(int row, int column){
        int row_lenght = 10;
        return ((row - 1) * row_lenght) + (column -1);
    }

    std::unordered_set<int> wells;
    std::set<pair<int,int>> wells_def = {{1,8},
                                     {2,5},
                                     {2,8},
                                     {2,12},
                                     {3,4},
                                     {4,12},
                                     {5,3},
                                     {5,13},
                                     {6,14},
                                     {7,2},
                                     {14,2},
                                     {14,14},
                                     {16,3},
                                     {16,13},
                                     {17,13},
                                     {18,4},
                                     {18,12},
                                     {19,7}};

    std::unordered_set<int> rivers;
    std::set<pair<int,int>> rivers_def = {{1,10},
                                         {2,10},
                                         {3,10},
                                         {4,10},
                                         {5,10},
                                         {6,10},
                                         {7,10},
                                         {8,10},
                                         {9,10},
                                         {10,10},
                                         {11,10}};

    std::unordered_set<int> rivers2;
    std::set<pair<int,int>> rivers_def2 = {{1,10},
                                          {2,10},
                                          {3,10},
                                          {4,10},
                                          {5,10},
                                          {6,10},
                                          {7,10},
                                          {8,10},
                                          {9,10},
                                          {10,10},
                                          {11,10},
                                          {12,10},
                                          {13,10},
                                          {14,10},
                                          {15,10},
                                          {16,10},
                                          {17,10},
                                          {18,10},
                                          {19,10},
                                          {20,10},
                                          {21,10},
                                          {22,10}};

    std::unordered_set<int> ghb;
    std::set<pair<int,int>> ghb_def = {{1,1},
                                       {2,1},
                                       {3,1},
                                       {4,1},
                                       {5,1},
                                       {6,1},
                                       {7,1},
                                       {8,1},
                                       {9,1},
                                       {10,1},
                                       {11,1},
    };

    std::unordered_set<int> ghb2;
    std::set<pair<int,int>> ghb_def2 = {{1,1},
                                       {2,1},
                                       {3,1},
                                       {4,1},
                                       {5,1},
                                       {6,1},
                                       {7,1},
                                       {8,1},
                                       {9,1},
                                       {10,1},
                                       {11,1},
                                       {12,1},
                                       {13,1},
                                       {14,1},
                                       {15,1},
                                       {16,1},
                                       {17,1},
                                       {18,1},
                                       {19,1},
                                       {20,1},
                                       {21,1},
                                       {22,1}};

    /**
    * Read defintions from MODFLOW like grid file
    */
    size_t readNodesFromFile(string path, size_t layers, int col, int row){
        ifstream in(path);
        size_t i = 0;

        if(!in){
            cerr << "Cannot open file \n";
            return 0;
        }

        for(auto pe : wells_def)
            wells.insert(calcPosition(pe.first,pe.second));

        double K = 0;
        int layer = 0;
        double value;
        vector<int> no_flow;
        for(size_t y = 0; y < row * layers + layers-1 ; y++){
            for (size_t x = 0; x < col; ++x) {
                in >> value;
                if(value == -99){
                    //new layer
                    layer++;
                    break;
                }
                if(value == -1){
                    nodes.get()->emplace_back(new StaticHeadNode(nodes, i, 4000000 * si::square_meter));
                    if(layer == 0){
                        nodes.get()->at(i)->setK_direct(hydConRows * si::meter / day);
                        nodes.get()->at(i)->setK_vertical(hydConVertial * si::meter / day);
                        nodes.get()->at(i)->setSStorage(specificSorage * perMeter);
                    }
                    if(layer == 1){
                        nodes.get()->at(i)->setK_direct(hydConRows2 * si::meter / day);
                        nodes.get()->at(i)->setK_vertical(hydConVertial2 * si::meter / day);
                        nodes.get()->at(i)->setSStorage(specificStorage2 * perMeter);
                    }
                    if(layer == 2){
                        nodes.get()->at(i)->setK_direct(hydConRows3 * si::meter / day);
                        nodes.get()->at(i)->setK_vertical(hydConVertial3 * si::meter / day);
                        nodes.get()->at(i)->setSStorage(specificStorage3 * perMeter);
                    }
                    if(layer == 3){
                        nodes.get()->at(i)->setK_direct(hydroCondRows4 * si::meter / day);
                        nodes.get()->at(i)->setK_vertical(hydConVetical4 * si::meter / day);
                        nodes.get()->at(i)->setSStorage(specficStorage4 * perMeter);
                    }
                    nodes.get()->at(i)->setHead_direct(initialHead);
                    nodes.get()->at(i)->setYield(yield);
                    nodes.get()->at(i)->setSimpleDistance();
                    nodes.get()->at(i)->setSimpleK();
                    switch (layer){
                        case 0 :
                            nodes.get()->at(i)->setElevation__simple(150 * si::meter);
                            nodes.get()->at(i)->setDepth(100 * si::meter);
                            break;
                        case 1 :
                            nodes.get()->at(i)->setElevation__simple(50 * si::meter);
                            nodes.get()->at(i)->setDepth(150 * si::meter);
                            break;
                        case 2:
                            nodes.get()->at(i)->setElevation__simple(-100 * si::meter);
                            nodes.get()->at(i)->setDepth(50 * si::meter);
                            break;
                        case 3:
                            nodes.get()->at(i)->setElevation__simple(-150 * si::meter);
                            nodes.get()->at(i)->setDepth(200 * si::meter);
                            break;
                        default:
                            break;
                    }
                    i++;
                    continue;
                }
                if(value == 0){
                    nodes.get()->emplace_back(new StaticHeadNode(nodes, i, 1 * si::square_meter));
                    nodes.get()->at(i)->setHead_direct(0);
                    i++;
                    continue;
                }

                K = hydConRows;
                nodes.get()->emplace_back(new StandardNode(nodes,x,y,4000000 * si::square_meters,i,i,K * (si::meter / day)));

                if(layer == 0){
                    nodes.get()->at(i)->setK_direct(hydConRows * si::meter / day);
                    nodes.get()->at(i)->setK_vertical(hydConVertial * si::meter / day);
                    nodes.get()->at(i)->setSStorage(specificSorage * perMeter);
                }
                if(layer == 1){
                    nodes.get()->at(i)->setK_direct(hydConRows2 * si::meter / day);
                    nodes.get()->at(i)->setK_vertical(hydConVertial2 * si::meter / day);
                    nodes.get()->at(i)->setSStorage(specificStorage2 * perMeter);
                }
                if(layer == 2){
                    nodes.get()->at(i)->setK_direct(hydConRows3 * si::meter / day);
                    nodes.get()->at(i)->setK_vertical(hydConVertial3 * si::meter / day);
                    nodes.get()->at(i)->setSStorage(specificStorage3 * perMeter);
                }
                if(layer == 3){
                    nodes.get()->at(i)->setK_direct(hydroCondRows4 * si::meter / day);
                    nodes.get()->at(i)->setK_vertical(hydConVetical4 * si::meter / day);
                    nodes.get()->at(i)->setSStorage(specficStorage4 * perMeter);
                }
                nodes.get()->at(i)->setHead_direct(initialHead);
                nodes.get()->at(i)->setYield(yield);
                nodes.get()->at(i)->setSimpleDistance();
                nodes.get()->at(i)->setSimpleK();

                switch (layer){
                    case 0 :
                        nodes.get()->at(i)->setElevation__simple(150 * si::meter);
                        nodes.get()->at(i)->setDepth(100 * si::meter);
                        break;
                    case 1 :
                        nodes.get()->at(i)->setElevation__simple(50 * si::meter);
                        nodes.get()->at(i)->setDepth(150 * si::meter);
                        break;
                    case 2:
                        nodes.get()->at(i)->setElevation__simple(-100 * si::meter);
                        nodes.get()->at(i)->setDepth(50 * si::meter);
                        break;
                    case 3:
                        nodes.get()->at(i)->setElevation__simple(-150 * si::meter);
                        nodes.get()->at(i)->setDepth(200 * si::meter);
                        break;
                    default:
                        break;
                }
                i++;
            }
        }

        auto checkNoFlow = [this](int place) -> bool{
            try {
                if(nodes->at(place)->getArea() == 1* si::square_meter){
                    return false;
                } else
                    return true;
            }catch (...){
                return false;
            }
        };


        for (int l = 0; l < layers; ++l) {
            for (int r = 0; r < row; ++r) {
                for (int c = 0; c < col; ++c) {
                    double currentPos = l*(row * col) + col*r + c;
                    if(nodes->at(currentPos)->getArea() == 1* si::square_meter)
                        continue;
                    //LEFT
                    if(c > 0 and checkNoFlow(currentPos - 1))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - 1,LEFT);
                    //RIGHT
                    if(c + 1 < col and checkNoFlow(currentPos + 1))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + 1,RIGHT);
                    //FRONT
                    if(r > 0 and checkNoFlow(currentPos - col))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - col,FRONT);
                    //BACK
                    if(r + 1 < row and checkNoFlow(currentPos + col))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + col,BACK);
                    //TOP
                    if(l > 0 and checkNoFlow((currentPos - (col * row))))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - (col * row),TOP);
                    //DOWN
                    if(l + 1 < layers and checkNoFlow(currentPos + (col * row)))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + (col * row),DOWN);

                    const bool is_in = wells.find(currentPos) != wells.end();
                    if(is_in){
                        nodes.get()->at(currentPos)->addExternalFlow(RECHARGE,0 * si::meter,2200,0);
                    }
                }
            }
        }

        return i;
    }

    /**
    * Read defintions from MODFLOW like grid file
     * FIXME !!!! Redundant code!!
    */
    size_t readNodesFromFile2(string path, size_t layers, int col, int row){
        ifstream in(path);
        size_t i = 0;

        if(!in){
            cerr << "Cannot open file \n";
            return 0;
        }

        for(auto pe : rivers_def)
            rivers.insert(calcPosition2(pe.first,pe.second));

        for(auto pe : ghb_def)
            ghb.insert(calcPosition2(pe.first,pe.second));

        double K = 0;
        int layer = 0;
        double value;
        vector<int> no_flow;
        for(size_t y = 0; y < row * layers + layers-1 ; y++){
            for (size_t x = 0; x < col; ++x) {
                in >> value;


                nodes->emplace_back(new StandardNode(nodes,x,y,1000 * 1000 * si::square_meters,i,i,hydroConRows_s * (si::meter / day)));
                nodes->at(i)->setK_direct(hydroConRows_s * si::meter / day);
                nodes->at(i)->setK_vertical(vertCon_s * si::meter / day);

                //TODO
                nodes->at(i)->setSStorage(specificSorage * perMeter);

                nodes->at(i)->setHead_direct(initalHead_s);

                //TODO
                nodes->at(i)->setYield(yield);

                nodes->at(i)->setSimpleDistance();
                nodes->at(i)->setSimpleK();
                nodes->at(i)->setElevation__simple(0 * si::meter);
                nodes->at(i)->setDepth(100 * si::meter);
                i++;
            }
        }

        auto checkNoFlow = [this](int place) -> bool{
          try {
              if(nodes->at(place)->getArea() == 1* si::square_meter){
                  return false;
              } else
                  return true;
          }catch (...){
              return false;
          }
        };


        for (int l = 0; l < layers; ++l) {
            for (int r = 0; r < row; ++r) {
                for (int c = 0; c < col; ++c) {
                    double currentPos = l*(row * col) + col*r + c;
                    if(nodes->at(currentPos)->getArea() == 1* si::square_meter)
                        continue;
                    //LEFT
                    if(c > 0 and checkNoFlow(currentPos - 1))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - 1,LEFT);
                    //RIGHT
                    if(c + 1 < col and checkNoFlow(currentPos + 1))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + 1,RIGHT);
                    //FRONT
                    if(r > 0 and checkNoFlow(currentPos - col))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - col,FRONT);
                    //BACK
                    if(r + 1 < row and checkNoFlow(currentPos + col))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + col,BACK);
                    //TOP
                    if(l > 0 and checkNoFlow((currentPos - (col * row))))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - (col * row),TOP);
                    //DOWN
                    if(l + 1 < layers and checkNoFlow(currentPos + (col * row)))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + (col * row),DOWN);

                    //stage = 0 conduc= 0.2818E-01 bot-el=-100
                    const bool is_in = rivers.find(currentPos) != rivers.end();
                    if(is_in){
                        nodes->at(currentPos)->addExternalFlow(RIVER,0 * si::meter,0.02818,-100 * si::meter);
                    }
                    //Stage 200 conduc= 0.2818E-01
                    const bool is_in2 = ghb.find(currentPos) != ghb.end();
                    if(is_in2){
                        nodes->at(currentPos)->addExternalFlow(GENERAL_HEAD_BOUNDARY,200 * si::meter,0.02818,0 * si::meter);
                    }
                }
            }
        }

        return i;
    }

    /**
    * Read defintions from MODFLOW like grid file
     * FIXME !!!! Redundant code!!
    */
    size_t readNodesFromFile3(string path, size_t layers, int col, int row){
        ifstream in(path);
        size_t i = 0;

        if(!in){
            cerr << "Cannot open file \n";
            return 0;
        }

        for(auto pe : rivers_def2)
            rivers2.insert(calcPosition2(pe.first,pe.second));

        for(auto pe : ghb_def2)
            ghb2.insert(calcPosition2(pe.first,pe.second));

        double K = 0;
        int layer = 0;
        double value;
        vector<int> no_flow;
        for(size_t y = 0; y < row * layers + layers-1 ; y++){
            for (size_t x = 0; x < col; ++x) {
                in >> value;

                if(value == -99){
                    //new layer
                    layer++;
                    break;
                }

                nodes->emplace_back(new StandardNode(nodes,x,y,1000 * 1000 * si::square_meters,i,i,hydroConRows_s * (si::meter / day)));
                nodes->at(i)->setK_direct(hydroConRows_s * si::meter / day);
                nodes->at(i)->setK_vertical(vertCon_s * si::meter / day);

                //TODO
                nodes->at(i)->setSStorage(specificSorage * perMeter);

                nodes->at(i)->setHead_direct(initalHead_s);

                //TODO
                nodes->at(i)->setYield(yield);

                nodes->at(i)->setSimpleDistance();
                nodes->at(i)->setSimpleK();
                if(layer == 0)
                    nodes->at(i)->setElevation__simple(0 * si::meter);
                if(layer == 1)
                    nodes->at(i)->setElevation__simple(-50 * si::meter);
                nodes->at(i)->setDepth(50 * si::meter);
                i++;
            }
        }

        auto checkNoFlow = [this](int place) -> bool{
          try {
              if(nodes->at(place)->getArea() == 1* si::square_meter){
                  return false;
              } else
                  return true;
          }catch (...){
              return false;
          }
        };


        for (int l = 0; l < layers; ++l) {
            for (int r = 0; r < row; ++r) {
                for (int c = 0; c < col; ++c) {
                    double currentPos = l*(row * col) + col*r + c;
                    if(nodes->at(currentPos)->getArea() == 1* si::square_meter)
                        continue;
                    //LEFT
                    if(c > 0 and checkNoFlow(currentPos - 1))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - 1,LEFT);
                    //RIGHT
                    if(c + 1 < col and checkNoFlow(currentPos + 1))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + 1,RIGHT);
                    //FRONT
                    if(r > 0 and checkNoFlow(currentPos - col))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - col,FRONT);
                    //BACK
                    if(r + 1 < row and checkNoFlow(currentPos + col))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + col,BACK);
                    //TOP
                    if(l > 0 and checkNoFlow((currentPos - (col * row))))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos - (col * row),TOP);
                    //DOWN
                    if(l + 1 < layers and checkNoFlow(currentPos + (col * row)))
                        nodes.get()->at(currentPos).get()->setNeighbour(currentPos + (col * row),DOWN);

                    const bool is_in = rivers2.find(currentPos) != rivers2.end();
                    if(is_in){
                        nodes->at(currentPos)->addExternalFlow(RIVER,0 * si::meter,0.01409,-50 * si::meter);
                    }
                    //Stage 200 conduc= 0.2818E-01
                    const bool is_in2 = ghb2.find(currentPos) != ghb2.end();
                    if(is_in2){
                        nodes->at(currentPos)->addExternalFlow(GENERAL_HEAD_BOUNDARY,200 * si::meter,0.01409,0 * si::meter);
                    }
                }
            }
        }

        return i;
    }
};


#endif //STANDALONE_DATAREADERFIXTURE_HPP