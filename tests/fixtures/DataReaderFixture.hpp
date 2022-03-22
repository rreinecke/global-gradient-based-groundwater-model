//
// Created by robert on 21.01.16.
//

#ifndef DATAREADERFIXTURE_HPP
#define DATAREADERFIXTURE_HPP

#include <gtest/gtest.h>
#include "../../src/Node.hpp"
#include "../../src/DataReader.hpp"

using namespace std;

class DataReaderFixture : public ::testing::Test {
public:
    NodeVector nodes;

    void SetUp(){
    };
    void TearDown(){
    };

    int buildNeighbourMap__TEST(){

        NodeVector ptr(new vector<unique_ptr<NodeInterface>>);
        nodes = std::move(ptr);

        DataReader reader;
        int numberOfTOPNodes = 12;
        nodes->reserve(numberOfTOPNodes);

        reader.readLandMask(nodes,"tests/files/motherDefSIMPLE.csv",12,0);
        reader.buildNeighbourMap(nodes, numberOfTOPNodes, 1);

        int out = 0;
        for(const auto& node : *nodes){
            if(node->isStaticNode())
                out++;
        }
        //return number of static nodes
        return out;
    }

    int buildNeighbourBigMap__TEST(){

      NodeVector ptr(new vector<unique_ptr<NodeInterface>>);
      nodes = std::move(ptr);

      DataReader reader;
      int numberOfTOPNodes = 22;
      nodes->reserve(numberOfTOPNodes);

      reader.readLandMask(nodes,"tests/files/motherDefSIMPLE2.csv",numberOfTOPNodes,0);
      reader.buildNeighbourMap(nodes, numberOfTOPNodes, 1);

      int out = 0;
      for(const auto& node : *nodes){
          if(node->isStaticNode())
              out++;
      }
      //return number of static nodes
      return out;
    }

    int buildNeighbourBigMapMultipleLayers__TEST(){

      NodeVector ptr(new vector<unique_ptr<NodeInterface>>);
      nodes = std::move(ptr);

      DataReader reader;
      int numberOfTOPNodes = 22;
      nodes->reserve(numberOfTOPNodes * 3);

      reader.readLandMask(nodes,"tests/files/motherDefSIMPLE2.csv",numberOfTOPNodes,0);
      reader.buildBottomLayers(nodes,3);
      reader.buildNeighbourMap(nodes, numberOfTOPNodes, 3);

      int out = 0;
      for(const auto& node : *nodes){
          if(node->isStaticNode())
              out++;
      }
      //return number of static nodes
      return out;
    }

    void roundValue__TEST() {

    }

};

#endif //DATAREADERFIXTURE_HPP