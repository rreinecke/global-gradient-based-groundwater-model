#include <gtest/gtest.h>
#include "../../fixtures/DataReaderFixture.hpp"

// NEW Neighbouring
TEST_F (DataReaderFixture, SmallNeighbouringTests) {
buildNeighbourMap_NEW__TEST();
// center node with 4 existing neighbours in landmask
//ASSERT_EQ(nodes->at(6)->getNeighbour(Model::EAST)->getSpatialID(), 7); // alternatively internal id as value to compare
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::EAST)->getSpatialID(), 39);
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::WEST)->getSpatialID(), 37);
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::NORTH)->getSpatialID(), 2);
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::SOUTH)->getSpatialID(), 74);

// center node with 1 existing neighbour in landmask
ASSERT_EQ(nodes->at(15)->getNeighbour(Model::EAST)->getSpatialID(), 591);
ASSERT_THROW(nodes->at(15)->getNeighbour(Model::WEST)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
ASSERT_THROW(nodes->at(15)->getNeighbour(Model::NORTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
ASSERT_THROW(nodes->at(15)->getNeighbour(Model::SOUTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// ASSERT_TRUE(nodes->at(15)->hasTypeOfExternalFlow(Model::GENERAL_HEAD_BOUNDARY)); // Check for Ocean Boundary, only existence of minimum of 1 border, not count of borders
}

TEST_F (DataReaderFixture, SmallNeighbouringTestsCorners) {
buildNeighbourMap_NEW__TEST();
// 4 grid corners:
// upper left (North-West) corner of global grid
ASSERT_EQ(nodes->at(0)->getNeighbour(Model::EAST)->getSpatialID(), 2);
ASSERT_EQ(nodes->at(0)->getNeighbour(Model::WEST)->getSpatialID(), 36);
ASSERT_EQ(nodes->at(0)->getNeighbour(Model::SOUTH)->getSpatialID(), 37);
ASSERT_THROW(nodes->at(0)->getNeighbour(Model::NORTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// upper right (North-East) corner of global grid
ASSERT_EQ(nodes->at(4)->getNeighbour(Model::EAST)->getSpatialID(), 1);
ASSERT_EQ(nodes->at(4)->getNeighbour(Model::WEST)->getSpatialID(), 35);
ASSERT_EQ(nodes->at(4)->getNeighbour(Model::SOUTH)->getSpatialID(), 72);
ASSERT_THROW(nodes->at(4)->getNeighbour(Model::NORTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// lower left (South-West) corner of global grid
ASSERT_EQ(nodes->at(18)->getNeighbour(Model::EAST)->getSpatialID(), 614);
ASSERT_EQ(nodes->at(18)->getNeighbour(Model::WEST)->getSpatialID(), 648);
ASSERT_EQ(nodes->at(18)->getNeighbour(Model::NORTH)->getSpatialID(), 577);
ASSERT_THROW(nodes->at(18)->getNeighbour(Model::SOUTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// lower right (South-East) corner of global grid
ASSERT_EQ(nodes->at(22)->getNeighbour(Model::EAST)->getSpatialID(), 613);
ASSERT_EQ(nodes->at(22)->getNeighbour(Model::WEST)->getSpatialID(), 647);
ASSERT_EQ(nodes->at(22)->getNeighbour(Model::NORTH)->getSpatialID(), 612);
ASSERT_THROW(nodes->at(22)->getNeighbour(Model::SOUTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);

}

TEST_F (DataReaderFixture, SmallNeighbouringTestsFirstLast) {
buildNeighbourMap_NEW__TEST();
// first and last row and column
// first column (west)
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::EAST)->getSpatialID(), 38);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::WEST)->getSpatialID(), 72);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::NORTH)->getSpatialID(), 1);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::SOUTH)->getSpatialID(), 73);
// last column (east)
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::EAST)->getSpatialID(), 37);
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::WEST)->getSpatialID(), 71);
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::NORTH)->getSpatialID(), 36);
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::SOUTH)->getSpatialID(), 108);
// first row (north)
ASSERT_EQ(nodes->at(1)->getNeighbour(Model::EAST)->getSpatialID(), 3);
ASSERT_EQ(nodes->at(1)->getNeighbour(Model::WEST)->getSpatialID(), 1);
ASSERT_EQ(nodes->at(1)->getNeighbour(Model::SOUTH)->getSpatialID(), 38);
ASSERT_THROW(nodes->at(1)->getNeighbour(Model::NORTH)->getSpatialID(), GlobalFlow::Model::StandardNode::NodeNotFoundException);
// last row (south)
ASSERT_EQ(nodes->at(19)->getNeighbour(Model::EAST)->getSpatialID(), 615);
ASSERT_EQ(nodes->at(19)->getNeighbour(Model::WEST)->getSpatialID(), 613);
ASSERT_EQ(nodes->at(19)->getNeighbour(Model::NORTH)->getSpatialID(), 578);
ASSERT_THROW(nodes->at(19)->getNeighbour(Model::SOUTH)->getSpatialID(), GlobalFlow::Model::StandardNode::NodeNotFoundException);
}

TEST_F (DataReaderFixture, MultiLayerNeighbouringTests) {
buildNeighbourBigMapMultipleLayers_NEW__TEST();
// center node with 4 existing neighbours in landmask
//ASSERT_EQ(nodes->at(6)->getNeighbour(Model::EAST)->getSpatialID(), 7); // alternatively internal id as value to compare
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::EAST)->getSpatialID(), 39);
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::WEST)->getSpatialID(), 37);
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::NORTH)->getSpatialID(), 2);
ASSERT_EQ(nodes->at(6)->getNeighbour(Model::SOUTH)->getSpatialID(), 74);

// center node with 1 existing neighbour in landmask
ASSERT_EQ(nodes->at(15)->getNeighbour(Model::EAST)->getSpatialID(), 591);
ASSERT_THROW(nodes->at(15)->getNeighbour(Model::WEST)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
ASSERT_THROW(nodes->at(15)->getNeighbour(Model::NORTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
ASSERT_THROW(nodes->at(15)->getNeighbour(Model::SOUTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// ASSERT_TRUE(nodes->at(15)->hasTypeOfExternalFlow(Model::GENERAL_HEAD_BOUNDARY)); // Check for Ocean Boundary, only existence of minimum of 1 border, not count of borders

// Test for nonTop Layers
// layer with in between others
// center node with 4 existing neighbours in landmask
ASSERT_EQ(nodes->at(29)->getSpatialID(), 38);

ASSERT_EQ(nodes->at(29)->getNeighbour(Model::TOP)->getSpatialID(), 38);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::TOP)->getID(), 6);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::DOWN)->getSpatialID(), 38);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::DOWN)->getID(), 52);

ASSERT_EQ(nodes->at(29)->getNeighbour(Model::EAST)->getID(), 30);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::WEST)->getID(), 28);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::NORTH)->getID(), 24);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::SOUTH)->getID(), 34);
}


TEST_F (DataReaderFixture, MultiLayerNeighbouringTestsCorners) {
buildNeighbourBigMapMultipleLayers_NEW__TEST();
// 4 grid corners:
// upper left (North-West) corner of global grid
ASSERT_EQ(nodes->at(0)->getNeighbour(Model::EAST)->getSpatialID(), 2);
ASSERT_EQ(nodes->at(0)->getNeighbour(Model::WEST)->getSpatialID(), 36);
ASSERT_EQ(nodes->at(0)->getNeighbour(Model::SOUTH)->getSpatialID(), 37);
ASSERT_THROW(nodes->at(0)->getNeighbour(Model::NORTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// upper right (North-East) corner of global grid
ASSERT_EQ(nodes->at(4)->getNeighbour(Model::EAST)->getSpatialID(), 1);
ASSERT_EQ(nodes->at(4)->getNeighbour(Model::WEST)->getSpatialID(), 35);
ASSERT_EQ(nodes->at(4)->getNeighbour(Model::SOUTH)->getSpatialID(), 72);
ASSERT_THROW(nodes->at(4)->getNeighbour(Model::NORTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// lower left (South-West) corner of global grid
ASSERT_EQ(nodes->at(18)->getNeighbour(Model::EAST)->getSpatialID(), 614);
ASSERT_EQ(nodes->at(18)->getNeighbour(Model::WEST)->getSpatialID(), 648);
ASSERT_EQ(nodes->at(18)->getNeighbour(Model::NORTH)->getSpatialID(), 577);
ASSERT_THROW(nodes->at(18)->getNeighbour(Model::SOUTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
// lower right (South-East) corner of global grid
ASSERT_EQ(nodes->at(22)->getNeighbour(Model::EAST)->getSpatialID(), 613);
ASSERT_EQ(nodes->at(22)->getNeighbour(Model::WEST)->getSpatialID(), 647);
ASSERT_EQ(nodes->at(22)->getNeighbour(Model::NORTH)->getSpatialID(), 612);
ASSERT_THROW(nodes->at(22)->getNeighbour(Model::SOUTH)->getSpatialID(),
        GlobalFlow::Model::StandardNode::NodeNotFoundException);
}

TEST_F (DataReaderFixture, MultiLayerNeighbouringTestsFirstLast) {
buildNeighbourBigMapMultipleLayers_NEW__TEST();
// first and last row and column
// first column (west)
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::EAST)->getSpatialID(), 38);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::WEST)->getSpatialID(), 72);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::NORTH)->getSpatialID(), 1);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::SOUTH)->getSpatialID(), 73);
// last column (east)
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::EAST)->getSpatialID(), 37);
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::WEST)->getSpatialID(), 71);
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::NORTH)->getSpatialID(), 36);
ASSERT_EQ(nodes->at(9)->getNeighbour(Model::SOUTH)->getSpatialID(), 108);
// first row (north)
ASSERT_EQ(nodes->at(1)->getNeighbour(Model::EAST)->getSpatialID(), 3);
ASSERT_EQ(nodes->at(1)->getNeighbour(Model::WEST)->getSpatialID(), 1);
ASSERT_EQ(nodes->at(1)->getNeighbour(Model::SOUTH)->getSpatialID(), 38);
ASSERT_THROW(nodes->at(1)->getNeighbour(Model::NORTH)->getSpatialID(), GlobalFlow::Model::StandardNode::NodeNotFoundException);
// last row (south)
ASSERT_EQ(nodes->at(19)->getNeighbour(Model::EAST)->getSpatialID(), 615);
ASSERT_EQ(nodes->at(19)->getNeighbour(Model::WEST)->getSpatialID(), 613);
ASSERT_EQ(nodes->at(19)->getNeighbour(Model::NORTH)->getSpatialID(), 578);
ASSERT_THROW(nodes->at(19)->getNeighbour(Model::SOUTH)->getSpatialID(), GlobalFlow::Model::StandardNode::NodeNotFoundException);
}



// OLD Neighbouring
TEST_F (DataReaderFixture, SmallNeighbouringTests_OLD) {
buildNeighbourMap__TEST();
ASSERT_EQ(nodes->at(2)->getNeighbour(Model::EAST)->getSpatialID(), 4);
ASSERT_EQ(nodes->at(2)->getNeighbour(Model::WEST)->getSpatialID(), 2);
ASSERT_EQ(nodes->at(11)->getNeighbour(Model::WEST)->getSpatialID(), 11);
ASSERT_EQ(nodes->at(7)->getNeighbour(Model::NORTH)->getSpatialID(), 6);
ASSERT_EQ(nodes->at(7)->getNeighbour(Model::SOUTH)->getSpatialID(), 12);
}

TEST_F (DataReaderFixture, BigNeighbouringTests_OLD) {
buildNeighbourBigMap__TEST();
ASSERT_EQ(nodes->at(7)->getNeighbour(Model::SOUTH)->getSpatialID(), 13);
ASSERT_EQ(nodes->at(7)->getNeighbour(Model::NORTH)->getSpatialID(), 3);
ASSERT_EQ(nodes->at(8)->getNeighbour(Model::EAST)->getSpatialID(), 10);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::WEST)->getSpatialID(), 5);
ASSERT_EQ(nodes->at(20)->getNeighbour(Model::NORTH)->getSpatialID(), 15);
ASSERT_EQ(nodes->at(21)->getNeighbour(Model::WEST)->getSpatialID(), 21);
}

TEST_F (DataReaderFixture, MultiLayerNeighbouringTests_OLD) {
buildNeighbourBigMapMultipleLayers__TEST();
ASSERT_EQ(nodes->at(7)->getNeighbour(Model::SOUTH)->getSpatialID(), 13);
ASSERT_EQ(nodes->at(7)->getNeighbour(Model::NORTH)->getSpatialID(), 3);
ASSERT_EQ(nodes->at(8)->getNeighbour(Model::EAST)->getSpatialID(), 10);
ASSERT_EQ(nodes->at(5)->getNeighbour(Model::WEST)->getSpatialID(), 5);
ASSERT_EQ(nodes->at(20)->getNeighbour(Model::NORTH)->getSpatialID(), 15);
ASSERT_EQ(nodes->at(21)->getNeighbour(Model::WEST)->getSpatialID(), 21);

ASSERT_EQ(nodes->at(29)->getSpatialID(), 8);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::NORTH)->getSpatialID(), 3);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::SOUTH)->getSpatialID(), 13);

ASSERT_EQ(nodes->at(29)->getNeighbour(Model::TOP)->getSpatialID(), 8);
ASSERT_EQ(nodes->at(29)->getNeighbour(Model::DOWN)->getID(), 51);
}
