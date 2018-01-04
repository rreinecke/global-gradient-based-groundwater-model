#include <gtest/gtest.h>
#include "../../fixtures/DataReaderFixture.hpp"

TEST_F (DataReaderFixture, SmallNeighbouringTests) {
    ASSERT_EQ(buildNeighbourMap__TEST(),24);
    ASSERT_EQ(nodes->at(2)->getNeighbour(RIGHT)->getARCID(), 4);
    ASSERT_EQ(nodes->at(2)->getNeighbour(LEFT)->getARCID(), 2);
    ASSERT_EQ(nodes->at(11)->getNeighbour(LEFT)->getARCID(), 11);
    ASSERT_EQ(nodes->at(7)->getNeighbour(FRONT)->getARCID(), 6);
    ASSERT_EQ(nodes->at(7)->getNeighbour(BACK)->getARCID(), 12);
}

TEST_F (DataReaderFixture, BigNeighbouringTests) {
    ASSERT_EQ(buildNeighbourBigMap__TEST(),34);
    ASSERT_EQ(nodes->at(7)->getNeighbour(BACK)->getARCID(), 13);
    ASSERT_EQ(nodes->at(7)->getNeighbour(FRONT)->getARCID(), 3);
    ASSERT_EQ(nodes->at(8)->getNeighbour(RIGHT)->getARCID(), 10);
    ASSERT_EQ(nodes->at(5)->getNeighbour(LEFT)->getARCID(), 5);
    ASSERT_EQ(nodes->at(20)->getNeighbour(FRONT)->getARCID(), 15);
    ASSERT_EQ(nodes->at(21)->getNeighbour(LEFT)->getARCID(), 21);
}

TEST_F (DataReaderFixture, MultiLayerNeighbouringTests) {
    ASSERT_EQ(buildNeighbourBigMapMultipleLayers__TEST(),102);
    ASSERT_EQ(nodes->at(7)->getNeighbour(BACK)->getARCID(), 13);
    ASSERT_EQ(nodes->at(7)->getNeighbour(FRONT)->getARCID(), 3);
    ASSERT_EQ(nodes->at(8)->getNeighbour(RIGHT)->getARCID(), 10);
    ASSERT_EQ(nodes->at(5)->getNeighbour(LEFT)->getARCID(), 5);
    ASSERT_EQ(nodes->at(20)->getNeighbour(FRONT)->getARCID(), 15);
    ASSERT_EQ(nodes->at(21)->getNeighbour(LEFT)->getARCID(), 21);

    ASSERT_EQ(nodes->at(29)->getARCID(), 8);
    ASSERT_EQ(nodes->at(29)->getNeighbour(FRONT)->getARCID(), 3);
    ASSERT_EQ(nodes->at(29)->getNeighbour(BACK)->getID(), 34);

    ASSERT_EQ(nodes->at(29)->getNeighbour(TOP)->getARCID(), 8);
    ASSERT_EQ(nodes->at(29)->getNeighbour(DOWN)->getID(), 51);
}
