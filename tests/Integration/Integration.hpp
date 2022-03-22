#include <gtest/gtest-death-test.h>
#include "fixtures/SimulationFixture.hpp"

TEST_F (SimulationFixture, GridReadTest) {
    ASSERT_EQ(Test__readNodesLayerFromFile(), 1200);
    ASSERT_EQ(Test__runEquationStep(),0);
    ASSERT_EQ(Test__runSimpleWithrivers(),0);
    ASSERT_EQ(Test__runSimpleWithrivers3D(),0);
}
