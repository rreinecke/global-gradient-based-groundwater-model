#ifndef TESTING_OPTIONS_MOCK_HPP
#define TESTING_OPTIONS_MOCK_HPP

#include <gmock/gmock.h>
#include "../../../src/Simulation/Options.hpp"

using ::testing::Return;

class MockOptions : public Options {
public:
    MOCK_METHOD0(getMaxIterations, int());

    MOCK_METHOD0(getConverganceCriteria, double());

    MOCK_METHOD0(getInitialHead, double());

    MOCK_METHOD0(getMaxHeadChange, double());

    MOCK_METHOD0(isDampingEnabled, bool());
};

#endif