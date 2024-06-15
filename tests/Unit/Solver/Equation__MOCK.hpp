#ifndef TESTING_EQUATION_MOCK_HPP
#define TESTING_EQUATION_MOCK_HPP

#include <gmock/gmock.h>
#include "../../../src/Solver/Equation.hpp"
using namespace GlobalFlow::Solver;

class MockEquation : public Equation{
public:
    MockEquation(NodeVector nodes, GlobalFlow::Simulation::Options op) : Equation(nodes, op) {}
    MOCK_METHOD1(updateStepSize, void(double));
};

#endif //TESTING_EQUATION_MOCK_HPP
