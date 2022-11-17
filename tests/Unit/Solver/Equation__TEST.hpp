#include <gtest/gtest.h>
#include "../Simulation/Options__MOCK.hpp"
#include "../../../src/Model/Node.hpp"
#include "../../../src/Solver/Equation.hpp"

using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

class EquationFixture : public ::testing::Test {
public:
    NodeVector nodes;
    Equation *eq;
    MockOptions options;
    DensityProperties densityProperties;

    void SetUp() {
        NodeVector ptr(new std::vector<unique_ptr<GlobalFlow::Model::NodeInterface>>);

        nodes = std::move(ptr);
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 0, 0, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, densityProperties
                ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 1, 1, 0.2 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, densityProperties
                ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 1, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 2, 2, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, densityProperties
                ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 1, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, densityProperties
                ));

        nodes->at(0)->setNeighbour(1, RIGHT);
        nodes->at(1)->setNeighbour(0, LEFT);

        nodes->at(0)->setNeighbour(2, DOWN);
        nodes->at(2)->setNeighbour(0, TOP);

        nodes->at(2)->setNeighbour(3, DOWN);
        nodes->at(3)->setNeighbour(2, TOP);

        eq = new Equation(4, nodes, options);
    }

    using p_node = std::unique_ptr<GlobalFlow::Model::NodeInterface>;

    p_node &at(int pos) { return nodes->at(pos); }
};

TEST_F(EquationFixture, toggleSteadyState) {
    ASSERT_TRUE(eq->toggleSteadyState());
    ASSERT_FALSE(eq->toggleSteadyState());
}

TEST_F(EquationFixture, updateStepSize) {
    ASSERT_EQ((at(0)->getProperties().get<t_dim, StepModifier>().value()), 1);
    eq->updateStepSize(10);
    ASSERT_EQ((at(0)->getProperties().get<t_dim, StepModifier>().value()), 10);
}

/**
 * Also tests getError() and getItter()
 * TODO this should fail not a real test
 */
TEST_F(EquationFixture, solve) {
    ON_CALL(options, getMaxIterations()).WillByDefault(Return(5));
    ON_CALL(options, getConverganceCriteria()).WillByDefault(Return(.1));
    ON_CALL(options, getInitialHead()).WillByDefault(Return(1));
    ON_CALL(options, getMaxHeadChange()).WillByDefault(Return(1));
    ON_CALL(options, isDampingEnabled()).WillByDefault(Return(false));
    eq = new Equation(4, nodes, options);
    eq->solve();
    // todo ASSERT...
}

TEST_F(EquationFixture, getResiduals) {
    eq->solve();
    std::cout << eq->getResiduals();
    // todo ASSERT...
}

/*
TEST_F(EquationFixture, updateClosingCrit) {
    FAIL();
}

TEST_F(EquationFixture, getResults) {
    FAIL();
}

TEST_F(EquationFixture, coutOperator) {
    FAIL();
}
*/
