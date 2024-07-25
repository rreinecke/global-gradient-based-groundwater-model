#include <gtest/gtest.h>
#include "../../../src/Model/Node.hpp"

using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

class StandardNodeFixture : public ::testing::Test {
public:
    NodeVector nodes;

    void SetUp() {
        NodeVector ptr(new std::vector <std::unique_ptr<GlobalFlow::Model::NodeInterface>>);
        nodes = std::move(ptr);
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 0, 0, 0.1 * si::meter / day,
                1 * si::meter, 10, 1, 0.2, 0.1, true, false, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 1, 1, 0.2 * si::meter / day,
                1 * si::meter, 10, 1, 0.2, 0.1, true, false, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 1, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 2, 2, 0.1 * si::meter / day,
                1 * si::meter, 10, 1, 0.2, 0.1, true, false, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 1, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day,
                1 * si::meter, 10, 1, 0.2, 0.1, true, false, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
        ));

        nodes->at(0)->setNeighbour(1, RIGHT);
        nodes->at(1)->setNeighbour(0, LEFT);

        nodes->at(0)->setNeighbour(2, DOWN);
        nodes->at(2)->setNeighbour(0, TOP);

        nodes->at(2)->setNeighbour(3, DOWN);
        nodes->at(3)->setNeighbour(2, TOP);
    };

    using p_node = std::unique_ptr<GlobalFlow::Model::NodeInterface>;

    p_node &at(int pos) { return nodes->at(pos); }
};

TEST_F(StandardNodeFixture, setElevation_allLayers) {
    at(0)->setElevation_allLayers(10 * si::meter);
    ASSERT_EQ((at(0)->getProperties().get<t_meter, Elevation>().value()), 10);
    ASSERT_EQ((at(2)->getProperties().get<t_meter, Elevation>().value()), 0);
    ASSERT_EQ((at(3)->getProperties().get<t_meter, Elevation>().value()), -10);
    ASSERT_EQ((at(3)->getProperties().get<t_meter, TopElevation>().value()), 10);
}

TEST_F(StandardNodeFixture, setEfold) {
    at(0)->setEfold(10);
    ASSERT_EQ((at(0)->getProperties().get<t_meter, EFolding>().value()), 10);
    ASSERT_EQ((at(2)->getProperties().get<t_meter, EFolding>().value()), 10);
}

TEST_F(StandardNodeFixture, setEqHead) {
    at(0)->setElevation_allLayers(10 * si::meter);
    at(0)->setEqHead_allLayers(5 * si::meter);
    ASSERT_EQ((at(0)->getProperties().get<t_meter, EQHead>().value()), 5);
    ASSERT_EQ((at(2)->getProperties().get<t_meter, Head>().value()), 5);
}

TEST_F(StandardNodeFixture, getEqFlow) {
    at(0)->setElevation_allLayers(10 * si::meter);
    at(0)->setEqHead_allLayers(5 * si::meter);
    at(1)->setElevation_allLayers(20 * si::meter);
    at(1)->setEqHead_allLayers(5 * si::meter);
    ASSERT_NEAR(at(0)->getEqFlow().value(), 13.33, 0.1);
    ASSERT_NEAR(at(1)->getEqFlow().value(), -13.33, 0.1);
}

TEST_F(StandardNodeFixture, getK) {
    ASSERT_NEAR(at(0)->getK().value(), 0.1, 0.001);
    at(0)->getProperties().set<t_meter, EFolding>(0.1 * si::meter);
    ASSERT_NEAR(at(0)->getK().value(), 0.1, 0.0001);
    at(0)->getProperties().set<t_dim, StepSize>(0.1);
    ASSERT_NEAR(at(0)->getK().value(), 0.01, 0.001);
}

TEST_F(StandardNodeFixture, setK) {
    at(0)->setK(0.2 * si::meter / day);
    ASSERT_NEAR((at(0)->getK().value()), 0.2, 1e-7);
    ASSERT_NEAR((at(2)->getK().value()), 0.2, 1e-7);
}

TEST_F(StandardNodeFixture, getStorageCapacity) {
    at(0)->setElevation_allLayers(10 * si::meter);
    at(0)->setHead(5 * si::meter);
    ASSERT_EQ(at(0)->getStorageCapacity().value(), 0.2);
    at(0)->setElevation_allLayers(5 * si::meter);
    ASSERT_EQ(at(0)->getStorageCapacity().value(), 0.2);
}

TEST_F(StandardNodeFixture, AddAndRemoveExternalFlow) {
    int id = at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    ASSERT_EQ(id, 1);
    at(0)->removeExternalFlow(RECHARGE);
    id = at(0)->addExternalFlow(RIVER, 1 * si::meter, 0.1, 1 * si::meter);
    ASSERT_EQ(id, 1);
}

TEST_F(StandardNodeFixture, getExternalFlowByName) {
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    //at(0)->addExternalFlow(RECHARGE,0,5,0);
    ASSERT_EQ((at(0)->getExternalFlowByName(RECHARGE).getRecharge().value()), 50);
    ASSERT_EQ((at(0)->getExternalFlowByName(RECHARGE).getType()), RECHARGE);
    at(0)->removeExternalFlow(RECHARGE);
    //Currently no support for mutliple recharge types
    //ASSERT_EQ((at(0)->getExternalFlowByName(RECHARGE).getRecharge().value()), 50);
    //at(0)->removeExternalFlow(RECHARGE);
    try {
        at(0)->getExternalFlowByName(RECHARGE);
        FAIL();
    } catch (std::out_of_range const &err) { SUCCEED(); }
}

TEST_F(StandardNodeFixture, calculateExternalFlowVolume) {
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    ASSERT_EQ((at(0)->calculateExternalFlowVolume(at(0)->getExternalFlowByName(RECHARGE)).value()), 50);
}

TEST_F(StandardNodeFixture, getExternalFlowVolumeByName) {
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    ASSERT_EQ((at(0)->getExternalFlowVolumeByName(RECHARGE).value()), 50);
}

TEST_F(StandardNodeFixture, getNeighbour) {
    ASSERT_EQ((at(0)->getNeighbour(DOWN)->getID()), 2);
    ASSERT_EQ((at(2)->getNeighbour(TOP)->getID()), 0);
    ASSERT_EQ((at(1)->getNeighbour(LEFT)->getID()), 0);
}

TEST_F(StandardNodeFixture, calculateDewateredFlow) {
    at(2)->setHead(1 * si::meter);
    at(0)->setElevation_allLayers(20 * si::meter);
    at(0)->setHead(5 * si::meter);
    ASSERT_EQ((at(2)->calculateDewateredFlow().value()), 0);
    ASSERT_EQ((at(0)->calculateDewateredFlow().value()), 0);
}

TEST_F(StandardNodeFixture, updateUniqueFlow) {
    at(0)->updateUniqueFlow(5);
    ASSERT_EQ((at(0)->getExternalFlowVolumeByName(RECHARGE).value()), 5);
    at(0)->removeExternalFlow(RECHARGE);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->updateUniqueFlow(5);
    ASSERT_EQ((at(0)->getExternalFlowVolumeByName(RECHARGE).value()), 5);
}

TEST_F(StandardNodeFixture, getQ) {
    at(0)->setHead(1 * si::meter);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getQ().value()), 150);
    at(0)->setHead(-50 * si::meter);
    ASSERT_NE((at(0)->getQ().value()), 1);
}

TEST_F(StandardNodeFixture, getP) {
    at(0)->setHead(1 * si::meter);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getP_aboveFlowBottom().value()), 0);
    at(0)->setHead(-50 * si::meter);
    ASSERT_EQ((at(0)->getP_aboveFlowBottom().value()), 0);
}

TEST_F(StandardNodeFixture, calculateNotHeadDependandFlows) {
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->setHead(100 * si::meter);
    ASSERT_EQ(at(0)->getP_belowFlowBottom().value(), 0);
    at(0)->setHead(0 * si::meter);
    ASSERT_EQ(at(0)->getP_belowFlowBottom().value(), -50);
}

TEST_F(StandardNodeFixture, getConductance) {
    ASSERT_NEAR((at(0)->getMatrixEntries()[0].value()), -1.54, 0.01);
    ASSERT_NEAR((at(0)->getMatrixEntries()[1].value()), 1.33, 0.01);
}

TEST_F(StandardNodeFixture, getRHSConstantDensity) {
    at(0)->setHead(1 * si::meter);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getRHS().value()), -50);
}

TEST_F(StandardNodeFixture, getRHS) {
    at(0)->setHead(1 * si::meter);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getRHS().value()), -50);
}



