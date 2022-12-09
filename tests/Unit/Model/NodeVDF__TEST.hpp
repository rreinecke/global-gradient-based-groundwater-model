#include <gtest/gtest.h>
#include "../../../src/Model/Node.hpp"

using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

class StandardNodeVDFFixture : public ::testing::Test {
public:
    NodeVector nodes;
    DensityProperties densityProperties;

    void SetUp() {
        densityProperties = GlobalFlow::Model::DensityProperties::setDensityProperties(true, {1000.0, 1012.5, 1025.0}, 0.2, 0.2);

        NodeVector ptr(new vector <unique_ptr<GlobalFlow::Model::NodeInterface>>);
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
                nodes, 0, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 2, 2, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, densityProperties
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, densityProperties
        ));

        at(0)->setElevation(10 * si::meter);
        at(1)->setElevation(10 * si::meter);

        // set nodes 0 and 1 as horizontal neighbours
        nodes->at(0)->setNeighbour(1, RIGHT);
        nodes->at(1)->setNeighbour(0, LEFT);

        // set nodes 0 and 2 as vertical neighbours (0 on top of 2)
        nodes->at(0)->setNeighbour(2, DOWN);
        nodes->at(2)->setNeighbour(0, TOP);

        // set nodes 1 and 3 as vertical neighbours (1 on top of 3)
        nodes->at(1)->setNeighbour(3, DOWN);
        nodes->at(3)->setNeighbour(1, TOP);

        // set nodes 2 and 3 as horizontal neighbours
        nodes->at(2)->setNeighbour(3, RIGHT);
        nodes->at(3)->setNeighbour(2, LEFT);

        // zeta surface 0 (all at top of nodes)
        nodes->at(0)->addInitialZeta(10.0 * si::meter);
        nodes->at(1)->addInitialZeta(10.0 * si::meter);
        nodes->at(2)->addInitialZeta(0.0 * si::meter);
        nodes->at(3)->addInitialZeta(0.0 * si::meter);

        // zeta surface 1 (in nodes 0 & 1 between, in nodes 2 & 3 at top)
        nodes->at(0)->addInitialZeta(7.5 * si::meter);
        nodes->at(1)->addInitialZeta(2.5 * si::meter);
        nodes->at(2)->addInitialZeta(0.0 * si::meter);
        nodes->at(3)->addInitialZeta(0.0 * si::meter);

        // zeta surface 2 (in nodes 0 & 1 at bottom, in nodes 2 & 3 between)
        nodes->at(0)->addInitialZeta(0.0 * si::meter);
        nodes->at(1)->addInitialZeta(0.0 * si::meter);
        nodes->at(2)->addInitialZeta(-2.5 * si::meter);
        nodes->at(3)->addInitialZeta(-7.5 * si::meter);

        // zeta surface 3 (all at bottom of nodes)
        nodes->at(0)->addInitialZeta(0.0 * si::meter);
        nodes->at(1)->addInitialZeta(0.0 * si::meter);
        nodes->at(2)->addInitialZeta(-10.0 * si::meter);
        nodes->at(3)->addInitialZeta(-10.0 * si::meter);
    }

    using p_node = std::unique_ptr<GlobalFlow::Model::NodeInterface>;

    p_node &at(int pos) { return nodes->at(pos); }
};

TEST_F(StandardNodeVDFFixture, setEffectivePorosity) {
    at(0)->setEffectivePorosity(0.2 * si::si_dimensionless);
    ASSERT_EQ((at(0)->getProperties().get<t_dim, EffectivePorosity>().value()), 0.2);
}

TEST_F(StandardNodeVDFFixture, setZeta) {
    at(0)->setZeta(0, 9 * si::meter);
    ASSERT_EQ(at(0)->getZeta(0).value(),9);
}

TEST_F(StandardNodeVDFFixture, setZetaChange) {
    at(0)->setZetaChange(0, -1 * si::meter);
    ASSERT_EQ(at(0)->getZetaChange(0).value(),-1);
}

TEST_F(StandardNodeVDFFixture, setTopZetaToHead) {
    // for confined node nothing should be done
    at(0)->setZeta(0, 20 * si::meter);
    at(0)->setTopZetaToHead();
    ASSERT_EQ(at(0)->getZeta(0).value(),20);

    // add and test unconfined node
    nodes->emplace_back(new GlobalFlow::Model::StandardNode(
            nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day, 1, 10, 1,
            0.2, 0.1, false, densityProperties
    ));
    at(4)->setElevation(10 * si::meter);
    at(4)->setHead_direct(10);
    // - zeta height ABOVE head (should be set to head)
    nodes->at(4)->addInitialZeta(20.0 * si::meter);
    at(4)->setTopZetaToHead();
    ASSERT_EQ(at(4)->getZeta(0).value(),10);
    // - zeta height BELOW head (should be set to head)
    at(4)->setZeta(0, 5 * si::meter);
    at(4)->setTopZetaToHead();
    ASSERT_EQ(at(4)->getZeta(0).value(),10);
}

TEST_F(StandardNodeVDFFixture, setZetaPosInNode) {
    at(0)->setZetaPosInNode(0);
    at(0)->setZetaPosInNode(1);
    at(0)->setZetaPosInNode(2);
    at(0)->setZetaPosInNode(3);
    ASSERT_EQ((at(0)->getZetaPosInNode(0)), "between");
    ASSERT_EQ((at(0)->getZetaPosInNode(1)), "between");
    ASSERT_EQ((at(0)->getZetaPosInNode(2)), "bottom");
    ASSERT_EQ((at(0)->getZetaPosInNode(3)), "bottom");
}

TEST_F(StandardNodeVDFFixture, getNusTop) {
    at(0)->setZetaPosInNode(0);
    at(0)->setZetaPosInNode(1);
    at(0)->setZetaPosInNode(2);
    at(0)->setZetaPosInNode(3);
    ASSERT_EQ((at(0)->getNusTop().value()), 0.0);
    at(2)->setZetaPosInNode(0);
    at(2)->setZetaPosInNode(1);
    at(2)->setZetaPosInNode(2);
    at(2)->setZetaPosInNode(3);
    ASSERT_EQ((at(2)->getNusTop().value()), 0.0125);
}

TEST_F(StandardNodeVDFFixture, getNusBot) {
    at(0)->setZetaPosInNode(0);
    at(0)->setZetaPosInNode(1);
    at(0)->setZetaPosInNode(2);
    at(0)->setZetaPosInNode(3);
    ASSERT_EQ((at(0)->getNusBot().value()), 0.0125);
    at(2)->setZetaPosInNode(0);
    at(2)->setZetaPosInNode(1);
    at(2)->setZetaPosInNode(2);
    at(2)->setZetaPosInNode(3);
    ASSERT_EQ((at(2)->getNusBot().value()), 0.025);
}

TEST_F(StandardNodeVDFFixture, getRHSConstantDensity) {
    at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getRHSConstantDensity().value()), -50);
}

TEST_F(StandardNodeVDFFixture, getZoneConductance) {
    unordered_map<NeighbourPosition, large_num> neighbourList = at(0)->getListOfNeighbours();
    std::unordered_map<NeighbourPosition, large_num>::const_iterator got = neighbourList.find(RIGHT);
    for (int nodeID = 0; nodeID <= 1; nodeID++){
        for(int zetaID = 0; zetaID <= 3; zetaID++) {
            at(nodeID)->setZetaPosInNode(zetaID);
        }
    }
    ASSERT_NEAR((at(0)->getZoneConductances(got)[0].value()), 0.666, 0.01);
    //ASSERT_NEAR((at(0)->getZoneConductances(got)[1].value()), 0.666, 0.01); // todo
    ASSERT_NEAR((at(0)->getZoneConductances(got)[2].value()), 0, 0.01);
}

TEST_F(StandardNodeVDFFixture, getZoneConductanceCum) {
    unordered_map<NeighbourPosition, large_num> neighbourList = at(0)->getListOfNeighbours();
    std::unordered_map<NeighbourPosition, large_num>::const_iterator got = neighbourList.find(RIGHT);
    for (int nodeID = 0; nodeID <= 1; nodeID++){
        for(int zetaID = 0; zetaID <= 3; zetaID++) {
            at(nodeID)->setZetaPosInNode(zetaID);
        }
    }
    // todo
}

TEST_F(StandardNodeVDFFixture, getRHS) {
    at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);

    //ASSERT_EQ((at(0)->getRHS().value()), 0); // todo calculate correct result
}

TEST_F(StandardNodeVDFFixture, getZetaRHS) {
    //at(0)->getZetaRHS(1);
    // todo
}

TEST_F(StandardNodeVDFFixture, getEffectivePorosityTerm) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getVDFMatrixEntries) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getSourceTermBelowZeta) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getTipToeFlow) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getFlowPseudoSource) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getVerticalFluxCorrection) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getFluxCorrTop) {
    // todo
}

TEST_F(StandardNodeVDFFixture, getFluxCorrDown) {
    // todo
}

TEST_F(StandardNodeVDFFixture, verticalZetaMovement) {
    // todo
}

TEST_F(StandardNodeVDFFixture, horizontalZetaMovement) {
    // todo
}

TEST_F(StandardNodeVDFFixture, clipInnerZetas) {
    for(int localZetaID = 0; localZetaID < 4; localZetaID++) { at(0)->setZetaPosInNode(localZetaID); }
    at(0)->setZeta(1, 20 * si::meter); // set zeta outside the upper bound (Zetas.front())
    at(0)->clipInnerZetas();
    ASSERT_EQ(at(0)->getZeta(1).value(), 10);

    for(int localZetaID = 0; localZetaID < 4; localZetaID++) { at(0)->setZetaPosInNode(localZetaID); }
    at(0)->setZeta(1, -20 * si::meter);
    at(0)->clipInnerZetas();
    ASSERT_EQ(at(0)->getZeta(3).value(), 0);
}

TEST_F(StandardNodeVDFFixture, correctCrossingZetas) {
    // todo
}

TEST_F(StandardNodeVDFFixture, preventZetaLocking) {
    // todo
}

TEST_F(StandardNodeVDFFixture, setZoneOfSinksAndSources) {
    // todo
}

