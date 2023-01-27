#include <gtest/gtest.h>
#include "../../../src/Model/Node.hpp"

using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

class StandardNodeVDFFixture : public ::testing::Test {
public:
    NodeVector nodes;

    void SetUp() {

        NodeVector ptr(new vector <unique_ptr<GlobalFlow::Model::NodeInterface>>);
        nodes = std::move(ptr);
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 0, 0, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, true
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 1, 1, 0.2 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, true
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 2, 2, 0.05 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, true
        ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day, 1, 10, 1,
                0.2, 0.1, true, true
        ));

        double densityFresh = 1000.0;
        vector<quantity<Dimensionless>> densityZones = {1000.0, 1010.0, 1020.0,1030.0};
        vector<quantity<Dimensionless>> nusInZones;
        vector<quantity<Dimensionless>> delnus;

        for (int id = 0; id < densityZones.size(); id++) {
            // nus of zones is equal to nus of zeta surface below
            nusInZones.push_back(((densityZones[id] - densityFresh) / densityFresh) * si::si_dimensionless);
            if (id == 0) {
                delnus.push_back(nusInZones[id]); // density difference in top zone
            } else {
                delnus.push_back((nusInZones[id] - nusInZones[id - 1]));
            }
        }

        large_num numberOfNodes = 4;
        for (large_num k = 0; k < numberOfNodes; ++k) {
            nodes->at(k)->setDelnus(delnus);
            nodes->at(k)->setNusInZones(nusInZones);
            nodes->at(k)->setMaxTipToeSlope(0.2 * si::si_dimensionless);
            nodes->at(k)->setSlopeAdjFactor(0.1 * si::si_dimensionless);
            nodes->at(k)->setMinDepthFactor(0.1 * si::si_dimensionless);
            nodes->at(k)->setVDFLock(0.001 * si::meter);
            nodes->at(k)->setEffectivePorosity(0.2 * si::si_dimensionless);
            nodes->at(k)->setZoneOfSinksAndSources(0, 3, 4);
        }

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

        at(0)->setElevation(10 * si::meter);
        at(1)->setElevation(10 * si::meter);
        at(0)->setHead_direct(10); // Question: what does that cause?
        at(1)->setHead_direct(10);
        at(2)->setHead_direct(0);
        at(3)->setHead_direct(0);

        // zeta surface 0 (all at top of nodes)
        nodes->at(0)->addInitialZeta(10.0 * si::meter);
        nodes->at(1)->addInitialZeta(10.0 * si::meter);
        nodes->at(2)->addInitialZeta(0.0 * si::meter);
        nodes->at(3)->addInitialZeta(0.0 * si::meter);

        // zeta surface 1 (in nodes 0 & 1 between, in nodes 2 & 3 at top)
        nodes->at(0)->addInitialZeta(9 * si::meter);
        nodes->at(1)->addInitialZeta(8 * si::meter);
        nodes->at(2)->addInitialZeta(0.0 * si::meter);
        nodes->at(3)->addInitialZeta(0.0 * si::meter);

        // zeta surface 2 (in nodes 0 & 1 between, in nodes 2 & 3 at top)
        nodes->at(0)->addInitialZeta(7.5 * si::meter);
        nodes->at(1)->addInitialZeta(2.5 * si::meter);
        nodes->at(2)->addInitialZeta(0.0 * si::meter);
        nodes->at(3)->addInitialZeta(0.0 * si::meter);

        // zeta surface 3 (in nodes 0 & 1 at bottom, in nodes 2 & 3 between)
        nodes->at(0)->addInitialZeta(0.0 * si::meter);
        nodes->at(1)->addInitialZeta(0.0 * si::meter);
        nodes->at(2)->addInitialZeta(-2.5 * si::meter);
        nodes->at(3)->addInitialZeta(-7.5 * si::meter);

        // zeta surface 4 (all at bottom of nodes)
        nodes->at(0)->addInitialZeta(0.0 * si::meter);
        nodes->at(1)->addInitialZeta(0.0 * si::meter);
        nodes->at(2)->addInitialZeta(-10.0 * si::meter);
        nodes->at(3)->addInitialZeta(-10.0 * si::meter);

        for (int nodeID = 0; nodeID <= 3; nodeID++){
            for(int zetaID = 0; zetaID <= 4; zetaID++) {
                at(nodeID)->setZetaPosInNode(zetaID);
            }
        }
    }

    using p_node = std::unique_ptr<GlobalFlow::Model::NodeInterface>;

    p_node &at(int pos) { return nodes->at(pos); }
};

TEST_F(StandardNodeVDFFixture, getEffectivePorosity) {
    ASSERT_EQ((at(0)->getProperties().get<t_dim, EffectivePorosity>().value()), 0.2);
}

TEST_F(StandardNodeVDFFixture, getDelnus) {
    auto delnus_node_0 = at(0)->getProperties().get<vector<t_dim>, Delnus>();
    ASSERT_EQ(delnus_node_0[0].value(),0.0);
    ASSERT_EQ(delnus_node_0[1].value(),0.01);
    ASSERT_EQ(delnus_node_0[2].value(),0.01);
}

TEST_F(StandardNodeVDFFixture, setZeta) {
    at(0)->setZeta(0, 9 * si::meter);
    ASSERT_EQ(at(0)->getZeta(0).value(),9);
}

TEST_F(StandardNodeVDFFixture, setTopZetaToHead) {
    // for confined node nothing should be done
    at(0)->setZeta(0, 20 * si::meter);
    at(0)->setTopZetaToHead();
    ASSERT_EQ(at(0)->getZeta(0).value(),20);

    // add and test unconfined node
    nodes->emplace_back(new GlobalFlow::Model::StandardNode(
            nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day, 1, 10, 1,
            0.2, 0.1, false, true
    ));
    at(4)->setElevation(10 * si::meter);
    at(4)->setHead_direct(10);

    // - zeta height ABOVE head (should be reset to head)
    nodes->at(4)->addInitialZeta(20.0 * si::meter);
    at(4)->setTopZetaToHead();
    ASSERT_EQ(at(4)->getZeta(0).value(),10);

    // - zeta height BELOW head (should be set to head)
    at(4)->setZeta(0, 5 * si::meter);
    at(4)->setTopZetaToHead();
    ASSERT_EQ(at(4)->getZeta(0).value(),10);
}

TEST_F(StandardNodeVDFFixture, setZetaChange) {
    at(0)->setZetaChange(0, 9 * si::meter);
    ASSERT_EQ(at(0)->getZetaChange(0).value(),-1);
}

TEST_F(StandardNodeVDFFixture, getZetaPosInNode) {
    ASSERT_EQ((at(0)->getZetaPosInNode(0)), "top");
    ASSERT_EQ((at(0)->getZetaPosInNode(1)), "between");
    ASSERT_EQ((at(0)->getZetaPosInNode(2)), "between");
    ASSERT_EQ((at(0)->getZetaPosInNode(3)), "bottom");
    ASSERT_EQ((at(0)->getZetaPosInNode(4)), "bottom");
}

TEST_F(StandardNodeVDFFixture, getNusTop) {
    ASSERT_EQ((at(0)->getNusTop().value()), 0.0);
    ASSERT_EQ((at(2)->getNusTop().value()), 0.02);
}

TEST_F(StandardNodeVDFFixture, getNusBot) {
    ASSERT_EQ((at(0)->getNusBot().value()), 0.02);
    ASSERT_EQ((at(2)->getNusBot().value()), 0.03);
}

TEST_F(StandardNodeVDFFixture, getZoneConductances) {
    std::unordered_map<NeighbourPosition, large_num>::const_iterator got;
    unordered_map<NeighbourPosition, large_num> neighbourList;

    neighbourList = at(0)->getListOfNeighbours();
    got = neighbourList.find(RIGHT);
    ASSERT_NEAR((at(0)->getZoneConductances(got)[0].value()), 0.199, 0.01);
    ASSERT_NEAR((at(0)->getZoneConductances(got)[1].value()), 0.466, 0.01);
    ASSERT_NEAR((at(0)->getZoneConductances(got)[2].value()), 0.0, 0.01); // zone is across layers -> 0
    ASSERT_NEAR((at(0)->getZoneConductances(got)[3].value()), 0, 0.01); // zone height in node = 0

    neighbourList = at(2)->getListOfNeighbours();
    got = neighbourList.find(RIGHT);
    ASSERT_NEAR((at(2)->getZoneConductances(got)[0].value()), 0, 0.01); // zone height in node = 0
    ASSERT_NEAR((at(2)->getZoneConductances(got)[1].value()), 0, 0.01); // zone height in node = 0
    ASSERT_NEAR((at(2)->getZoneConductances(got)[2].value()), 0, 0.01); // zone is across layers -> 0
    ASSERT_NEAR((at(2)->getZoneConductances(got)[3].value()), 0.333, 0.01);
}

TEST_F(StandardNodeVDFFixture, getZoneConductanceCum) {
    std::unordered_map<NeighbourPosition, large_num>::const_iterator got;
    unordered_map<NeighbourPosition, large_num> neighbourList;

    neighbourList = at(0)->getListOfNeighbours();
    got = neighbourList.find(RIGHT);
    ASSERT_NEAR((at(0)->getZoneConductanceCum(0, at(0)->getZoneConductances(got)).value()), 0.666, 0.01);
    ASSERT_NEAR((at(0)->getZoneConductanceCum(1, at(0)->getZoneConductances(got)).value()), 0.466, 0.01);
    ASSERT_NEAR((at(0)->getZoneConductanceCum(2, at(0)->getZoneConductances(got)).value()), 0.0, 0.01);
    neighbourList = at(2)->getListOfNeighbours();
    got = neighbourList.find(RIGHT);
    ASSERT_NEAR((at(2)->getZoneConductanceCum(0, at(2)->getZoneConductances(got)).value()), 0.333, 0.01);
    ASSERT_NEAR((at(2)->getZoneConductanceCum(2, at(2)->getZoneConductances(got)).value()), 0.333, 0.01);
}

TEST_F(StandardNodeVDFFixture, getRHSConstantDensity) {
    //at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getRHSConstantDensity().value()), -150); // with head set to 1, this is -50
}

TEST_F(StandardNodeVDFFixture, getPseudoSourceNode) {
    ASSERT_NEAR((at(0)->getPseudoSourceNode().value()), 0.00466666,0.0000001); // in zone 0: delnus=0, in zone 2,3: zoneConductance=0
    ASSERT_NEAR((at(1)->getPseudoSourceNode().value()), -0.00466666,0.0000001); // in zone 0: delnus=0, in zone 2,3: zoneConductance=0

    ASSERT_NEAR((at(2)->getPseudoSourceNode().value()), 0.01666666, 0.0000001); // zone 0,1,2: same zeta (at top)
    ASSERT_NEAR((at(3)->getPseudoSourceNode().value()), -0.01666666, 0.0000001); // zone 0,1,2: same zeta (at top)
}

TEST_F(StandardNodeVDFFixture, getVerticalFluxCorrection) {
    ASSERT_NEAR((at(2)->getVerticalFluxCorrection().value()), 0.00110000, 0.0000001);
    ASSERT_NEAR((at(3)->getVerticalFluxCorrection().value()), 0.00140000, 0.0000001);
    // todo test unconfined node
}

TEST_F(StandardNodeVDFFixture, getFluxCorrTop) {
    ASSERT_EQ((at(0)->getFluxCorrTop().value()), 0.0);
    ASSERT_EQ((at(1)->getFluxCorrTop().value()), 0.0);
    ASSERT_NEAR((at(2)->getFluxCorrTop().value()), -0.06776666, 0.0000001);
    ASSERT_NEAR((at(3)->getFluxCorrTop().value()), -0.13473333, 0.0000001);
    // todo test for unconfined node, are there SWI2 changes for confined nodes?
}

TEST_F(StandardNodeVDFFixture, getFluxCorrDown) {
    ASSERT_NEAR((at(0)->getFluxCorrDown().value()), 0.06776666, 0.0000001);
    ASSERT_NEAR((at(1)->getFluxCorrDown().value()), 0.13473333, 0.0000001);
    ASSERT_EQ((at(2)->getFluxCorrDown().value()), 0.0);
    ASSERT_EQ((at(3)->getFluxCorrDown().value()), 0.0);
    // todo test for unconfined node, are there SWI2 changes for confined nodes?
}

TEST_F(StandardNodeVDFFixture, getRHS) {
    //at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_NEAR((at(0)->getRHS().value()), -149.9275666, 0.0000001);
    // with getRHSConstantDensity = -150; pseudoSourceNode = 0.00466666;
    // fluxCorrectionTop = 0; fluxCorrectionDown = 0.06776666

    at(1)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(1)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(1)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_NEAR((at(1)->getRHS().value()), -149.8699333, 0.0000001);
    // with getRHSConstantDensity = -150; pseudoSourceNode = 0.00466666;
    // fluxCorrectionTop = 0; fluxCorrectionDown = 0.134733333
}

TEST_F(StandardNodeVDFFixture, getEffectivePorosityTerm) {
    ASSERT_EQ((at(0)->getEffectivePorosityTerm().value()), 0.2);
    ASSERT_EQ((at(0)->getEffectivePorosityTerm().value()), 0.2);
    ASSERT_EQ((at(0)->getEffectivePorosityTerm().value()), 0.2);
    ASSERT_EQ((at(0)->getEffectivePorosityTerm().value()), 0.2);
}

TEST_F(StandardNodeVDFFixture, getZoneOfSources) {
    ASSERT_EQ((at(0)->getZoneOfSources()), 3);
}

TEST_F(StandardNodeVDFFixture, getZoneOfSinks) {
    ASSERT_EQ((at(0)->getZoneOfSinks()), 0);
}

TEST_F(StandardNodeVDFFixture, getSourcesBelowZeta) {
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);

    ASSERT_EQ((at(0)->getSourcesBelowZeta(0).value()), 852); // sink of 852 m³/d (currently not actually used)
    // Question: why is getSourcesBelowZeta(0) used according to the SWI2 code? (lines 3523-3569, more precisely 3555)
    ASSERT_EQ((at(0)->getSourcesBelowZeta(1).value()), 0);
    ASSERT_EQ((at(0)->getSourcesBelowZeta(2).value()), 0);
    ASSERT_EQ((at(0)->getSourcesBelowZeta(3).value()), 0);
    ASSERT_EQ((at(0)->getSourcesBelowZeta(4).value()), 0);
    // todo test for unconfined node? (changes computation of HCOF)
}

TEST_F(StandardNodeVDFFixture, getPseudoSourceBelowZeta) {
    ASSERT_NEAR((at(0)->getPseudoSourceBelowZeta(0).value()), 0.0046666, 0.0000001);
    ASSERT_EQ((at(0)->getPseudoSourceBelowZeta(1).value()), 0);
    ASSERT_NEAR((at(0)->getPseudoSourceBelowZeta(2).value()), 0.0046666, 0.0000001);
    ASSERT_EQ((at(0)->getPseudoSourceBelowZeta(3).value()), 0); // at bottom -> should be 0
    ASSERT_EQ((at(0)->getPseudoSourceBelowZeta(4).value()), 0); // at bottom -> should be 0
    // todo test head part with two horizontal neighbor nodes that have different heads
}

TEST_F(StandardNodeVDFFixture, getTipToeFlow) {
    ASSERT_EQ((at(0)->getTipToeFlow(0).value()), 0); // no tip and toe tracking for the 0th surface
    ASSERT_NEAR((at(0)->getTipToeFlow(1).value()), 0.0677666, 0.0000001); // only vertical
    ASSERT_NEAR((at(0)->getTipToeFlow(2).value()), 0.0677666, 0.0000001); // only vertical
    ASSERT_EQ((at(0)->getTipToeFlow(3).value()), 0); // not active in node 0
    ASSERT_EQ((at(0)->getTipToeFlow(4).value()), 0); // not active in node 0
    // todo
}

TEST_F(StandardNodeVDFFixture, getZetaRHS) {
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);

    ASSERT_NEAR((at(0)->getZetaRHS(0).value()), -853.9953333, 0.0000001);
    ASSERT_NEAR((at(0)->getZetaRHS(1).value()), -1.7322333, 0.0000001);
    ASSERT_NEAR((at(0)->getZetaRHS(2).value()), -1.4275666, 0.0000001);
    ASSERT_NEAR((at(0)->getZetaRHS(3).value()), 0, 0.0000001);
    ASSERT_NEAR((at(0)->getZetaRHS(4).value()), 0, 0.0000001);
}

TEST_F(StandardNodeVDFFixture, getVDFMatrixEntries) {
    // node 0
    // - localZetaID = 0
    ASSERT_NEAR((at(0)->getVDFMatrixEntries(0)[0].value()), -0.2, 0.0000001);
    ASSERT_NEAR((at(0)->getVDFMatrixEntries(0)[1].value()), 0, 0.0000001);
    // - localZetaID = 1
    ASSERT_NEAR((at(0)->getVDFMatrixEntries(1)[0].value()), -0.2046666, 0.0000001);
    ASSERT_NEAR((at(0)->getVDFMatrixEntries(1)[1].value()), 0.0046666, 0.0000001);
    // node 1 with localZetaID = 2
    ASSERT_NEAR((at(1)->getVDFMatrixEntries(2)[1].value()), -0.2, 0.0000001);
    ASSERT_NEAR((at(1)->getVDFMatrixEntries(2)[0].value()), 0, 0.0000001);
    // node 2 with localZetaID = 3
    ASSERT_NEAR((at(2)->getVDFMatrixEntries(3)[2].value()), -0.2033333, 0.0000001);
    ASSERT_NEAR((at(2)->getVDFMatrixEntries(3)[3].value()), 0.0033333, 0.0000001);
}

TEST_F(StandardNodeVDFFixture, verticalZetaMovement) {
    at(2)->setZeta(3, 0 * si::meter);
    at(2)->setZetaPosInNode(3);
    at(2)->verticalZetaMovement();
    ASSERT_EQ(at(2)->getZeta(0).value(), 0); // surface 0 is at the top of node 0 and 2
    ASSERT_EQ(at(2)->getZeta(1).value(), 0); // surface 1 is between in node 0 and at the top of node 2
    ASSERT_EQ(at(2)->getZeta(2).value(), 0); // surface 2 is at the bottom of node 0 and top of node 2 -> move
    ASSERT_NEAR(at(0)->getZeta(3).value(), 0, 0.0000001); // getFluxCorrTop() < 0 -> do not move zeta 3 in node 0
    ASSERT_NEAR(at(2)->getZeta(3).value(), -0.3388333, 0.0000001); // getFluxCorrTop() < 0 -> lower zeta 3 in node 2
    ASSERT_EQ(at(2)->getZeta(4).value(), -10); // remains unchanged
}

TEST_F(StandardNodeVDFFixture, horizontalZetaMovement) {
    at(3)->setZeta(3, 0 * si::meter);
    at(3)->setZetaPosInNode(3);
    at(2)->horizontalZetaMovement();
    ASSERT_NEAR(at(2)->getZeta(3).value(), -2.49, 0.0000001); // moved zeta 3 in node 2 up
    ASSERT_NEAR(at(3)->getZeta(3).value(), -0.01, 0.0000001); // moved zeta 3 in node 3 down

}

TEST_F(StandardNodeVDFFixture, clipInnerZetas) {
    at(0)->setZeta(1, 20 * si::meter); // set zeta outside the upper bound (Zetas.front())
    at(0)->clipInnerZetas();
    ASSERT_EQ(at(0)->getZeta(1).value(), 10);

    at(0)->setZeta(1, -20 * si::meter);
    at(0)->clipInnerZetas();
    ASSERT_EQ(at(0)->getZeta(3).value(), 0);
}

TEST_F(StandardNodeVDFFixture, correctCrossingZetas) {
    at(0)->setZeta(1, 5 * si::meter); // set zeta 1 below zeta 2 (=7.5m)
    at(0)->correctCrossingZetas();
    ASSERT_EQ(at(0)->getZeta(1).value(), 6.25); // both zetas will be set to their average value
    ASSERT_EQ(at(0)->getZeta(2).value(), 6.25); // both zetas will be set to their average value
}

TEST_F(StandardNodeVDFFixture, preventZetaLocking) {
    at(2)->setZeta(3, -10 * si::meter);
    at(3)->setZeta(3, 0 * si::meter);
    at(2)->setZetaPosInNode(3);
    at(3)->setZetaPosInNode(3);
    at(2)->preventZetaLocking();
    ASSERT_NEAR(at(2)->getZeta(3).value(), -9.9995, 0.0000001); // moved zeta 3 in node 2 up
    ASSERT_NEAR(at(3)->getZeta(3).value(), -0.0005, 0.0000001); // moved zeta 3 in node 3 down
}


