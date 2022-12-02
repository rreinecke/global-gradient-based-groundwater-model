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

TEST_F(StandardNodeVDFFixture, setZetaPosInZone) {

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
    ASSERT_EQ((at(0)->getZoneConductances(got)[0].value()), 0); // todo compute correct result
}

TEST_F(StandardNodeVDFFixture, getRHS) {
    at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);

    ASSERT_EQ((at(0)->getRHS().value()), 0); // todo calculate correct result
}

TEST_F(StandardNodeVDFFixture, getZetaRHS) {
    at(0)->getZetaRHS(1);
}

TEST_F(StandardNodeVDFFixture, getEffectivePorosityTerm) {

}

TEST_F(StandardNodeVDFFixture, getZetaConductannce) {

}

TEST_F(StandardNodeVDFFixture, getNusTop) {

}

TEST_F(StandardNodeVDFFixture, getNusBot) {

}

TEST_F(StandardNodeVDFFixture, getSourceTermBelowZeta) {

}

TEST_F(StandardNodeVDFFixture, getTipToeFlow) {

}

TEST_F(StandardNodeVDFFixture, getZoneConductanceCum) {

}

TEST_F(StandardNodeVDFFixture, getFlowPseudoSource) {

}

TEST_F(StandardNodeVDFFixture, getVerticalFluxCorrection) {

}

TEST_F(StandardNodeVDFFixture, getFluxCorrTop) {

}

TEST_F(StandardNodeVDFFixture, getFluxCorrDown) {

}

TEST_F(StandardNodeVDFFixture, clipTopZetasToHeads) {

}

TEST_F(StandardNodeVDFFixture, verticalZetaMovement) {

}

TEST_F(StandardNodeVDFFixture, horizontalZetaMovement) {

}

TEST_F(StandardNodeVDFFixture, clipZetaHeights) {

}

TEST_F(StandardNodeVDFFixture, correctCrossingZetas) {

}

TEST_F(StandardNodeVDFFixture, preventZetaLocking) {

}

TEST_F(StandardNodeVDFFixture, addInitialZeta) {

}

TEST_F(StandardNodeVDFFixture, setZeta) {

}

TEST_F(StandardNodeVDFFixture, getZetas) {

}

TEST_F(StandardNodeVDFFixture, getZetasChange) {

}

TEST_F(StandardNodeVDFFixture, updateZetaChange) {

}

TEST_F(StandardNodeVDFFixture, setZoneOfSinksAndSources) {

}

TEST_F(StandardNodeVDFFixture, setEffectivePorosity) {

}

TEST_F(StandardNodeVDFFixture, setEffectivePorosity_direct) {

}