#include <gtest/gtest.h>
#include "../../../src/Model/Node.hpp"

using NodeVector = std::shared_ptr<std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>>;

class StandardNodeFixture : public ::testing::Test {
public:
    NodeVector nodes;

    void SetUp() {
        DensityProperties densityProperties = GlobalFlow::Model::DensityProperties::setDensityProperties(true, {1000.0, 1012.5, 1025.0}, 0.2, 0.2);

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
        nodes->at(0)->addInitialZeta(10.0 * si::meter, 0);
        nodes->at(1)->addInitialZeta(10.0 * si::meter, 0);
        nodes->at(2)->addInitialZeta(0.0 * si::meter, 0);
        nodes->at(3)->addInitialZeta(0.0 * si::meter, 0);

        // zeta surface 1 (in nodes 0 & 1 between, in nodes 2 & 3 at top)
        nodes->at(0)->addInitialZeta(7.5 * si::meter, 0);
        nodes->at(1)->addInitialZeta(2.5 * si::meter, 1);
        nodes->at(2)->addInitialZeta(0.0 * si::meter, 2);
        nodes->at(3)->addInitialZeta(0.0 * si::meter, 3);

        // zeta surface 2 (in nodes 0 & 1 at bottom, in nodes 2 & 3 between)
        nodes->at(0)->addInitialZeta(0.0 * si::meter, 0);
        nodes->at(1)->addInitialZeta(0.0 * si::meter, 1);
        nodes->at(2)->addInitialZeta(-2.5 * si::meter, 2);
        nodes->at(3)->addInitialZeta(-7.5 * si::meter, 3);

        // zeta surface 3 (all at bottom of nodes)
        nodes->at(0)->addInitialZeta(0.0 * si::meter, 0);
        nodes->at(1)->addInitialZeta(0.0 * si::meter, 1);
        nodes->at(2)->addInitialZeta(-10.0 * si::meter, 2);
        nodes->at(3)->addInitialZeta(-10.0 * si::meter, 3);
    }

    using p_node = std::unique_ptr<GlobalFlow::Model::NodeInterface>;

    p_node &at(int pos) { return nodes->at(pos); }
};

TEST_F(StandardNodeFixture, setZetaPosInZone) {

}

TEST_F(StandardNodeFixture, getRHSConstantDensity) {
    at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);
    ASSERT_EQ((at(0)->getRHSConstantDensity().value()), -50);
}

TEST_F(StandardNodeFixture, getRHS) {
    at(0)->setHead_direct(1);
    at(0)->addExternalFlow(RECHARGE, 0, 50, 0);
    at(0)->addExternalFlow(RIVER, 1 * si::meter, 50, 1 * si::meter);
    at(0)->addExternalFlow(WETLAND, 1 * si::meter, 50, 1 * si::meter);

    // todo calculate correct result
    ASSERT_EQ((at(0)->getRHS().value()), 0);
}

TEST_F(StandardNodeFixture, getZetaRHS) {

}

TEST_F(StandardNodeFixture, getEffectivePorosityTerm) {

}

TEST_F(StandardNodeFixture, getZetaConductannce) {

}

TEST_F(StandardNodeFixture, getNusTop) {

}

TEST_F(StandardNodeFixture, getNusBot) {

}

TEST_F(StandardNodeFixture, getSourceTermBelowZeta) {

}

TEST_F(StandardNodeFixture, getTipToeFlow) {

}

TEST_F(StandardNodeFixture, getZoneConductance) {

}

TEST_F(StandardNodeFixture, getZoneConductanceCum) {

}

TEST_F(StandardNodeFixture, getFlowPseudoSource) {

}

TEST_F(StandardNodeFixture, getVerticalFluxCorrection) {

}

TEST_F(StandardNodeFixture, getFluxCorrTop) {

}

TEST_F(StandardNodeFixture, getFluxCorrDown) {

}

TEST_F(StandardNodeFixture, clipTopZetasToHeads) {

}

TEST_F(StandardNodeFixture, verticalZetaMovement) {

}

TEST_F(StandardNodeFixture, horizontalZetaMovement) {

}

TEST_F(StandardNodeFixture, clipZetaHeights) {

}

TEST_F(StandardNodeFixture, correctCrossingZetas) {

}

TEST_F(StandardNodeFixture, preventZetaLocking) {

}

TEST_F(StandardNodeFixture, addInitialZeta) {

}

TEST_F(StandardNodeFixture, setZeta) {

}

TEST_F(StandardNodeFixture, getZetas) {

}

TEST_F(StandardNodeFixture, getZetasChange) {

}

TEST_F(StandardNodeFixture, updateZetaChange) {

}

TEST_F(StandardNodeFixture, setZoneOfSinksAndSources) {

}

TEST_F(StandardNodeFixture, setEffectivePorosity) {

}

TEST_F(StandardNodeFixture, setEffectivePorosity_direct) {

}