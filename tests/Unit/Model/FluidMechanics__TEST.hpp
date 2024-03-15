#include <gtest/gtest.h>
#include "../../../src/Model/FluidMechanics.hpp"
#include <tuple>

TEST(FluidMechanics, estimateConductance) {
    FluidMechanics m = FluidMechanics();
    //t_vel K, t_meter length, t_meter width, t_meter Daq, t_meter G, depth
    ASSERT_DOUBLE_EQ(m.estimateConductance(1 * si::meter / day, 100 * si::meter, 10 * si::meter, 100 * si::meter,
                                           1000 * si::meter, 10 * si::meter).value(), 32.855614163202617);
    /*ASSERT_DOUBLE_EQ(m.estimateConductance(0.01 * si::meter / day, 1000 * si::meter, 100 * si::meter, 100 * si::meter,
                                           1000 * si::meter, 10 * si::meter).value(), 10.949907990410804);*/
    testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(m.estimateConductance(0.01 * si::meter / day, 1000 * si::meter, 100 * si::meter, 100 * si::meter,
                                       1000 * si::meter, 1000 * si::meter).value(),
                 "Depth of river and thickness of cell are not fit for this equation");
    ASSERT_DEATH(m.estimateConductance(0 * si::meter / day, 0 * si::meter, 0 * si::meter, 0 * si::meter,
                                       0 * si::meter, 0 * si::meter).value(), "Inputs can't be 0!");
    ASSERT_DEATH(m.estimateConductance(0 * si::meter / day, 0 * si::meter, 0 * si::meter, 1 * si::meter,
                                       1 * si::meter, 0 * si::meter).value(), "Inputs can't be 0!");
}

TEST(FluidMechanics, calcDeltaV) {
    FluidMechanics m = FluidMechanics();
    //quantity<Meter> head, quantity<Meter> elevation, quantity<Meter> depth
    ASSERT_EQ(m.calcDeltaV(11 * si::meter, 20 * si::meter, 10 * si::meter).value(), 1);
    ASSERT_EQ(m.calcDeltaV(5 * si::meter, 2 * si::meter, 5 * si::meter).value(), 5);
    ASSERT_EQ(m.calcDeltaV(5 * si::meter, 20 * si::meter, 10 * si::meter).value(), 0);
}

TEST(FluidMechanics, calculateHarmonicMeanConductance) {
    FluidMechanics m = FluidMechanics();
    quantity<Velocity> k_neig = 0.1 * si::meter / day;
    quantity<Velocity> k_self = 0.1 * si::meter / day;
    quantity<Meter> edgeLength_neig = 1 * si::meter;
    quantity<Meter> edgeLength_self = 1 * si::meter;
    quantity<Meter> edgeWidth_self = 1 * si::meter;
    quantity<Meter> head_neig = 1 * si::meter;
    quantity<Meter> head_self = 1 * si::meter;
    quantity<Meter> ele_neig = 1 * si::meter;
    quantity<Meter> ele_self = 1 * si::meter;
    quantity<Meter> deltaV_neig = 1 * si::meter;
    quantity<Meter> deltaV_self = 1 * si::meter;
    bool confined = false;

    FlowInputHor t = std::make_tuple(k_neig, k_self, edgeLength_neig, edgeLength_self, edgeWidth_self,
                                     head_neig, head_self, ele_neig, ele_self, deltaV_neig,
                                     deltaV_self, confined);

    //TODO add more complex examples
    ASSERT_DOUBLE_EQ(m.calculateHarmonicMeanConductance(t).value(), 0.1);
}

TEST(FluidMechanics, calculateVerticalConductance) {
    FluidMechanics m = FluidMechanics();
    quantity<Velocity> k_vert_neig = 0.1 * si::meter / day;
    quantity<Velocity> k_vert_self = 0.1 * si::meter / day;
    quantity<Meter> verticalSize_self = 1 * si::meter;
    quantity<Meter> verticalSize_neig = 1 * si::meter;
    quantity<Meter> head_self = 1 * si::meter;
    quantity<Meter> elevation_self = 1 * si::meter;
    quantity<SquareMeter> area_self = 1 * si::square_meter;
    quantity<Meter> elevation_neig = 1 * si::meter;
    quantity<Meter> depth_neig = 1 * si::meter;
    quantity<Meter> head_neig = 1 * si::meter;
    bool confined = false;

    FlowInputVert t = std::make_tuple(k_vert_neig, k_vert_self, verticalSize_self, verticalSize_neig, head_self,
                                      head_neig, elevation_self, elevation_neig, area_self, confined);
    ASSERT_DOUBLE_EQ(m.calculateVerticalConductance(t).value(), 0.1);
}

