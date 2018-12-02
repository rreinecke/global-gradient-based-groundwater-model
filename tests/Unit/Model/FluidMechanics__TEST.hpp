#include <gtest/gtest.h>
#include "../../../src/Model/FluidMechanics.hpp"

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
    quantity<Meter> head_neig = 1 * si::meter;
    quantity<Meter> head_self = 1 * si::meter;
    quantity<Meter> ele_neig = 1 * si::meter;
    quantity<Meter> ele_self = 1 * si::meter;
    quantity<Meter> deltaV_neig = 1 * si::meter;
    quantity<Meter> deltaV_self = 1 * si::meter;
    bool confined = false;

    FlowInputHor t = std::make_tuple(k_neig, k_self, edgeLength_neig, edgeLength_self,
                                     head_neig, head_self, ele_neig, ele_self, deltaV_neig,
                                     deltaV_self, confined);

    //TODO add complexer examples
    ASSERT_DOUBLE_EQ(m.calculateHarmonicMeanConductance(t).value(), 0.1);
}

TEST(FluidMechanics, smoothFunction__NWT) {
    //quantity<Meter> elevation, quantity<Meter> verticalSize, quantity<Meter> head
    FluidMechanics m = FluidMechanics();
    //TODO more coverage
    ASSERT_DOUBLE_EQ((m.smoothFunction__NWT(1 * si::meter, 1 * si::meter, 1 * si::meter)), 1);
    ASSERT_DOUBLE_EQ((m.smoothFunction__NWT(10 * si::meter, 1 * si::meter, 1 * si::meter)), 10e-21);
}

TEST(FluidMechanics, getHCOF) {
    FluidMechanics m = FluidMechanics();
    bool steadyState = true;
    quantity<Dimensionless> stepModifier = 1 * si::si_dimensionless;
    quantity<SquareMeter> storageCapacity = 1 * si::square_meter;
    quantity<MeterSquaredPerTime> P = 2 * si::square_meter / day;

    ASSERT_DOUBLE_EQ((m.getHCOF(steadyState, stepModifier, storageCapacity, P).value()), 2);
    steadyState = false;
    ASSERT_DOUBLE_EQ((m.getHCOF(steadyState, stepModifier, storageCapacity, P).value()), 1);
}

TEST(FluidMechanics, calculateVerticalConductance) {
    FluidMechanics m = FluidMechanics();
    quantity<Velocity> k_vert_neig = 0.1 * si::meter / day;
    quantity<Velocity> k_vert_self = 0.1 * si::meter / day;
    quantity<Meter> verticalSize_self = 1 * si::meter;
    quantity<Meter> head_self = 1 * si::meter;
    quantity<Meter> elevation_self = 1 * si::meter;
    quantity<SquareMeter> area_self = 1 * si::square_meter;
    quantity<Meter> elevation_neig = 1 * si::meter;
    quantity<Meter> depth_neig = 1 * si::meter;
    quantity<Meter> head_neig = 1 * si::meter;
    bool confined = false;


    FlowInputVert t = std::make_tuple(k_vert_neig, k_vert_self, verticalSize_self, head_self, elevation_self, area_self,
                                     elevation_neig, depth_neig, head_neig,
                                     confined);
    ASSERT_DOUBLE_EQ(m.calculateVerticalConductance(t).value(), 0.1);
}

TEST(FluidMechanics, getDerivate__NWT) {
    FluidMechanics m = FluidMechanics();
    ASSERT_DOUBLE_EQ(m.getDerivate__NWT(10 * si::meter, 5 * si::meter, 5 * si::meter), 0);
    ASSERT_DOUBLE_EQ(m.getDerivate__NWT(10 * si::meter, 5 * si::meter, 10 * si::meter), 0.2);
}
