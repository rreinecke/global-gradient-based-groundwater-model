#ifndef TESTING_NUMERICS_TEST_HPP
#define TESTING_NUMERICS_TEST_HPP

#include <gtest/gtest.h>
#include "../../../src/Solver/Numerics.hpp"

//Required due to long double matrix support maybe remove in future
using dvector = Eigen::Matrix<long double, Dynamic, 1>;
using longMatrix = Eigen::Matrix<long double, -1, 1, 0, -1, 1>;

class NumericsFixture : public ::testing::Test {
public:
    AdaptiveDamping adp;
    dvector residuals;
    dvector x;

    NumericsFixture() {
        longMatrix m = longMatrix::Random(3, 3);
        residuals = Map<dvector>(m.data(), m.cols() * m.rows());
        m.transposeInPlace();
        x = Map<dvector>(m.data(), m.cols() * m.rows());;
    }
};

TEST_F(NumericsFixture, ConstructorDeath) {
    ASSERT_DEATH(adp.getDamping(residuals, x, false);, "Damping hasn't been properly initialized");
}

TEST_F(NumericsFixture, ComplexConstructor) {
    adp = AdaptiveDamping(0, 0, 0, x);
    dvector ret = adp.getDamping(residuals, x, false);
    ASSERT_TRUE(x.isApprox(x + ret));
}

TEST_F(NumericsFixture, getDamping) {
    FAIL();
}

#endif
