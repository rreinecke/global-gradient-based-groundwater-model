//#include "Integration.cpp"
#include "Unit/UnitTests.hpp"

int main(int argc, char **argv) {
    ::testing::InitGoogleMock(&argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    auto result(RUN_ALL_TESTS());
    ::testing::FLAGS_gtest_death_test_style = "fast";
    //::testing::internal::TimeInMillis elapsed(
    //        ::testing::UnitTest::GetInstance()->elapsed_time());
    //std::cerr << "Elapsed Time: " << elapsed << "\n";
    //ASSERT_LT(elapsed, measurePerf ? 180 * 1000 : 215 * 1000);
    return result;
}
