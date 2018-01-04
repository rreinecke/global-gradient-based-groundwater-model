//#include "Integration.cpp"
#include "Unit/UnitTests.hpp"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    auto result(RUN_ALL_TESTS());
    //::testing::internal::TimeInMillis elapsed(
    //        ::testing::UnitTest::GetInstance()->elapsed_time());
    //std::cerr << "Elapsed Time: " << elapsed << "\n";
    //ASSERT_LT(elapsed, measurePerf ? 180 * 1000 : 215 * 1000);
    return result;
}
