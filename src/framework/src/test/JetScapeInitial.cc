// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...
#include "../InitialCondition.h"
#include "gtest/gtest.h"

TEST(JetscapeInitialTest, TEST_SAMPLE){
    auto ini = JetScapeInitial("auau200", 0, 5, 10, 0.2);

    EXPECT_EQ(ini.entropy_density_distribution_.size(), 10000);
    double mul = ini.info_.total_entropy;
    ASSERT_TRUE(mul >= 100 && mul <= 150);
}
