// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...
#include "../InitialCondition.h"
#include "gtest/gtest.h"

TEST(JetscapeInitialTest, TEST_SAMPLE){
    double cent_low = 0.0;
    double cent_high = 10.0;
    double stored_slow = 112.36181125;
    double stored_shigh = 205.62898671;
    double grid_max = 10.0;
    double grid_step = 0.2;
    auto ini = JetScapeInitial("auau200", cent_low, cent_high,
                               grid_max, grid_step);

    EXPECT_EQ(ini.entropy_density_distribution_.size(), 10000);
    double mul = ini.info_.total_entropy;
    ASSERT_TRUE(mul >= stored_slow && mul <= stored_shigh);
}
