// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...
#include "../TrentoInitial.h"
#include "gtest/gtest.h"
#include <cstdlib>

using namespace Jetscape;

TEST(JetscapeInitialTest, TEST_SAMPLE){
    double cent_low = 30.0;
    double cent_high = 40.0;
    double stored_slow = 33.4898493386; 
    double stored_shigh = 51.671394987;
    double grid_max = 10.0;
    double grid_step = 0.2;
    for (int k = 0; k < 100; k++ ) {
        auto ini = TrentoInitial();
        unsigned random_seed = std::rand();
        ini.pre_defined("auau200", cent_low, cent_high,
                               grid_max, grid_step, random_seed);
        EXPECT_EQ(ini.entropy_density_distribution_.size(), 10000);
        double mul = ini.info_.total_entropy;
        ASSERT_TRUE(mul >= stored_slow && mul <= stored_shigh);

        EXPECT_EQ(ini.get_x_size(), 100);
        EXPECT_EQ(ini.get_y_size(), 100);
        EXPECT_EQ(ini.get_z_size(), 1);

        auto idx = ini.coord_from_idx(150);
        EXPECT_EQ(std::get<0>(idx), -10 + 50 * ini.get_x_step());
        EXPECT_EQ(std::get<1>(idx), -10 + 1 * ini.get_y_step());
        EXPECT_EQ(std::get<2>(idx), 0 * ini.get_z_step());
    }
}

