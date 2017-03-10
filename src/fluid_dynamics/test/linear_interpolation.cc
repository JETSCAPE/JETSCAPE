#include "../linear_interpolation.h"
#include "gtest/gtest.h"

TEST(LinearInterpolationTest, TEST_TRUE){
    EXPECT_EQ(0.5, linear_int(0.0, 1.0, 0.0, 1.0, 0.5));
}
