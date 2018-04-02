/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "../linear_interpolation.h"
#include "gtest/gtest.h"

using namespace Jetscape;

TEST(LinearInterpolationTest, TEST_TRUE){
    EXPECT_EQ(0.5, linear_int(0.0, 1.0, 0.0, 1.0, 0.5));
}

// test code when the type of y is int
TEST(LinearInterpolationTest, TEST_INT){
    EXPECT_EQ(0, linear_int(0.0, 1.0, 0, 1, 0.5));
}
