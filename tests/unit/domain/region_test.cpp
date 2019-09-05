//
// Created by genshen on 2019-04-15.
//

#include <gtest/gtest.h>
#include <domain/region.hpp>

TEST(region_ref_test, region_test) {
    comm::Region<double> region(-1, -2, -3, 1, 2, 3);
    EXPECT_EQ(region.x_low, -1);
    EXPECT_EQ(region.y_low, -2);
    EXPECT_EQ(region.z_low, -3);
    EXPECT_EQ(region.x_high, 1);
    EXPECT_EQ(region.y_high, 2);
    EXPECT_EQ(region.z_high, 3);
}

// test assignment of region class
TEST(region_assign_test, region_test) {
    comm::Region<int> region(-1, -2, -3, 1, 2, 3);
    comm::Region<int> region_ass;
    region_ass = region;
    EXPECT_EQ(region_ass.x_low, -1);
    EXPECT_EQ(region_ass.y_low, -2);
    EXPECT_EQ(region_ass.z_low, -3);
    EXPECT_EQ(region_ass.x_high, 1);
    EXPECT_EQ(region_ass.y_high, 2);
    EXPECT_EQ(region_ass.z_high, 3);

    region_ass = comm::Region<int>(1, 2, 3, -1, -2, -3);
    EXPECT_EQ(region_ass.x_low, 1);
    EXPECT_EQ(region_ass.y_low, 2);
    EXPECT_EQ(region_ass.z_low, 3);
    EXPECT_EQ(region_ass.x_high, -1);
    EXPECT_EQ(region_ass.y_high, -2);
    EXPECT_EQ(region_ass.z_high, -3);
}
