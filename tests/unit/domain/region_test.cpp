//
// Created by genshen on 2019-04-15.
//

#include <gtest/gtest.h>
#include <comm/domain/region.hpp>

TEST(region_ref_test, region_test) {
    comm::Region<double> region(-1, -2, -3, 1, 2, 3);
    EXPECT_EQ(region.low[0], -1);
    EXPECT_EQ(region.low[1], -2);
    EXPECT_EQ(region.low[2], -3);
    EXPECT_EQ(region.high[0], 1);
    EXPECT_EQ(region.high[1], 2);
    EXPECT_EQ(region.high[2], 3);
}

// test assignment of region class
TEST(region_assign_test, region_test) {
    comm::Region<int> region(-1, -2, -3, 1, 2, 3);
    comm::Region<int> region_ass;
    region_ass = region;
    EXPECT_EQ(region_ass.low[0], -1);
    EXPECT_EQ(region_ass.low[1], -2);
    EXPECT_EQ(region_ass.low[2], -3);
    EXPECT_EQ(region_ass.high[0], 1);
    EXPECT_EQ(region_ass.high[1], 2);
    EXPECT_EQ(region_ass.high[2], 3);

    EXPECT_EQ(&(region_ass.x_low), &(region_ass.low[0]));
    EXPECT_EQ(&(region_ass.y_low), &(region_ass.low[1]));
    EXPECT_EQ(&(region_ass.z_low), &(region_ass.low[2]));
    EXPECT_EQ(&(region_ass.x_high), &(region_ass.high[0]));
    EXPECT_EQ(&(region_ass.y_high), &(region_ass.high[1]));
    EXPECT_EQ(&(region_ass.z_high), &(region_ass.high[2]));

    region_ass = comm::Region<int>(1, 2, 3, -1, -2, -3);
    EXPECT_EQ(region_ass.low[0], 1);
    EXPECT_EQ(region_ass.low[1], 2);
    EXPECT_EQ(region_ass.low[2], 3);
    EXPECT_EQ(region_ass.high[0], -1);
    EXPECT_EQ(region_ass.high[1], -2);
    EXPECT_EQ(region_ass.high[2], -3);

    EXPECT_EQ(&(region_ass.x_low), &(region_ass.low[0]));
    EXPECT_EQ(&(region_ass.y_low), &(region_ass.low[1]));
    EXPECT_EQ(&(region_ass.z_low), &(region_ass.low[2]));
    EXPECT_EQ(&(region_ass.x_high), &(region_ass.high[0]));
    EXPECT_EQ(&(region_ass.y_high), &(region_ass.high[1]));
    EXPECT_EQ(&(region_ass.z_high), &(region_ass.high[2]));
}

TEST(region_assign_scope_test, region_test) {
    comm::Region<int> region;
    {
        region = comm::Region<int>{-1, -2, -3, 1, 2, 3};
    }

    EXPECT_EQ(region.low[0], -1);
    EXPECT_EQ(region.low[1], -2);
    EXPECT_EQ(region.low[2], -3);
    EXPECT_EQ(region.high[0], 1);
    EXPECT_EQ(region.high[1], 2);
    EXPECT_EQ(region.high[2], 3);
}

std::vector<comm::Region<int> > returnRegionTest() {
    return {
            {-1, -2, -3, 1, 2, 3}
    };
}

TEST(region_assign_scope2_test, region_test) {
    std::vector<comm::Region<int> > regions = returnRegionTest();

    EXPECT_EQ(regions[0].x_low, regions[0].low[0]);
    EXPECT_EQ(regions[0].y_low, regions[0].low[1]);
    EXPECT_EQ(regions[0].z_low, regions[0].low[2]);
    EXPECT_EQ(regions[0].x_high, regions[0].high[0]);
    EXPECT_EQ(regions[0].y_high, regions[0].high[1]);
    EXPECT_EQ(regions[0].z_high, regions[0].high[2]);

    EXPECT_EQ(&(regions[0].x_low), &(regions[0].low[0]));
    EXPECT_EQ(&(regions[0].y_low), &(regions[0].low[1]));
    EXPECT_EQ(&(regions[0].z_low), &(regions[0].low[2]));
    EXPECT_EQ(&(regions[0].x_high), &(regions[0].high[0]));
    EXPECT_EQ(&(regions[0].y_high), &(regions[0].high[1]));
    EXPECT_EQ(&(regions[0].z_high), &(regions[0].high[2]));
}
