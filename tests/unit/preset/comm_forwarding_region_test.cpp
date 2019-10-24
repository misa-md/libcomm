//
// Created by genshen on 2019/10/24.
//

#include <gtest/gtest.h>
#include <comm/domain/region.hpp>
#include <comm/types_define.h>
#include <comm/preset/comm_forwarding_region.h>

TEST(fw_region_volume_test, fw_region_test) {
    const comm::_type_lattice_size gx = 2, gy = 3, gz = 4;
    const comm::_type_lattice_size ghost_size[comm::DIMENSION_SIZE] = {gx, gy, gz};
    // define a region volume is 10*11*12
    const comm::_type_lattice_size bx = 10, by = 11, bz = 12;
    const comm::Region<comm::_type_lattice_coord> local_box_region{gx, gy, gz, gx + bx, gy + by, gz + bz};

    comm::Region<comm::_type_lattice_size> region_x_l = comm::fwCommLocalRegion(
            ghost_size, local_box_region, comm::DIM_X, comm::DIR_LOWER);
    comm::Region<comm::_type_lattice_size> region_x_h = comm::fwCommLocalRegion(
            ghost_size, local_box_region, comm::DIM_X, comm::DIR_HIGHER);
    EXPECT_EQ(region_x_l.volume(), region_x_h.volume());
    EXPECT_EQ(region_x_l.volume(), gx * by * bz);

    comm::Region<comm::_type_lattice_size> region_y_l = comm::fwCommLocalRegion(
            ghost_size, local_box_region, comm::DIM_Y, comm::DIR_LOWER);
    comm::Region<comm::_type_lattice_size> region_y_h = comm::fwCommLocalRegion(
            ghost_size, local_box_region, comm::DIM_Y, comm::DIR_HIGHER);
    EXPECT_EQ(region_y_l.volume(), region_y_h.volume());
    EXPECT_EQ(region_y_l.volume(), (bx + 2 * gx) * gy * bz);

    comm::Region<comm::_type_lattice_size> region_z_l = comm::fwCommLocalRegion(
            ghost_size, local_box_region, comm::DIM_Z, comm::DIR_LOWER);
    comm::Region<comm::_type_lattice_size> region_z_h = comm::fwCommLocalRegion(
            ghost_size, local_box_region, comm::DIM_Z, comm::DIR_HIGHER);
    EXPECT_EQ(region_z_l.volume(), region_z_h.volume());
    EXPECT_EQ(region_z_l.volume(), (bx + 2 * gx) * (by + 2 * gy) * gz);
}
