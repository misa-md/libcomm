//
// Created by genshen on 2019/10/24.
//

#include <comm/domain/region.hpp>
#include <comm/preset/comm_forwarding_region.h>
#include <comm/types_define.h>
#include <gtest/gtest.h>

TEST(fw_send_region_volume_test, fw_region_test) {
  const comm::_type_lattice_size gx = 2, gy = 3, gz = 4;
  const comm::_type_lattice_size ghost_size[comm::DIMENSION_SIZE] = {gx, gy, gz};
  // define a region volume is 10*11*12
  const comm::_type_lattice_size bx = 10, by = 11, bz = 12;
  const comm::Region<comm::_type_lattice_coord> local_box_region{gx, gy, gz, gx + bx, gy + by, gz + bz};

  comm::Region<comm::_type_lattice_size> region_x_l =
      comm::fwCommLocalSendRegion(ghost_size, local_box_region, comm::DIM_X, comm::DIR_LOWER);
  comm::Region<comm::_type_lattice_size> region_x_h =
      comm::fwCommLocalSendRegion(ghost_size, local_box_region, comm::DIM_X, comm::DIR_HIGHER);
  EXPECT_EQ(region_x_l.volume(), region_x_h.volume());
  EXPECT_EQ(region_x_l.volume(), gx * by * bz);

  comm::Region<comm::_type_lattice_size> region_y_l =
      comm::fwCommLocalSendRegion(ghost_size, local_box_region, comm::DIM_Y, comm::DIR_LOWER);
  comm::Region<comm::_type_lattice_size> region_y_h =
      comm::fwCommLocalSendRegion(ghost_size, local_box_region, comm::DIM_Y, comm::DIR_HIGHER);
  EXPECT_EQ(region_y_l.volume(), region_y_h.volume());
  EXPECT_EQ(region_y_l.volume(), (bx + 2 * gx) * gy * bz);

  comm::Region<comm::_type_lattice_size> region_z_l =
      comm::fwCommLocalSendRegion(ghost_size, local_box_region, comm::DIM_Z, comm::DIR_LOWER);
  comm::Region<comm::_type_lattice_size> region_z_h =
      comm::fwCommLocalSendRegion(ghost_size, local_box_region, comm::DIM_Z, comm::DIR_HIGHER);
  EXPECT_EQ(region_z_l.volume(), region_z_h.volume());
  EXPECT_EQ(region_z_l.volume(), (bx + 2 * gx) * (by + 2 * gy) * gz);
}

comm::Region<comm::_type_lattice_size>
TEST_receive_region_imp(const comm::_type_lattice_size ghost_size[comm::DIMENSION_SIZE],
                        const comm::Region<comm::_type_lattice_coord> local_box_region, const unsigned int dimension,
                        const unsigned int direction) {
  switch (dimension << 2 | direction) {
  case comm::DIM_X << 2 | comm::DIR_LOWER: { // x dimension, lower direction
    comm::_type_lattice_size xstart = local_box_region.x_high;
    comm::_type_lattice_size ystart = local_box_region.y_low;
    comm::_type_lattice_size zstart = local_box_region.z_low;
    comm::_type_lattice_size xstop = local_box_region.x_high + ghost_size[0];
    comm::_type_lattice_size ystop = local_box_region.y_high;
    comm::_type_lattice_size zstop = local_box_region.z_high;
    return comm::Region<comm::_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
  }
  case comm::DIM_X << 2 | comm::DIR_HIGHER: { // x dimension, higher direction
    comm::_type_lattice_size xstart = local_box_region.x_low - ghost_size[0];
    comm::_type_lattice_size ystart = local_box_region.y_low;
    comm::_type_lattice_size zstart = local_box_region.z_low;
    comm::_type_lattice_size xstop = local_box_region.x_low;
    comm::_type_lattice_size ystop = local_box_region.y_high;
    comm::_type_lattice_size zstop = local_box_region.z_high;
    return comm::Region<comm::_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
  }
  case comm::DIM_Y << 2 | comm::DIR_LOWER: { // y dimension, lower direction
    comm::_type_lattice_size xstart = local_box_region.x_low - ghost_size[0];
    comm::_type_lattice_size ystart = local_box_region.y_high;
    comm::_type_lattice_size zstart = local_box_region.z_low;
    comm::_type_lattice_size xstop = local_box_region.x_high + ghost_size[0];
    comm::_type_lattice_size ystop = local_box_region.y_high + ghost_size[1];
    comm::_type_lattice_size zstop = local_box_region.z_high;
    return comm::Region<comm::_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
  }
  case comm::DIM_Y << 2 | comm::DIR_HIGHER: { // y dimension, higher direction
    comm::_type_lattice_size xstart = local_box_region.x_low - ghost_size[0];
    comm::_type_lattice_size ystart = local_box_region.y_low - ghost_size[1];
    comm::_type_lattice_size zstart = local_box_region.z_low;
    comm::_type_lattice_size xstop = local_box_region.x_high + ghost_size[0];
    comm::_type_lattice_size ystop = local_box_region.y_low;
    comm::_type_lattice_size zstop = local_box_region.z_high;
    return comm::Region<comm::_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
  }
  case comm::DIM_Z << 2 | comm::DIR_LOWER: { // z dimension, lower direction
    comm::_type_lattice_size xstart = local_box_region.x_low - ghost_size[0];
    comm::_type_lattice_size ystart = local_box_region.y_low - ghost_size[1];
    comm::_type_lattice_size zstart = local_box_region.z_high;
    comm::_type_lattice_size xstop = local_box_region.x_high + ghost_size[0];
    comm::_type_lattice_size ystop = local_box_region.y_high + ghost_size[1];
    comm::_type_lattice_size zstop = local_box_region.z_high + ghost_size[2];
    return comm::Region<comm::_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
  }
  case comm::DIM_Z << 2 | comm::DIR_HIGHER: { // z dimension, higher direction
    comm::_type_lattice_size xstart = local_box_region.x_low - ghost_size[0];
    comm::_type_lattice_size ystart = local_box_region.y_low - ghost_size[1];
    comm::_type_lattice_size zstart = local_box_region.z_low - ghost_size[2];
    comm::_type_lattice_size xstop = local_box_region.x_high + ghost_size[0];
    comm::_type_lattice_size ystop = local_box_region.y_high + ghost_size[1];
    comm::_type_lattice_size zstop = local_box_region.z_low;
    return comm::Region<comm::_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
  }
  default:
    // this case is not allowed.
    assert(false);
  }
}

TEST(fw_recv_region_imp_comparing, fw_region_test) {
  const comm::_type_lattice_size gx = 2, gy = 3, gz = 4;
  const comm::_type_lattice_size ghost_size[comm::DIMENSION_SIZE] = {gx, gy, gz};
  // define a region volume is 10*11*12
  const comm::_type_lattice_size bx = 10, by = 11, bz = 12;
  const comm::Region<comm::_type_lattice_coord> local_box_region{gx, gy, gz, gx + bx, gy + by, gz + bz};

  // compare volume of each dimension and direction
  for (int dim = comm::DIM_X; dim <= comm::DIM_Z; dim++) {
    for (int dir = comm::DIR_LOWER; dir <= comm::DIR_HIGHER; dir++) {
      comm::Region<comm::_type_lattice_size> region_recv =
          comm::fwCommLocalRecvRegion(ghost_size, local_box_region, dim, dir);
      comm::Region<comm::_type_lattice_size> region_recv_table =
          TEST_receive_region_imp(ghost_size, local_box_region, dim, dir);
      EXPECT_EQ(region_recv.x_low, region_recv_table.x_low);
      EXPECT_EQ(region_recv.y_low, region_recv_table.y_low);
      EXPECT_EQ(region_recv.z_low, region_recv_table.z_low);
      EXPECT_EQ(region_recv.x_high, region_recv_table.x_high);
      EXPECT_EQ(region_recv.y_high, region_recv_table.y_high);
      EXPECT_EQ(region_recv.z_high, region_recv_table.z_high);
    }
  }
}

TEST(fw_region_volume_test, fw_region_test) {
  const comm::_type_lattice_size gx = 2, gy = 3, gz = 4;
  const comm::_type_lattice_size ghost_size[comm::DIMENSION_SIZE] = {gx, gy, gz};
  // define a region volume is 10*11*12
  const comm::_type_lattice_size bx = 10, by = 11, bz = 12;
  const comm::Region<comm::_type_lattice_coord> local_box_region{gx, gy, gz, gx + bx, gy + by, gz + bz};

  // compare volume of each dimension and direction
  for (int dim = comm::DIM_X; dim <= comm::DIM_Z; dim++) {
    for (int dir = comm::DIR_LOWER; dir <= comm::DIR_HIGHER; dir++) {
      comm::Region<comm::_type_lattice_size> region_send =
          comm::fwCommLocalSendRegion(ghost_size, local_box_region, dim, dir);
      comm::Region<comm::_type_lattice_size> region_recv =
          comm::fwCommLocalRecvRegion(ghost_size, local_box_region, dim, dir);
      EXPECT_EQ(region_send.volume(), region_recv.volume());
    }
  }
}
