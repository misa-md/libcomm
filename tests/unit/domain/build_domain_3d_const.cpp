//
// Created by genshen on 2026/6/20.
//

#include <comm/domain/domain.h>
#include <gtest/gtest.h>

TEST(domain_3d_lat_const_measured_test, domain_local_build_test) {
  const int grid_size[3] = {2, 2, 2};
  const int grid_coord[3] = {0, 0, 0};
  const std::array<int, 3> lat_size_3d = {50, 60, 72};
  const int64_t space[3] = {lat_size_3d[0] * grid_size[0], lat_size_3d[1] * grid_size[1],
                            lat_size_3d[2] * grid_size[2]};
  const std::array<double, 3> lattice_const = {1.0, 2.0, 5.0};

  constexpr double cutoff_radius_factor = 1.1421;
  constexpr double ghost_size = 2;

  comm::Domain *p_domain = comm::Domain::Builder()
                               .setPhaseSpace(space)
                               .setCutoffRadius(cutoff_radius_factor, lattice_const[0])
                               .setLatticeConst(lattice_const)
                               .setGhostSize(ghost_size)
                               .localBuild(grid_size, grid_coord);

  EXPECT_EQ(p_domain->sub_box_lattice_size[0], lat_size_3d[0]); // space in this sub-box
  EXPECT_EQ(p_domain->sub_box_lattice_size[1], lat_size_3d[1]); // space in this sub-box
  EXPECT_EQ(p_domain->sub_box_lattice_size[2], lat_size_3d[2]); // space in this sub-box

  // check the measure size.
  // the ghost externed box region for current sub-box
  EXPECT_EQ(p_domain->meas_ghost_ext_region.x_low, -ghost_size * lattice_const[0]);
  EXPECT_EQ(p_domain->meas_ghost_ext_region.y_low, -ghost_size * lattice_const[1]);
  EXPECT_EQ(p_domain->meas_ghost_ext_region.z_low, -ghost_size * lattice_const[2]);
  EXPECT_EQ(p_domain->meas_ghost_ext_region.x_high, (lat_size_3d[0] + ghost_size) * lattice_const[0]);
  EXPECT_EQ(p_domain->meas_ghost_ext_region.y_high, (lat_size_3d[1] + ghost_size) * lattice_const[1]);
  EXPECT_EQ(p_domain->meas_ghost_ext_region.z_high, (lat_size_3d[2] + ghost_size) * lattice_const[2]);

  // the global box length
  EXPECT_EQ(p_domain->meas_global_length[0], grid_size[0] * lat_size_3d[0] * lattice_const[0]);
  EXPECT_EQ(p_domain->meas_global_length[1], grid_size[1] * lat_size_3d[1] * lattice_const[1]);
  EXPECT_EQ(p_domain->meas_global_length[2], grid_size[2] * lat_size_3d[2] * lattice_const[2]);

  // the global box region
  EXPECT_EQ(p_domain->meas_global_region.x_low, 0.0);
  EXPECT_EQ(p_domain->meas_global_region.y_low, 0.0);
  EXPECT_EQ(p_domain->meas_global_region.z_low, 0.0);
  EXPECT_EQ(p_domain->meas_global_region.x_high, grid_size[0] * lat_size_3d[0] * lattice_const[0]);
  EXPECT_EQ(p_domain->meas_global_region.y_high, grid_size[1] * lat_size_3d[1] * lattice_const[1]);
  EXPECT_EQ(p_domain->meas_global_region.z_high, grid_size[2] * lat_size_3d[2] * lattice_const[2]);

  // ghost region length
  EXPECT_EQ(p_domain->meas_ghost_length[0], ghost_size * lattice_const[0]);
  EXPECT_EQ(p_domain->meas_ghost_length[1], ghost_size * lattice_const[1]);
  EXPECT_EQ(p_domain->meas_ghost_length[2], ghost_size * lattice_const[2]);

  // the region of the sub-box
  EXPECT_EQ(p_domain->meas_sub_box_region.x_low, 0);
  EXPECT_EQ(p_domain->meas_sub_box_region.y_low, 0);
  EXPECT_EQ(p_domain->meas_sub_box_region.z_low, 0);
  EXPECT_EQ(p_domain->meas_sub_box_region.x_high, 1 * lat_size_3d[0] * lattice_const[0]);
  EXPECT_EQ(p_domain->meas_sub_box_region.y_high, 1 * lat_size_3d[1] * lattice_const[1]);
  EXPECT_EQ(p_domain->meas_sub_box_region.z_high, 1 * lat_size_3d[2] * lattice_const[2]);

  delete p_domain;
}

// todo: add test for `domain_3d_lat_const_test` with MPI rank coord [1, 1, 1]