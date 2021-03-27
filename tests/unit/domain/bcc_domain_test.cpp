//
// Created by genshen on 2019-04-20.
//

#include <comm/domain/bcc_domain.h>
#include <comm/domain/domain.h>
#include <gtest/gtest.h>

TEST(bcc_domain_test_from_domain, bcc_domain_test) {
  const int grid_size[3] = {2, 2, 2};
  const int grid_coord[3] = {0, 0, 0};
  const int64_t space[3] = {50 * grid_size[0], 60 * grid_size[1], 72 * grid_size[2]};
  const double lattice_const = 0.86;
  const double cutoff_radius_factor = 1.1421;
  comm::Domain *p_domain = comm::Domain::Builder()
                               .setPhaseSpace(space)
                               .setCutoffRadius(cutoff_radius_factor)
                               .setLatticeConst(lattice_const)
                               .localBuild(grid_size, grid_coord);

  comm::Domain p_domain2 = *p_domain;
  // test domain assign
  EXPECT_EQ(p_domain2.cutoff_radius_factor, p_domain->cutoff_radius_factor);
  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    EXPECT_EQ(p_domain2.grid_size[d], p_domain->grid_size[d]);
  }

  comm::BccDomain p_bcc_domain(*p_domain);
  EXPECT_EQ(p_bcc_domain.cutoff_radius_factor, p_domain->cutoff_radius_factor);
  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    EXPECT_EQ(p_bcc_domain.grid_size[d], p_domain->grid_size[d]);
    EXPECT_EQ(p_bcc_domain.ghost_ext_lattice_region.low[d], p_domain->ghost_ext_lattice_region.low[d]);
    EXPECT_EQ(p_bcc_domain.ghost_ext_lattice_region.high[d], p_domain->ghost_ext_lattice_region.high[d]);
  }
}

TEST(bcc_domain_test_dbx, bcc_domain_test) {
  const int grid_size[3] = {2, 2, 2};
  const int grid_coord[3] = {0, 0, 0};
  const int64_t space[3] = {50 * grid_size[0], 60 * grid_size[1], 72 * grid_size[2]};
  const double lattice_const = 0.86;
  const double cutoff_radius_factor = 1.1421;
  comm::Domain *p_domain = comm::Domain::Builder()
                               .setPhaseSpace(space)
                               .setCutoffRadius(cutoff_radius_factor)
                               .setLatticeConst(lattice_const)
                               .localBuild(grid_size, grid_coord);

  comm::BccDomain p_bcc_domain(*p_domain);
  EXPECT_EQ(p_bcc_domain.cutoff_radius_factor, p_domain->cutoff_radius_factor);
  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    int times = d == 0 ? 2 : 1;
    EXPECT_EQ(p_bcc_domain.dbx_lattice_size_ghost[d], times * p_domain->lattice_size_ghost[d]);
    EXPECT_EQ(p_bcc_domain.dbx_ghost_extended_lattice_size[d], times * p_domain->ghost_extended_lattice_size[d]);
    EXPECT_EQ(p_bcc_domain.dbx_sub_box_lattice_size[d], times * p_domain->sub_box_lattice_size[d]);

    EXPECT_EQ(p_bcc_domain.dbx_ghost_ext_lattice_region.low[d], times * p_domain->ghost_ext_lattice_region.low[d]);
    EXPECT_EQ(p_bcc_domain.dbx_ghost_ext_lattice_region.high[d], times * p_domain->ghost_ext_lattice_region.high[d]);

    EXPECT_EQ(p_bcc_domain.dbx_local_ghost_ext_lattice_region.low[d],
              times * p_domain->local_ghost_ext_lattice_region.low[d]);
    EXPECT_EQ(p_bcc_domain.dbx_local_ghost_ext_lattice_region.high[d],
              times * p_domain->local_ghost_ext_lattice_region.high[d]);

    EXPECT_EQ(p_bcc_domain.dbx_local_sub_box_lattice_region.low[d],
              times * p_domain->local_sub_box_lattice_region.low[d]);
    EXPECT_EQ(p_bcc_domain.dbx_local_sub_box_lattice_region.high[d],
              times * p_domain->local_sub_box_lattice_region.high[d]);

    EXPECT_EQ(p_bcc_domain.dbx_sub_box_lattice_region.low[d], times * p_domain->sub_box_lattice_region.low[d]);
    EXPECT_EQ(p_bcc_domain.dbx_sub_box_lattice_region.high[d], times * p_domain->sub_box_lattice_region.high[d]);
  }
}
