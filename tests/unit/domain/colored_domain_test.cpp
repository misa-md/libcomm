//
// Created by genshen on 2019-08-11.
//

#include <gtest/gtest.h>
#include <comm/domain/colored_domain.h>

TEST(sector_size_test, colored_domain_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {50 * grid_size[0], 60 * grid_size[1], 71 * grid_size[2]};
    const double lattice_const = 0.86;
    const double cutoff_radius_factor = 1.1421;

    comm::ColoredDomain *p_domain = comm::ColoredDomain::Builder()
            .setPhaseSpace(space)
            .setGhostSize(5)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);

    // each domain is a 50*60*72 box split to 8 sectors.

    // test sector 0,0,0
    EXPECT_EQ(p_domain->sector_lattice_size[0][0], 50 / 2);
    EXPECT_EQ(p_domain->sector_lattice_size[0][1], 60 / 2);
    EXPECT_EQ(p_domain->sector_lattice_size[0][2], 71 / 2);

    // test sector 0,1,1 which is X_LOW|Y_HIGH|Z_HIGH = 6
    const int s = comm::X_LOW | comm::Y_HIGH | comm::Z_HIGH;
    EXPECT_EQ(p_domain->sector_lattice_size[s][0], 50 / 2);
    EXPECT_EQ(p_domain->sector_lattice_size[s][1], 60 / 2);
    EXPECT_EQ(p_domain->sector_lattice_size[s][2], 71 / 2 + 1);

    // test sector 0,0,0 and 1,1,1
    EXPECT_EQ(p_domain->sector_lattice_size[0][0] + p_domain->sector_lattice_size[7][0], 50);
    EXPECT_EQ(p_domain->sector_lattice_size[0][1] + p_domain->sector_lattice_size[7][1], 60);
    EXPECT_EQ(p_domain->sector_lattice_size[0][2] + p_domain->sector_lattice_size[7][2], 71);
}

TEST(local_sector_region_test, colored_domain_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {50 * grid_size[0], 60 * grid_size[1], 71 * grid_size[2]};
    const double lattice_const = 0.86;
    const double cutoff_radius_factor = 1.1421;
    const int ghost_size = 5;

    comm::ColoredDomain *p_domain = comm::ColoredDomain::Builder()
            .setPhaseSpace(space)
            .setGhostSize(ghost_size)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);

    // test sector 0,0,0
    EXPECT_EQ(p_domain->local_sector_region[0].x_low, ghost_size);
    EXPECT_EQ(p_domain->local_sector_region[0].y_low, ghost_size);
    EXPECT_EQ(p_domain->local_sector_region[0].z_low, ghost_size);
    EXPECT_EQ(p_domain->local_sector_region[0].x_high, ghost_size + 50 / 2); // from ghost_size + [0, 50/2)
    EXPECT_EQ(p_domain->local_sector_region[0].y_high, ghost_size + 60 / 2); // from ghost_size + [0, 60/2)
    EXPECT_EQ(p_domain->local_sector_region[0].z_high, ghost_size + 71 / 2); // from ghost_size + [0, 71/2)

    // test sector 0,1,1 which is X_LOW|Y_HIGH|Z_HIGH = 6
    const int s = comm::X_LOW | comm::Y_HIGH | comm::Z_HIGH;
    EXPECT_EQ(p_domain->local_sector_region[s].x_low, ghost_size);
    EXPECT_EQ(p_domain->local_sector_region[s].y_low, ghost_size + 60 / 2);
    EXPECT_EQ(p_domain->local_sector_region[s].z_low, ghost_size + 71 / 2);
    EXPECT_EQ(p_domain->local_sector_region[s].x_high, ghost_size + 50 / 2); // from ghost_size + [0, 50/2)
    EXPECT_EQ(p_domain->local_sector_region[s].y_high, ghost_size + 60); // from ghost_size + [60/2, 60)
    EXPECT_EQ(p_domain->local_sector_region[s].z_high, ghost_size + 71); // from ghost_size + [71/2, 71)

    // test sector 0,0,0 and 1,1,1: higher boundary of 0,0,0 equals lower boundary of 1,1,1
    EXPECT_EQ(p_domain->local_sector_region[0].x_high, p_domain->local_sector_region[7].x_low);
    EXPECT_EQ(p_domain->local_sector_region[0].y_high, p_domain->local_sector_region[7].y_low);
    EXPECT_EQ(p_domain->local_sector_region[0].z_high, p_domain->local_sector_region[7].z_low);
}
