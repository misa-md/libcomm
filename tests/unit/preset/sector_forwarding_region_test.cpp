//
// Created by genshen on 2019/9/5.
//
#include <gtest/gtest.h>
#include <comm/domain/colored_domain.h>
#include <comm/preset/sector_forwarding_region.h>

// test region area of forwarding communication of 8 sectors
TEST(sector_fw_region_area_test, sector_fw_region_test) {
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int local_space[3] = {50, 60, 70};
    const int64_t space[3] = {local_space[0] * grid_size[0],
                              local_space[1] * grid_size[1],
                              local_space[2] * grid_size[2]};
    const double lattice_const = 0.86;
    const double cutoff_radius_factor = 1.1421;
    const comm::_type_lattice_size ghost_size = 5;

    comm::ColoredDomain *p_domain = comm::ColoredDomain::Builder()
            .setPhaseSpace(space)
            .setGhostSize(ghost_size)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);

    const comm::_type_lattice_size g[comm::DIMENSION_SIZE] = {ghost_size, ghost_size, ghost_size};
    for (int s = 0; s < 8; s++) {
        for (int dim = 0; dim < 3; dim++) {
            comm::type_region_array regions_send = comm::fwCommSectorSendRegion(s, dim, g,
                                                                                p_domain->local_split_coord,
                                                                                p_domain->local_sub_box_lattice_region);
            comm::type_region_array regions_recv = comm::fwCommSectorRecvRegion(s, dim, g,
                                                                                p_domain->local_split_coord,
                                                                                p_domain->local_sub_box_lattice_region);

            EXPECT_EQ(regions_send.size(), regions_recv.size());
            comm::_type_lattice_size region_vol_total = 0;
            for (size_t i = 0; i < regions_send.size(); i++) {
                comm::_type_lattice_size region_vol_send = (regions_send[i].x_high - regions_send[i].x_low)
                                                           * (regions_send[i].y_high - regions_send[i].y_low)
                                                           * (regions_send[i].z_high - regions_send[i].z_low);
                comm::_type_lattice_size region_vol_recv = (regions_recv[i].x_high - regions_recv[i].x_low)
                                                           * (regions_recv[i].y_high - regions_recv[i].y_low)
                                                           * (regions_recv[i].z_high - regions_recv[i].z_low);
                EXPECT_EQ(region_vol_send, region_vol_recv);
                region_vol_total += region_vol_send;
            }
            // for convenience, we only test communication region volume if size of 8 sectors are all equal.
            if (local_space[0] % 2 == 0 && local_space[1] % 2 == 0 && local_space[2] % 2 == 0) {
                if (dim == comm::DIM_X) {
                    EXPECT_EQ(region_vol_total, ghost_size * (local_space[1] / 2 + 2 * ghost_size) *
                                                (local_space[2] / 2 + 2 * ghost_size));
                } else if (dim == comm::DIM_Y) {
                    EXPECT_EQ(region_vol_total, ghost_size * (local_space[0] / 2 + 2 * ghost_size) *
                                                (local_space[2] / 2 + 2 * ghost_size));
                } else {
                    EXPECT_EQ(region_vol_total, ghost_size * (local_space[0] / 2 + 2 * ghost_size) *
                                                (local_space[1] / 2 + 2 * ghost_size));
                }
            }
        }
    }
}
