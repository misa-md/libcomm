//
// Created by genshen on 2018-05-02.
//

#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <mpi.h>

#include "domain_test_utils.h"

// @MPI
TEST(domain_test_decomposition, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;

    comm::Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int grid_sum[3] = {0, 0, 0};
    int local_grid[3] = {(_domain->lattice_size_sub_box[0]),
                         (_domain->lattice_size_sub_box[1]),
                         (_domain->lattice_size_sub_box[2])};
    // update rank and all ranks.
    comm::mpi_process m_process{}; // not set comm
    MPI_Comm_size(MPI_COMM_WORLD, &m_process.all_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_process.own_rank);

    MPI_Reduce(local_grid, grid_sum, 3, MPI_INT, MPI_SUM, comm::COMM_MASTER, MPI_COMM_WORLD);

    grid_sum[0] /= m_process.all_ranks / _domain->grid_size[0];
    grid_sum[1] /= m_process.all_ranks / _domain->grid_size[1];
    grid_sum[2] /= m_process.all_ranks / _domain->grid_size[2];
    if (m_process.own_rank == comm::COMM_MASTER) {
        EXPECT_EQ(grid_sum[0], space[0]);
        EXPECT_EQ(grid_sum[1], space[1]);
        EXPECT_EQ(grid_sum[2], space[2]);
    }
    delete _domain;
}

// @MPI
TEST(domain_local_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    comm::Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    int nlocalx = floor(_domain->meas_sub_box_region.x_high / (lattice_const)) -
                  floor(_domain->meas_sub_box_region.x_low / (lattice_const));
    int nlocaly = floor(_domain->meas_sub_box_region.y_high / lattice_const) -
                  floor(_domain->meas_sub_box_region.y_low / lattice_const);
    int nlocalz = floor(_domain->meas_sub_box_region.z_high / lattice_const) -
                  floor(_domain->meas_sub_box_region.z_low / lattice_const);

    EXPECT_EQ(nlocalx, _domain->lattice_size_sub_box[0]);
    EXPECT_EQ(nlocaly, _domain->lattice_size_sub_box[1]);
    EXPECT_EQ(nlocalz, _domain->lattice_size_sub_box[2]);

    delete _domain;
}

// @MPI
TEST(domain_ghost_lattice_size, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    comm::Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    // ceil(_domain->meas_ghost_length[0] / lattice_const) equals to ceil(_cutoff_radius_factor)
    EXPECT_EQ(_domain->lattice_size_ghost[0], ceil(cutoff_radius_factor));
    EXPECT_EQ(_domain->lattice_size_ghost[1], ceil(cutoff_radius_factor));
    EXPECT_EQ(_domain->lattice_size_ghost[2], ceil(cutoff_radius_factor));

    int nghostx = _domain->lattice_size_sub_box[0] + 2 * ceil(_domain->meas_ghost_length[0] / lattice_const);
    int nghosty = _domain->lattice_size_sub_box[1] + 2 * ceil(_domain->meas_ghost_length[1] / lattice_const);
    int nghostz = _domain->lattice_size_sub_box[2] + 2 * ceil(_domain->meas_ghost_length[2] / lattice_const);

    EXPECT_EQ(nghostx, _domain->lattice_size_ghost_extended[0]);
    EXPECT_EQ(nghosty, _domain->lattice_size_ghost_extended[1]);
    EXPECT_EQ(nghostz, _domain->lattice_size_ghost_extended[2]);

    delete _domain;
}

// @MPI
TEST(domain_local_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    comm::Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);
    // test delta values

    EXPECT_EQ(_domain->lattice_size_sub_box[0],
              _domain->local_sub_box_lattice_coord_region.x_high - _domain->local_sub_box_lattice_coord_region.x_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[1],
              _domain->local_sub_box_lattice_coord_region.y_high - _domain->local_sub_box_lattice_coord_region.y_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[2],
              _domain->local_sub_box_lattice_coord_region.z_high - _domain->local_sub_box_lattice_coord_region.z_low);

    EXPECT_EQ(_domain->lattice_size_ghost_extended[0],
              _domain->local_ghost_lattice_coord_region.x_high - _domain->local_ghost_lattice_coord_region.x_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[1],
              _domain->local_ghost_lattice_coord_region.y_high - _domain->local_ghost_lattice_coord_region.y_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[2],
              _domain->local_ghost_lattice_coord_region.z_high - _domain->local_ghost_lattice_coord_region.z_low);
    // test double x dalta values
    EXPECT_EQ(_domain->dbx_lattice_size_sub_box[0],
              _domain->dbx_local_sub_box_lattice_coord_region.x_high -
              _domain->dbx_local_sub_box_lattice_coord_region.x_low);
    EXPECT_EQ(_domain->dbx_lattice_size_sub_box[1],
              _domain->dbx_local_sub_box_lattice_coord_region.y_high -
              _domain->dbx_local_sub_box_lattice_coord_region.y_low);
    EXPECT_EQ(_domain->dbx_lattice_size_sub_box[2],
              _domain->dbx_local_sub_box_lattice_coord_region.z_high -
              _domain->dbx_local_sub_box_lattice_coord_region.z_low);

    EXPECT_EQ(_domain->dbx_lattice_size_ghost_extended[0],
              _domain->dbx_local_ghost_lattice_coord_region.x_high -
              _domain->dbx_local_ghost_lattice_coord_region.x_low);
    EXPECT_EQ(_domain->dbx_lattice_size_ghost_extended[1],
              _domain->dbx_local_ghost_lattice_coord_region.y_high -
              _domain->dbx_local_ghost_lattice_coord_region.y_low);
    EXPECT_EQ(_domain->dbx_lattice_size_ghost_extended[2],
              _domain->dbx_local_ghost_lattice_coord_region.z_high -
              _domain->dbx_local_ghost_lattice_coord_region.z_low);

}

// @MPI
TEST(domain_global_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    comm::Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of local sub-box
    int lolocalx = floor(_domain->meas_sub_box_region.x_low / lattice_const);
    int lolocaly = floor(_domain->meas_sub_box_region.y_low / lattice_const);
    int lolocalz = floor(_domain->meas_sub_box_region.z_low / lattice_const);

    EXPECT_EQ(lolocalx, _domain->lattice_coord_sub_box_region.x_low);
    EXPECT_EQ(lolocaly, _domain->lattice_coord_sub_box_region.y_low);
    EXPECT_EQ(lolocalz, _domain->lattice_coord_sub_box_region.z_low);

    // upper boundary of lattice coordinate of local sub-box
    int uplocalx = floor(_domain->meas_sub_box_region.x_high / lattice_const);
    int uplocaly = floor(_domain->meas_sub_box_region.y_high / lattice_const);
    int uplocalz = floor(_domain->meas_sub_box_region.z_high / lattice_const);

    EXPECT_EQ(uplocalx, _domain->lattice_coord_sub_box_region.x_high);
    EXPECT_EQ(uplocaly, _domain->lattice_coord_sub_box_region.y_high);
    EXPECT_EQ(uplocalz, _domain->lattice_coord_sub_box_region.z_high);

    // test delta values
    EXPECT_EQ(_domain->lattice_size_sub_box[0],
              _domain->lattice_coord_sub_box_region.x_high - _domain->lattice_coord_sub_box_region.x_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[1],
              _domain->lattice_coord_sub_box_region.y_high - _domain->lattice_coord_sub_box_region.y_low);
    EXPECT_EQ(_domain->lattice_size_sub_box[2],
              _domain->lattice_coord_sub_box_region.z_high - _domain->lattice_coord_sub_box_region.z_low);
    delete _domain;
}

// @MPI
TEST(domain_ghost_lattice_coord, domain_test) {
    int64_t space[3] = {50, 60, 72};
    double lattice_const = 0.86;
    double cutoff_radius_factor = 1.1421;
    comm::Domain *_domain = getDomainInstance(space, lattice_const, cutoff_radius_factor);

    // lower boundary of lattice coordinate of ghost
    int loghostx = _domain->lattice_coord_sub_box_region.x_low - ceil(cutoff_radius_factor / lattice_const);
    int loghosty = _domain->lattice_coord_sub_box_region.y_low - ceil(cutoff_radius_factor / lattice_const);
    int loghostz = _domain->lattice_coord_sub_box_region.z_low - ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(loghostx, _domain->lattice_coord_ghost_region.x_low);
    EXPECT_EQ(loghosty, _domain->lattice_coord_ghost_region.y_low);
    EXPECT_EQ(loghostz, _domain->lattice_coord_ghost_region.z_low);

    // upper boundary of lattice coordinate of ghost
    int upghostx = _domain->lattice_coord_sub_box_region.x_high + ceil(cutoff_radius_factor / lattice_const);
    int upghosty = _domain->lattice_coord_sub_box_region.y_high + ceil(cutoff_radius_factor / lattice_const);
    int upghostz = _domain->lattice_coord_sub_box_region.z_high + ceil(cutoff_radius_factor / lattice_const);

    EXPECT_EQ(upghostx, _domain->lattice_coord_ghost_region.x_high);
    EXPECT_EQ(upghosty, _domain->lattice_coord_ghost_region.y_high);
    EXPECT_EQ(upghostz, _domain->lattice_coord_ghost_region.z_high);

    // test delta values
    EXPECT_EQ(_domain->lattice_size_ghost_extended[0],
              _domain->lattice_coord_ghost_region.x_high - _domain->lattice_coord_ghost_region.x_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[1],
              _domain->lattice_coord_ghost_region.y_high - _domain->lattice_coord_ghost_region.y_low);
    EXPECT_EQ(_domain->lattice_size_ghost_extended[2],
              _domain->lattice_coord_ghost_region.z_high - _domain->lattice_coord_ghost_region.z_low);
    delete _domain;
}

// @MPI
TEST(domain_local_build_test, domain_test) {
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

    EXPECT_EQ(p_domain->lattice_size_sub_box[0], 50); // space in this sub-box
    EXPECT_EQ(p_domain->lattice_size_sub_box[1], 60); // space in this sub-box
    EXPECT_EQ(p_domain->lattice_size_sub_box[2], 72); // space in this sub-box
    delete p_domain;
}
