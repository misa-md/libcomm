// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.

#include <algorithm>
#include <logs/logs.h>

#include "domain.h"
#include "utils/mpi_domain.h"
#include "utils/mpi_data_types.h"
#include "pack/pack.h"


Domain::Domain(const std::array<u_int64_t, DIMENSION> _phase_space,
               const double _lattice_const, const double _cutoff_radius_factor)
        : lattice_const(_lattice_const), cutoff_radius_factor(_cutoff_radius_factor),
          cut_lattice(static_cast<int>(ceil(_cutoff_radius_factor))), phase_space(_phase_space),
//      todo _grid_size(0),
//      todo _meas_global_length(0.0),
        /**initialize following references */
          meas_global_length(_meas_global_length),
          grid_size(_grid_size),
          meas_global_box_coord_region(_meas_global_box_coord_region),
          grid_coord_sub_box(_grid_coord_sub_box),
          rank_id_neighbours(_rank_id_neighbours),
          meas_sub_box_region(_meas_sub_box_region),
          meas_ghost_length(_meas_ghost_length),
          meas_ghost_region(_meas_ghost_region),
          lattice_size_sub_box(_lattice_sub_box_size),
          lattice_size_ghost_extended(_lattice_size_ghost_extended),
          lattice_size_ghost(_lattice_size_ghost),
          lattice_coord_sub_box_region(_lattice_coord_sub_box_region),
          lattice_coord_ghost_region(_lattice_coord_ghost_region),
          local_sub_box_lattice_coord_region(_local_sub_box_lattice_coord_region),
          local_ghost_lattice_coord_region(_local_ghost_lattice_coord_region),
          dbx_lattice_size_sub_box(_dbx_lattice_sub_box_size),
          dbx_lattice_size_ghost_extended(_dbx_lattice_size_ghost_extended),
          dbx_lattice_size_ghost(_dbx_lattice_size_ghost),
          dbx_lattice_coord_sub_box_region(_dbx_lattice_coord_sub_box_region),
          dbx_lattice_coord_ghost_region(_dbx_lattice_coord_ghost_region),
          dbx_local_sub_box_lattice_coord_region(_dbx_local_sub_box_lattice_coord_region),
          dbx_local_ghost_lattice_coord_region(_dbx_local_ghost_lattice_coord_region) {}

Domain::Builder &Domain::Builder::setPhaseSpace(const int64_t phaseSpace[DIMENSION]) {
    for (int i = 0; i < DIMENSION; i++) {
        _phase_space[i] = phaseSpace[i]; // todo type not match
    }
    return *this;
}

Domain::Builder &Domain::Builder::setLatticeConst(const double latticeConst) {
    _lattice_const = latticeConst;
    return *this;
}

Domain::Builder &Domain::Builder::setCutoffRadius(const double cutoff_radius_factor) {
    _cutoff_radius_factor = cutoff_radius_factor;
    return *this;
}

Domain::Builder &Domain::Builder::setComm(kiwi::mpi_process mpi_process, MPI_Comm *comm) {
    _mpi_pro = mpi_process;
    _p_comm = comm;
    return *this;
}

void Domain::Builder::decomposition(Domain &domain) {
    // Assume N can be decomposed as N = N_x * N_y * N_z,
    // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
    // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION-1 equals N.
    MPI_Dims_create(_mpi_pro.all_ranks, DIMENSION,
                    domain._grid_size); // fixme origin code: (int *) &domain._grid_size
    kiwi::logs::i(MASTER_PROCESSOR, "decomposition", "MPI grid dimensions: {0},{1},{2}\n",
                  domain._grid_size[0], domain._grid_size[1], domain._grid_size[2]);

    int period[DIMENSION];
    // 3维拓扑
    for (int d = 0; d < DIMENSION; d++) {
        period[d] = 1;
    }
    // sort the processors to fit 3D cartesian topology.
    // the rank id may change.
    MPI_Cart_create(_mpi_pro.comm, DIMENSION, domain._grid_size, period, true, _p_comm);

    kiwi::RID new_rank;
    MPI_Comm_rank(*_p_comm, &new_rank);
    // get cartesian coordinate of current processor.
    MPI_Cart_coords(*_p_comm, new_rank, DIMENSION,
                    domain._grid_coord_sub_box);

    // get the rank ids of contiguous processors of current processor.
    for (int d = 0; d < DIMENSION; d++) {
        MPI_Cart_shift(*_p_comm, d, 1, &domain._rank_id_neighbours[d][LOWER], &domain._rank_id_neighbours[d][HIGHER]);
    }
}

void Domain::Builder::createGlobalDomain(Domain &domain) {
    for (int d = 0; d < DIMENSION; d++) {
        //phaseSpace个单位长度(单位长度即latticeconst)
        domain._meas_global_length[d] = _phase_space[d] * _lattice_const;
        domain._meas_global_box_coord_region.low[d] = 0; // lower bounding is set to 0 by default.
        domain._meas_global_box_coord_region.high[d] = domain._meas_global_length[d];
    }
}

void Domain::Builder::buildLatticeDomain(Domain &domain) {
    // set lattice size of sub-box.
    for (int d = 0; d < DIMENSION; d++) {
        domain._lattice_sub_box_size[d] = _phase_space[d] / domain._grid_size[d] +
                                          (domain._grid_coord_sub_box[d] < (_phase_space[d] % domain._grid_size[d])
                                           ? 1 : 0);
    }

    // set ghost lattice size.
    for (int d = 0; d < DIMENSION; d++) {
        // i * ceil(x) >= ceil(i*x) for all x ∈ R and i ∈ Z
        domain._lattice_size_ghost[d] = domain.cut_lattice;
        domain._lattice_size_ghost_extended[d] = domain._lattice_sub_box_size[d] + 2 * domain._lattice_size_ghost[d];
    }

    // set lattice coordinate boundary of sub-box.
    for (int d = 0; d < DIMENSION; d++) {
        // floor equals to "/" if all operation number >=0.
        domain._lattice_coord_sub_box_region.low[d] =
                domain._grid_coord_sub_box[d] * (_phase_space[d] / domain._grid_size[d]) +
                std::min(domain._grid_coord_sub_box[d], static_cast<int>(_phase_space[d]) % domain._grid_size[d]);
        // todo set measure coord = lower*lattice_const.
        domain._lattice_coord_sub_box_region.high[d] =
                domain._lattice_coord_sub_box_region.low[d] + domain._lattice_sub_box_size[d];
        // (domain._grid_coord_sub_box[d] + 1) * _phase_space[d] / domain._grid_size[d];
    }

    // set lattice coordinate boundary for ghost.
    for (int d = 0; d < DIMENSION; d++) {
        domain._lattice_coord_ghost_region.low[d] = domain._lattice_coord_sub_box_region.low[d] -
                                                    domain._lattice_size_ghost[d];
        domain._lattice_coord_ghost_region.high[d] = domain._lattice_coord_sub_box_region.high[d] +
                                                     domain._lattice_size_ghost[d];
    }

    // set local lattice coordinate.
    for (int d = 0; d < DIMENSION; d++) {
        domain._local_ghost_lattice_coord_region.low[d] = 0;
        domain._local_ghost_lattice_coord_region.high[d] = domain._lattice_size_ghost_extended[d];
        domain._local_sub_box_lattice_coord_region.low[d] = domain._lattice_size_ghost[d];
        domain._local_sub_box_lattice_coord_region.high[d] =
                domain._lattice_size_ghost[d] + domain._lattice_sub_box_size[d];
    }

    // todo double x
    domain._dbx_lattice_sub_box_size[0] = 2 * domain._lattice_sub_box_size[0];
    domain._dbx_lattice_sub_box_size[1] = domain._lattice_sub_box_size[1];
    domain._dbx_lattice_sub_box_size[2] = domain._lattice_sub_box_size[2];

    domain._dbx_lattice_size_ghost[0] = 2 * domain._lattice_size_ghost[0];
    domain._dbx_lattice_size_ghost[1] = domain._lattice_size_ghost[1];
    domain._dbx_lattice_size_ghost[2] = domain._lattice_size_ghost[2];

    domain._dbx_lattice_size_ghost_extended[0] = 2 * domain._lattice_size_ghost_extended[0];
    domain._dbx_lattice_size_ghost_extended[1] = domain._lattice_size_ghost_extended[1];
    domain._dbx_lattice_size_ghost_extended[2] = domain._lattice_size_ghost_extended[2];

    domain._dbx_lattice_coord_sub_box_region.x_low = 2 * domain._lattice_coord_sub_box_region.x_low;
    domain._dbx_lattice_coord_sub_box_region.y_low = domain._lattice_coord_sub_box_region.y_low;
    domain._dbx_lattice_coord_sub_box_region.z_low = domain._lattice_coord_sub_box_region.z_low;
    domain._dbx_lattice_coord_sub_box_region.x_high = 2 * domain._lattice_coord_sub_box_region.x_high;
    domain._dbx_lattice_coord_sub_box_region.y_high = domain._lattice_coord_sub_box_region.y_high;
    domain._dbx_lattice_coord_sub_box_region.z_high = domain._lattice_coord_sub_box_region.z_high;

    domain._dbx_lattice_coord_ghost_region.x_low = 2 * domain._lattice_coord_ghost_region.x_low;
    domain._dbx_lattice_coord_ghost_region.y_low = domain._lattice_coord_ghost_region.y_low;
    domain._dbx_lattice_coord_ghost_region.z_low = domain._lattice_coord_ghost_region.z_low;
    domain._dbx_lattice_coord_ghost_region.x_high = 2 * domain._lattice_coord_ghost_region.x_high;
    domain._dbx_lattice_coord_ghost_region.y_high = domain._lattice_coord_ghost_region.y_high;
    domain._dbx_lattice_coord_ghost_region.z_high = domain._lattice_coord_ghost_region.z_high;

    domain._dbx_local_ghost_lattice_coord_region.x_low = 2 * domain._local_ghost_lattice_coord_region.x_low;
    domain._dbx_local_ghost_lattice_coord_region.y_low = domain._local_ghost_lattice_coord_region.y_low;
    domain._dbx_local_ghost_lattice_coord_region.z_low = domain._local_ghost_lattice_coord_region.z_low;
    domain._dbx_local_ghost_lattice_coord_region.x_high = 2 * domain._local_ghost_lattice_coord_region.x_high;
    domain._dbx_local_ghost_lattice_coord_region.y_high = domain._local_ghost_lattice_coord_region.y_high;
    domain._dbx_local_ghost_lattice_coord_region.z_high = domain._local_ghost_lattice_coord_region.z_high;

    domain._dbx_local_sub_box_lattice_coord_region.x_low = 2 * domain._local_sub_box_lattice_coord_region.x_low;
    domain._dbx_local_sub_box_lattice_coord_region.y_low = domain._local_sub_box_lattice_coord_region.y_low;
    domain._dbx_local_sub_box_lattice_coord_region.z_low = domain._local_sub_box_lattice_coord_region.z_low;
    domain._dbx_local_sub_box_lattice_coord_region.x_high = 2 * domain._local_sub_box_lattice_coord_region.x_high;
    domain._dbx_local_sub_box_lattice_coord_region.y_high = domain._local_sub_box_lattice_coord_region.y_high;
    domain._dbx_local_sub_box_lattice_coord_region.z_high = domain._local_sub_box_lattice_coord_region.z_high;
}

void Domain::Builder::buildMeasuredDomain(Domain &domain) {
    // calculate measured length in each dimension.
    for (int d = 0; d < DIMENSION; d++) {
        // the lower and upper bounding of current sub-box.
        domain._meas_sub_box_region.low[d] = domain._meas_global_box_coord_region.low[d] +
                                             domain._lattice_coord_sub_box_region.low[d] * _lattice_const;
        domain._meas_sub_box_region.high[d] = domain._meas_global_box_coord_region.low[d] +
                                              domain._lattice_coord_sub_box_region.high[d] * _lattice_const;

        domain._meas_ghost_length[d] = _cutoff_radius_factor * _lattice_const; // ghost length todo

        domain._meas_ghost_region.low[d] = domain._meas_sub_box_region.low[d] - domain._meas_ghost_length[d];
        domain._meas_ghost_region.high[d] = domain._meas_sub_box_region.high[d] + domain._meas_ghost_length[d];
    }
}

Domain *Domain::Builder::build() {
    Domain *p_domain = new Domain(_phase_space, _lattice_const, _cutoff_radius_factor);
    decomposition(*p_domain);
    createGlobalDomain(*p_domain);
    buildLatticeDomain(*p_domain);
    buildMeasuredDomain(*p_domain);
    return p_domain;
}

Domain *Domain::Builder::localBuild(const int _grid_size[DIMENSION], const int _grid_coord[DIMENSION]) {
    Domain *p_domain = new Domain(_phase_space, _lattice_const, _cutoff_radius_factor);
    for (int i = 0; i < 3; i++) {
        p_domain->_grid_size[i] = _grid_size[i];
        p_domain->_grid_coord_sub_box[i] = _grid_coord[i];
    }
    createGlobalDomain(*p_domain);
    buildLatticeDomain(*p_domain);
    buildMeasuredDomain(*p_domain);
    return p_domain;
}
