//
// Created by genshen on 2019-04-19.
//

#include "builder.h"


template<typename B, typename D>
B &Builder<B, D>::setPhaseSpace(const int64_t phaseSpace[DIMENSION_SIZE]) {
    for (int i = 0; i < DIMENSION_SIZE; i++) {
        _phase_space[i] = phaseSpace[i]; // todo type not match
    }
    return *static_cast<B *>(this);
}

template<typename B, typename D>
B &Builder<B, D>::setLatticeConst(const double latticeConst) {
    _lattice_const = latticeConst;
    return *static_cast<B *>(this);
}

template<typename B, typename D>
B &Builder<B, D>::setGhostSize(const unsigned int ghost_size) {
    _ghost_size = ghost_size;
    return *static_cast<B *>(this);
}

template<typename B, typename D>
B &Builder<B, D>::setCutoffRadius(const double cutoff_radius_factor) {
    _cutoff_radius_factor = cutoff_radius_factor;
    if (_ghost_size == 0) {
        // if ghost size not set, we set it as cut_lattice
        _ghost_size = static_cast<int>(ceil(_cutoff_radius_factor));
    }
    return *static_cast<B *>(this);
}

template<typename B, typename D>
B &Builder<B, D>::setComm(comm::mpi_process mpi_process, MPI_Comm *comm) {
    _mpi_pro = mpi_process;
    _p_comm = comm;
    return *static_cast<B *>(this);
}

template<typename B, typename D>
void Builder<B, D>::decomposition(D &domain) {
    // Assume N can be decomposed as N = N_x * N_y * N_z,
    // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
    // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION_SIZE-1 equals N.
    MPI_Dims_create(_mpi_pro.all_ranks, DIMENSION_SIZE, domain._grid_size);

    int period[DIMENSION_SIZE];
    // 3维拓扑
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        period[d] = 1;
    }
    // sort the processors to fit 3D cartesian topology.
    // the rank id may change.
    MPI_Cart_create(_mpi_pro.comm, DIMENSION_SIZE, domain._grid_size, period, true, _p_comm);

    comm::_MPI_Rank new_rank;
    MPI_Comm_rank(*_p_comm, &new_rank);
    // get cartesian coordinate of current processor.
    MPI_Cart_coords(*_p_comm, new_rank, DIMENSION_SIZE,
                    domain._grid_coord);

    // get the rank ids of contiguous processors of current processor.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        MPI_Cart_shift(*_p_comm, d, 1, &domain._rank_id_neighbours[d][DIR_LOWER],
                       &domain._rank_id_neighbours[d][DIR_HIGHER]);
    }
}

template<typename B, typename D>
void Builder<B, D>::createGlobalDomain(D &domain) {
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        //phaseSpace个单位长度(单位长度即latticeconst)
        domain._meas_global_length[d] = _phase_space[d] * _lattice_const;
        domain._meas_global_region.low[d] = 0; // lower bounding is set to 0 by default.
        domain._meas_global_region.high[d] = domain._meas_global_length[d];
    }
}

template<typename B, typename D>
void Builder<B, D>::buildLatticeDomain(D &domain) {
    // set lattice size of sub-box.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        domain._sub_box_lattice_size[d] = _phase_space[d] / domain._grid_size[d] +
                                          (domain._grid_coord[d] < (_phase_space[d] % domain._grid_size[d])
                                           ? 1 : 0);
    }

    // set ghost lattice size.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        // i * ceil(x) >= ceil(i*x) for all x ∈ R and i ∈ Z
        // add additional one lattice to make all neighbours can be fount in ghost area.
        domain._lattice_size_ghost[d] = _ghost_size;
        domain._ghost_extended_lattice_size[d] = domain._sub_box_lattice_size[d] + 2 * domain._lattice_size_ghost[d];
    }

    // set lattice coordinate boundary of sub-box.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        // floor equals to "/" if all operation number >=0.
        domain._sub_box_lattice_region.low[d] =
                domain._grid_coord[d] * (_phase_space[d] / domain._grid_size[d]) +
                std::min(domain._grid_coord[d], static_cast<int>(_phase_space[d]) % domain._grid_size[d]);
        // todo set measure coord = lower*lattice_const.
        domain._sub_box_lattice_region.high[d] =
                domain._sub_box_lattice_region.low[d] + domain._sub_box_lattice_size[d];
        // (domain._grid_coord[d] + 1) * _phase_space[d] / domain._grid_size[d];
    }

    // set lattice coordinate boundary for ghost.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        domain._ghost_ext_lattice_region.low[d] = domain._sub_box_lattice_region.low[d] -
                                                  domain._lattice_size_ghost[d];
        domain._ghost_ext_lattice_region.high[d] = domain._sub_box_lattice_region.high[d] +
                                                   domain._lattice_size_ghost[d];
    }

    // set local lattice coordinate.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        domain._local_ghost_ext_lattice_region.low[d] = 0;
        domain._local_ghost_ext_lattice_region.high[d] = domain._ghost_extended_lattice_size[d];
        domain._local_sub_box_lattice_region.low[d] = domain._lattice_size_ghost[d];
        domain._local_sub_box_lattice_region.high[d] =
                domain._lattice_size_ghost[d] + domain._sub_box_lattice_size[d];
    }
}

template<typename B, typename D>
void Builder<B, D>::buildMeasuredDomain(D &domain) {
    // calculate measured length in each dimension.
    for (int d = 0; d < DIMENSION_SIZE; d++) {
        // the lower and upper bounding of current sub-box.
        domain._meas_sub_box_region.low[d] = domain._meas_global_region.low[d] +
                                             domain._sub_box_lattice_region.low[d] * _lattice_const;
        domain._meas_sub_box_region.high[d] = domain._meas_global_region.low[d] +
                                              domain._sub_box_lattice_region.high[d] * _lattice_const;

        domain._meas_ghost_length[d] = domain._lattice_size_ghost[d] * _lattice_const; // ghost length fixme

        domain._meas_ghost_ext_region.low[d] = domain._meas_sub_box_region.low[d] - domain._meas_ghost_length[d];
        domain._meas_ghost_ext_region.high[d] = domain._meas_sub_box_region.high[d] + domain._meas_ghost_length[d];
    }
}
