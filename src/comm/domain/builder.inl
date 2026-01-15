//
// Created by genshen on 2019-04-19.
//

#include "builder.h"
#include "comm/topology/cart_3d_node_affinity.h"

#include <algorithm>
#include <stdexcept>
#include <stdio.h>
#include <string>

template <typename B, typename D> B &comm::Builder<B, D>::setPhaseSpace(const int64_t phaseSpace[DIMENSION_SIZE]) {
  for (int i = 0; i < DIMENSION_SIZE; i++) {
    _phase_space[i] = phaseSpace[i]; // todo type not match
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D> B &comm::Builder<B, D>::setLatticeConst(const double latticeConst) {
  _lattice_const = latticeConst;
  return *static_cast<B *>(this);
}

template <typename B, typename D> B &comm::Builder<B, D>::setGhostSize(const unsigned int ghost_size) {
  _ghost_size = ghost_size;
  return *static_cast<B *>(this);
}

template <typename B, typename D> B &comm::Builder<B, D>::setCutoffRadius(const double cutoff_radius_factor) {
  _cutoff_radius_factor = cutoff_radius_factor;
  if (_ghost_size == 0) {
    // if ghost size not set, we set it as cut_lattice
    _ghost_size = static_cast<int>(ceil(_cutoff_radius_factor));
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D> B &comm::Builder<B, D>::setComm(comm::mpi_process mpi_process, MPI_Comm *comm) {
  _mpi_pro = mpi_process;
  _p_comm = comm;
  return *static_cast<B *>(this);
}

template <typename B, typename D>
B &comm::Builder<B, D>::setMPIMap3dSubDim(const int mpi_map_3d_sub_dim[DIMENSION_SIZE]) {
  for (int i = 0; i < DIMENSION_SIZE; i++) {
    _mpi_map_3d_sub_dim[i] = mpi_map_3d_sub_dim[i];
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D> void comm::Builder<B, D>::decomposition(D &domain) {
  if (_mpi_map_3d_sub_dim[0] == 1 && _mpi_map_3d_sub_dim[1] == 1 && _mpi_map_3d_sub_dim[2] == 1) {
    decomposition_imp_v1(domain);
  } else {
    // note1: for inner node 1x1x1, it can also call imp_v2.
    // but for the given mpi process, it may produce different mpi rank when comparing with imp_v1.
    decomposition_imp_v2(domain);
  }
}

template <typename B, typename D> void comm::Builder<B, D>::decomposition_imp_v1(D &domain) {
  // Assume N can be decomposed as N = N_x * N_y * N_z,
  // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
  // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION_SIZE-1 equals N.
  MPI_Dims_create(_mpi_pro.all_ranks, DIMENSION_SIZE, domain._grid_size);

  int period[DIMENSION_SIZE];
  // 3维拓扑
  for (int &d : period) {
    d = 1;
  }
  // sort the processors to fit 3D cartesian topology.
  // the rank id may change.
  MPI_Cart_create(_mpi_pro.comm, DIMENSION_SIZE, domain._grid_size, period, true, _p_comm);

  comm::_MPI_Rank new_rank;
  MPI_Comm_rank(*_p_comm, &new_rank);
  // get cartesian coordinate of current processor.
  MPI_Cart_coords(*_p_comm, new_rank, DIMENSION_SIZE, domain._grid_coord);

  // get the rank ids of contiguous processors of current processor.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    MPI_Cart_shift(*_p_comm, d, 1, &domain._rank_id_neighbours[d][DIR_LOWER],
                   &domain._rank_id_neighbours[d][DIR_HIGHER]);
  }
}

template <typename B, typename D> void comm::Builder<B, D>::decomposition_imp_v2(D &domain) {
  // Assume N can be decomposed as N = N_x * N_y * N_z,
  // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
  // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION_SIZE-1 equals N.
  MPI_Dims_create(_mpi_pro.all_ranks, DIMENSION_SIZE, domain._grid_size);

  int node_dims[DIMENSION_SIZE];
  for (int d = 0; d < 3; d++) {
    if ((domain._grid_size[d] % _mpi_map_3d_sub_dim[d] != 0) || (domain._grid_size[d] < _mpi_map_3d_sub_dim[d])) {
      throw std::runtime_error("error: MPI map 3D sub-dimension does not match grid size in dimension " +
                               std::to_string(d) + ".");
      MPI_Abort(MPI_COMM_WORLD, -1);
    } else {
      node_dims[d] = domain._grid_size[d] / _mpi_map_3d_sub_dim[d];
    }
  }

  int period[DIMENSION_SIZE];
  // 3d with period boundary
  for (int &d : period) {
    d = 1;
  }

  comm::topology::Comm3dCart node_shared_cart_comm;
  // map node and ranks in a node into 3d topology
  comm::topology::LibComm_shared_Cart_3d_create(_mpi_pro.comm, DIMENSION_SIZE, node_dims, _mpi_map_3d_sub_dim, period,
                                                true, &node_shared_cart_comm);

  // get cartesian coordinate of current processor.
  comm::topology::LibComm_shared_Cart_3d_coords(&node_shared_cart_comm, node_dims, _mpi_map_3d_sub_dim,
                                                domain._grid_size, domain._grid_coord);

  // get rank is of neighbor mpi rank
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    comm::topology::LibComm_shared_Cart_3d_shift(domain._grid_size, d, 1, domain._grid_coord,
                                                 &(domain._rank_id_neighbours[d][DIR_LOWER]),
                                                 &(domain._rank_id_neighbours[d][DIR_HIGHER]));
  }

  comm::topology::LibComm_shared_Cart_3d_commit(_mpi_pro.comm, _p_comm, domain._grid_size, domain._grid_coord);

// #define LIBCOMM_DEBUG
#ifdef LIBCOMM_DEBUG
  int world_rank, world_size;
  MPI_Comm_rank(_mpi_pro.comm, &world_rank);
  MPI_Comm_size(_mpi_pro.comm, &world_size);

  // get process numbers and quantities within the node
  int node_id = 0, rank_in_node, node_size;
  MPI_Comm_rank(node_shared_cart_comm.node_comm, &rank_in_node);
  MPI_Comm_size(node_shared_cart_comm.node_comm, &node_size);
  if (rank_in_node == 0) {
    MPI_Comm_rank(node_shared_cart_comm.node_root_comm, &node_id);
  }
  MPI_Bcast(&node_id, 1, MPI_INT, 0, node_shared_cart_comm.node_comm);

  int new_rank;
  int new_size;
  MPI_Comm_rank(*_p_comm, &new_rank);
  MPI_Comm_size(*_p_comm, &new_size);
  MPI_Barrier(*_p_comm);
  for (int r = 0; r < world_size; r++) {
    if (new_rank == r) {
      char processor_name[MPI_MAX_PROCESSOR_NAME];
      int name_len;
      MPI_Get_processor_name(processor_name, &name_len);
      printf("[new-rank %2d | node %2d#%s] [total_size %d | inner_size %d] coord=(%d,%d,%d) "
             "neigh[x-%d x+%d y-%d y+%d z-%d z+%d] grid_size=(%d,%d,%d) inner_grid=(%d,%d,%d)\n",
             new_rank, node_id, processor_name, new_size, node_size, domain._grid_coord[0], domain._grid_coord[1],
             domain._grid_coord[2], domain._rank_id_neighbours[0][DIR_LOWER], domain._rank_id_neighbours[0][DIR_HIGHER],
             domain._rank_id_neighbours[1][DIR_LOWER], domain._rank_id_neighbours[1][DIR_HIGHER],
             domain._rank_id_neighbours[2][DIR_LOWER], domain._rank_id_neighbours[2][DIR_HIGHER], domain._grid_size[0],
             domain._grid_size[1], domain._grid_size[2], _mpi_map_3d_sub_dim[0], _mpi_map_3d_sub_dim[1],
             _mpi_map_3d_sub_dim[2]);
      fflush(stdout);
    }
    MPI_Barrier(*_p_comm);
  }
#endif
}

template <typename B, typename D> void comm::Builder<B, D>::createGlobalDomain(D &domain) {
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    // phaseSpace个单位长度(单位长度即latticeconst)
    domain._meas_global_length[d] = _phase_space[d] * _lattice_const;
    domain._meas_global_region.low[d] = 0; // lower bounding is set to 0 by default.
    domain._meas_global_region.high[d] = domain._meas_global_length[d];
  }
}

template <typename B, typename D> void comm::Builder<B, D>::buildLatticeDomain(D &domain) {
  // set lattice size of sub-box.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    domain._sub_box_lattice_size[d] = _phase_space[d] / domain._grid_size[d] +
                                      (domain._grid_coord[d] < (_phase_space[d] % domain._grid_size[d]) ? 1 : 0);
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
    domain._sub_box_lattice_region.high[d] = domain._sub_box_lattice_region.low[d] + domain._sub_box_lattice_size[d];
    // (domain._grid_coord[d] + 1) * _phase_space[d] / domain._grid_size[d];
  }

  // set lattice coordinate boundary for ghost.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    domain._ghost_ext_lattice_region.low[d] = domain._sub_box_lattice_region.low[d] - domain._lattice_size_ghost[d];
    domain._ghost_ext_lattice_region.high[d] = domain._sub_box_lattice_region.high[d] + domain._lattice_size_ghost[d];
  }

  // set local lattice coordinate.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    domain._local_ghost_ext_lattice_region.low[d] = 0;
    domain._local_ghost_ext_lattice_region.high[d] = domain._ghost_extended_lattice_size[d];
    domain._local_sub_box_lattice_region.low[d] = domain._lattice_size_ghost[d];
    domain._local_sub_box_lattice_region.high[d] = domain._lattice_size_ghost[d] + domain._sub_box_lattice_size[d];
  }
}

template <typename B, typename D> void comm::Builder<B, D>::buildMeasuredDomain(D &domain) {
  // calculate measured length in each dimension.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    // the lower and upper bounding of current sub-box.
    domain._meas_sub_box_region.low[d] =
        domain._meas_global_region.low[d] + domain._sub_box_lattice_region.low[d] * _lattice_const;
    domain._meas_sub_box_region.high[d] =
        domain._meas_global_region.low[d] + domain._sub_box_lattice_region.high[d] * _lattice_const;

    domain._meas_ghost_length[d] = domain._lattice_size_ghost[d] * _lattice_const; // ghost length fixme

    domain._meas_ghost_ext_region.low[d] = domain._meas_sub_box_region.low[d] - domain._meas_ghost_length[d];
    domain._meas_ghost_ext_region.high[d] = domain._meas_sub_box_region.high[d] + domain._meas_ghost_length[d];
  }
}
