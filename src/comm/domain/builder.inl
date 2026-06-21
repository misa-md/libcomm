//
// Created by genshen on 2019-04-19.
//

#include "builder.h"
#include "comm/topology/cart_3d_node_affinity.h"

#include <algorithm>

template <typename B, typename D> B &comm::Builder<B, D>::setPhaseSpace(const int64_t phaseSpace[DIMENSION_SIZE]) {
  for (int i = 0; i < DIMENSION_SIZE; i++) {
    _phase_space[i] = phaseSpace[i]; // todo type not match
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D> B &comm::Builder<B, D>::setLatticeConst(const double latticeConst) {
  _lattice_const[0] = latticeConst;
  _lattice_const[1] = latticeConst;
  _lattice_const[2] = latticeConst;
  return *static_cast<B *>(this);
}

template <typename B, typename D>
B &comm::Builder<B, D>::setLatticeConst(const std::array<double, DIMENSION_SIZE> lat_const) {
  _lattice_const[0] = lat_const[0];
  _lattice_const[1] = lat_const[1];
  _lattice_const[2] = lat_const[2];
  return *static_cast<B *>(this);
}

template <typename B, typename D> B &comm::Builder<B, D>::setGhostSize(const unsigned int ghost_size) {
  return setGhostSize({ghost_size, ghost_size, ghost_size});
}

template <typename B, typename D>
B &comm::Builder<B, D>::setGhostSize(const std::array<unsigned int, DIMENSION_SIZE> ghost_size) {
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    _ghost_lat_size[d] = ghost_size[d];
    // also set the measured ghost region length. Note: it can be overwritten by `ghost_measured_length`
    _ghost_meas_length[d] = _lattice_const[d] * ghost_size[d];
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D>
B &comm::Builder<B, D>::setGhostMeasLength(const std::array<double, DIMENSION_SIZE> ghost_measured_length) {
  // overwrite
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    _ghost_meas_length[d] = _lattice_const[d] * ghost_measured_length[d];
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D>
B &comm::Builder<B, D>::setCutoffRadius(const double cutoff_radius_factor, const double default_lat_const) {
  _cutoff_radius = cutoff_radius_factor * default_lat_const;
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    if (_ghost_lat_size[d] == 0) {
      // if ghost size not set, we set it as cut_lattice
      _ghost_lat_size[d] = static_cast<int>(ceil(cutoff_radius_factor));
    }
  }
  return *static_cast<B *>(this);
}

template <typename B, typename D>
B &comm::Builder<B, D>::setCutoffRadius_v2(const double cutoff_radius, const double default_lat_const) {
  _cutoff_radius = cutoff_radius;
  const double cutoff_radius_factor = cutoff_radius / default_lat_const;

  for (int d = 0; d < DIMENSION_SIZE; d++) {
    if (_ghost_lat_size[d] == 0) {
      // if ghost size not set, we set it as cut_lattice
      _ghost_lat_size[d] = static_cast<int>(ceil(cutoff_radius_factor));
    }
    if (_ghost_meas_length[d] == 0.0) {
      _ghost_meas_length[d] = _cutoff_radius;
    }
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
  domain.decomposition(_mpi_pro, _p_comm, _mpi_map_3d_sub_dim);
}

template <typename B, typename D> void comm::Builder<B, D>::createGlobalDomain(D &domain) {
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    // phaseSpace个单位长度(单位长度即latticeconst)
    domain._meas_global_length[d] = _phase_space[d] * _lattice_const[d];
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
    domain._lattice_size_ghost[d] = _ghost_lat_size[d];
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
        domain._meas_global_region.low[d] + domain._sub_box_lattice_region.low[d] * _lattice_const[d];
    domain._meas_sub_box_region.high[d] =
        domain._meas_global_region.low[d] + domain._sub_box_lattice_region.high[d] * _lattice_const[d];

    domain._meas_ghost_length[d] = _ghost_meas_length[d];

    domain._meas_ghost_ext_region.low[d] = domain._meas_sub_box_region.low[d] - domain._meas_ghost_length[d];
    domain._meas_ghost_ext_region.high[d] = domain._meas_sub_box_region.high[d] + domain._meas_ghost_length[d];
  }
}
