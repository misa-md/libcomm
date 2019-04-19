// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.

#include <algorithm>
#include <cmath>

#include "domain.h"


comm::Domain::Domain(const std::array<u_int64_t, DIMENSION_SIZE> _phase_space,
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

comm::Domain *comm::Domain::Builder::build() {
    Domain *p_domain = new Domain(_phase_space, _lattice_const, _cutoff_radius_factor);
    decomposition(*p_domain);
    createGlobalDomain(*p_domain);
    buildLatticeDomain(*p_domain);
    buildMeasuredDomain(*p_domain);
    return p_domain;
}

comm::Domain *comm::Domain::Builder::localBuild(const int _grid_size[DIMENSION_SIZE],
                                                const int _grid_coord[DIMENSION_SIZE]) {
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
