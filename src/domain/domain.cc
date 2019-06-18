// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.

#include <algorithm>
#include <cmath>

#include "domain.h"


comm::Domain::Domain(const std::array<u_int64_t, DIMENSION_SIZE> _phase_space,
                     const double _lattice_const, const double _cutoff_radius_factor)
        : lattice_const(_lattice_const), cutoff_radius_factor(_cutoff_radius_factor),
          cut_lattice(static_cast<int>(ceil(_cutoff_radius_factor))), phase_space(_phase_space) {}
//      todo _grid_size(0),
//      todo _meas_global_length(0.0),

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
