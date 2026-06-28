// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.

#include <algorithm>
#include <cmath>

#include "domain.h"

void comm::MeasuredDomain::rescale_measured(const std::array<double, DIMENSION_SIZE> scale_factors) {
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    this->_meas_global_length[d] *= scale_factors[d];
    this->_meas_global_region.low[d] *= scale_factors[d]; // lower bounding is set to 0 by default.
    this->_meas_global_region.high[d] *= scale_factors[d];
  }

  // calculate measured length in each dimension.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    // the lower and upper bounding of current sub-box.
    this->_meas_sub_box_region.low[d] *= scale_factors[d];
    this->_meas_sub_box_region.high[d] *= scale_factors[d];

    // ghost size is not scaled
    this->_meas_ghost_ext_region.low[d] = this->_meas_sub_box_region.low[d] - this->_meas_ghost_length[d];
    this->_meas_ghost_ext_region.high[d] = this->_meas_sub_box_region.high[d] + this->_meas_ghost_length[d];
  }
}

comm::LatticeDomain::LatticeDomain(const std::array<uint64_t, DIMENSION_SIZE> _phase_space,
                                   const std::array<double, DIMENSION_SIZE> _lattice_const, const double _cutoff_radius)
    : lattice_const(_lattice_const), cutoff_radius(_cutoff_radius),
      cut_lattice(static_cast<int>(ceil(cutoff_radius_factor()))), phase_space(_phase_space) {}
//  todo _grid_size(0),
//  todo _meas_global_length(0.0),

void comm::LatticeDomain::rescale(const double scale_factor) {
  const std::array<double, DIMENSION_SIZE> scale_factor_3d = {scale_factor, scale_factor, scale_factor};
  rescale(scale_factor_3d);
}

void comm::LatticeDomain::rescale(const std::array<double, DIMENSION_SIZE> scale_factors) {
  // rescale lattice const
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    this->lattice_const[d] *= scale_factors[d];
  }
  // rescale global box
  this->rescale_measured(scale_factors);
}

comm::LatticeDomain *comm::LatticeDomain::Builder::build() {
  LatticeDomain *p_domain = new LatticeDomain(_phase_space, _lattice_const, _cutoff_radius);
  decomposition(*p_domain);
  createGlobalDomain(*p_domain);
  buildLatticeDomain(*p_domain);
  buildMeasuredDomain(*p_domain);
  return p_domain;
}

comm::LatticeDomain *comm::LatticeDomain::Builder::localBuild(const int _grid_size[DIMENSION_SIZE],
                                                const int _grid_coord[DIMENSION_SIZE]) {
  LatticeDomain *p_domain = new LatticeDomain(_phase_space, _lattice_const, _cutoff_radius);
  for (int i = 0; i < 3; i++) {
    p_domain->_grid_size[i] = _grid_size[i];
    p_domain->_grid_coord[i] = _grid_coord[i];
  }
  createGlobalDomain(*p_domain);
  buildLatticeDomain(*p_domain);
  buildMeasuredDomain(*p_domain);
  return p_domain;
}
