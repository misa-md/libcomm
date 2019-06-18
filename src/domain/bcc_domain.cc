//
// Created by genshen on 2019-04-19.
//

#include "bcc_domain.h"


comm::BccDomain::BccDomain(const std::array<u_int64_t, 3> _phase_space, const double _lattice_const,
                           const double _cutoff_radius_factor)
        : Domain(_phase_space, _lattice_const, _cutoff_radius_factor) {}

comm::BccDomain::BccDomain(const comm::Domain &domain) : Domain(domain) {
    rescale(domain);
}

void comm::BccDomain::rescale(const comm::Domain &domain) {
    _dbx_lattice_sub_box_size[0] = 2 * domain.lattice_size_sub_box[0];
    _dbx_lattice_sub_box_size[1] = domain.lattice_size_sub_box[1];
    _dbx_lattice_sub_box_size[2] = domain.lattice_size_sub_box[2];

    _dbx_lattice_size_ghost[0] = 2 * domain.lattice_size_ghost[0];
    _dbx_lattice_size_ghost[1] = domain.lattice_size_ghost[1];
    _dbx_lattice_size_ghost[2] = domain.lattice_size_ghost[2];

    _dbx_lattice_size_ghost_extended[0] = 2 * domain.lattice_size_ghost_extended[0];
    _dbx_lattice_size_ghost_extended[1] = domain.lattice_size_ghost_extended[1];
    _dbx_lattice_size_ghost_extended[2] = domain.lattice_size_ghost_extended[2];

    _dbx_lattice_coord_sub_box_region.x_low = 2 * domain.lattice_coord_sub_box_region.x_low;
    _dbx_lattice_coord_sub_box_region.y_low = domain.lattice_coord_sub_box_region.y_low;
    _dbx_lattice_coord_sub_box_region.z_low = domain.lattice_coord_sub_box_region.z_low;
    _dbx_lattice_coord_sub_box_region.x_high = 2 * domain.lattice_coord_sub_box_region.x_high;
    _dbx_lattice_coord_sub_box_region.y_high = domain.lattice_coord_sub_box_region.y_high;
    _dbx_lattice_coord_sub_box_region.z_high = domain.lattice_coord_sub_box_region.z_high;

    _dbx_lattice_coord_ghost_region.x_low = 2 * domain.lattice_coord_ghost_region.x_low;
    _dbx_lattice_coord_ghost_region.y_low = domain.lattice_coord_ghost_region.y_low;
    _dbx_lattice_coord_ghost_region.z_low = domain.lattice_coord_ghost_region.z_low;
    _dbx_lattice_coord_ghost_region.x_high = 2 * domain.lattice_coord_ghost_region.x_high;
    _dbx_lattice_coord_ghost_region.y_high = domain.lattice_coord_ghost_region.y_high;
    _dbx_lattice_coord_ghost_region.z_high = domain.lattice_coord_ghost_region.z_high;

    _dbx_local_ghost_lattice_coord_region.x_low = 2 * domain.local_ghost_lattice_coord_region.x_low;
    _dbx_local_ghost_lattice_coord_region.y_low = domain.local_ghost_lattice_coord_region.y_low;
    _dbx_local_ghost_lattice_coord_region.z_low = domain.local_ghost_lattice_coord_region.z_low;
    _dbx_local_ghost_lattice_coord_region.x_high = 2 * domain.local_ghost_lattice_coord_region.x_high;
    _dbx_local_ghost_lattice_coord_region.y_high = domain.local_ghost_lattice_coord_region.y_high;
    _dbx_local_ghost_lattice_coord_region.z_high = domain.local_ghost_lattice_coord_region.z_high;

    _dbx_local_sub_box_lattice_coord_region.x_low = 2 * domain.local_sub_box_lattice_coord_region.x_low;
    _dbx_local_sub_box_lattice_coord_region.y_low = domain.local_sub_box_lattice_coord_region.y_low;
    _dbx_local_sub_box_lattice_coord_region.z_low = domain.local_sub_box_lattice_coord_region.z_low;
    _dbx_local_sub_box_lattice_coord_region.x_high = 2 * domain.local_sub_box_lattice_coord_region.x_high;
    _dbx_local_sub_box_lattice_coord_region.y_high = domain.local_sub_box_lattice_coord_region.y_high;
    _dbx_local_sub_box_lattice_coord_region.z_high = domain.local_sub_box_lattice_coord_region.z_high;
}

comm::BccDomain *comm::BccDomain::Builder::build() {
    BccDomain *p_domain = new BccDomain(_phase_space, _lattice_const, _cutoff_radius_factor);
    decomposition(*p_domain);
    createGlobalDomain(*p_domain);
    buildLatticeDomain(*p_domain);
    buildMeasuredDomain(*p_domain);
    p_domain->rescale(*p_domain);
    return p_domain;
}

comm::BccDomain *comm::BccDomain::Builder::localBuild(const int *_grid_size, const int *_grid_coord) {
    BccDomain *p_domain = new BccDomain(_phase_space, _lattice_const, _cutoff_radius_factor);
    for (int i = 0; i < 3; i++) {
        p_domain->_grid_size[i] = _grid_size[i];
        p_domain->_grid_coord_sub_box[i] = _grid_coord[i];
    }
    createGlobalDomain(*p_domain);
    buildLatticeDomain(*p_domain);
    buildMeasuredDomain(*p_domain);
    p_domain->rescale(*p_domain);
    return p_domain;
}
