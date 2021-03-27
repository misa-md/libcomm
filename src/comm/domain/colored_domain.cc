//
// Created by genshen on 2019-07-21.
//

#include "colored_domain.h"

void comm::ColoredDomain::splitSector(const comm::Domain &domain) {
  const _type_lattice_size sector_size_x_low = domain.sub_box_lattice_size[0] / 2;
  const _type_lattice_size sector_size_x_high = domain.sub_box_lattice_size[0] - domain.sub_box_lattice_size[0] / 2;
  const _type_lattice_size sector_size_y_low = domain.sub_box_lattice_size[1] / 2;
  const _type_lattice_size sector_size_y_high = domain.sub_box_lattice_size[1] - domain.sub_box_lattice_size[1] / 2;
  const _type_lattice_size sector_size_z_low = domain.sub_box_lattice_size[2] / 2;
  const _type_lattice_size sector_size_z_high = domain.sub_box_lattice_size[2] - domain.sub_box_lattice_size[2] / 2;

  // todo refactor: better implementation
  // set sector size in x direction for 8 sectors.
  sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0] = sector_size_x_low;
  sector_lattice_size[X_LOW | Y_HIGH | Z_LOW][0] = sector_size_x_low;
  sector_lattice_size[X_LOW | Y_LOW | Z_HIGH][0] = sector_size_x_low;
  sector_lattice_size[X_LOW | Y_HIGH | Z_HIGH][0] = sector_size_x_low;

  sector_lattice_size[X_HIGH | Y_LOW | Z_LOW][0] = sector_size_x_high;
  sector_lattice_size[X_HIGH | Y_HIGH | Z_LOW][0] = sector_size_x_high;
  sector_lattice_size[X_HIGH | Y_LOW | Z_HIGH][0] = sector_size_x_high;
  sector_lattice_size[X_HIGH | Y_HIGH | Z_HIGH][0] = sector_size_x_high;

  // set sector size in y direction for 8 sectors.
  sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1] = sector_size_y_low;
  sector_lattice_size[X_HIGH | Y_LOW | Z_LOW][1] = sector_size_y_low;
  sector_lattice_size[X_LOW | Y_LOW | Z_HIGH][1] = sector_size_y_low;
  sector_lattice_size[X_HIGH | Y_LOW | Z_HIGH][1] = sector_size_y_low;

  sector_lattice_size[X_LOW | Y_HIGH | Z_LOW][1] = sector_size_y_high;
  sector_lattice_size[X_HIGH | Y_HIGH | Z_LOW][1] = sector_size_y_high;
  sector_lattice_size[X_LOW | Y_HIGH | Z_HIGH][1] = sector_size_y_high;
  sector_lattice_size[X_HIGH | Y_HIGH | Z_HIGH][1] = sector_size_y_high;

  // set sector size in z direction for 8 sectors.
  sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2] = sector_size_z_low;
  sector_lattice_size[X_HIGH | Y_LOW | Z_LOW][2] = sector_size_z_low;
  sector_lattice_size[X_LOW | Y_HIGH | Z_LOW][2] = sector_size_z_low;
  sector_lattice_size[X_HIGH | Y_HIGH | Z_LOW][2] = sector_size_z_low;

  sector_lattice_size[X_LOW | Y_LOW | Z_HIGH][2] = sector_size_z_high;
  sector_lattice_size[X_HIGH | Y_LOW | Z_HIGH][2] = sector_size_z_high;
  sector_lattice_size[X_LOW | Y_HIGH | Z_HIGH][2] = sector_size_z_high;
  sector_lattice_size[X_HIGH | Y_HIGH | Z_HIGH][2] = sector_size_z_high;

  // set local index region of 8 sectors.
  local_sector_region[X_LOW | Y_LOW | Z_LOW] = Region<_type_lattice_coord>{
      0,
      0,
      0,
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
  };
  local_sector_region[X_HIGH | Y_LOW | Z_LOW] = Region<_type_lattice_coord>{
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      0,
      0,
      domain.sub_box_lattice_size[0],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
  };
  local_sector_region[X_LOW | Y_HIGH | Z_LOW] = Region<_type_lattice_coord>{
      0,
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      0,
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      domain.sub_box_lattice_size[1],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
  };
  local_sector_region[X_HIGH | Y_HIGH | Z_LOW] = Region<_type_lattice_coord>{
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      0,
      domain.sub_box_lattice_size[0],
      domain.sub_box_lattice_size[1],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
  };

  local_sector_region[X_LOW | Y_LOW | Z_HIGH] = Region<_type_lattice_coord>{
      0,
      0,
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      domain.sub_box_lattice_size[2],
  };
  local_sector_region[X_HIGH | Y_LOW | Z_HIGH] = Region<_type_lattice_coord>{
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0], 0,
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2], domain.sub_box_lattice_size[0],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1], domain.sub_box_lattice_size[2],
  };
  local_sector_region[X_LOW | Y_HIGH | Z_HIGH] = Region<_type_lattice_coord>{
      0,
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      domain.sub_box_lattice_size[1],
      domain.sub_box_lattice_size[2],
  };
  local_sector_region[X_HIGH | Y_HIGH | Z_HIGH] = Region<_type_lattice_coord>{
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1],
      sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2],
      domain.sub_box_lattice_size[0],
      domain.sub_box_lattice_size[1],
      domain.sub_box_lattice_size[2],
  };
  // add additional ghost offset to local_sector_region.
  // convert [0, sub_box_size/2) to [ghost_size, ghost_size + sub_box_size/2) or
  // [sub_box_size/2,  sub_box_size) to [ghost_size + sub_box_size/2,  ghost_size + sub_box_size).
  for (unsigned int i = (X_LOW | Y_LOW | Z_LOW); i <= (X_HIGH | Y_HIGH | Z_HIGH); i++) {
    local_sector_region[i].x_low = local_sector_region[i].x_low + lattice_size_ghost[0];
    local_sector_region[i].y_low = local_sector_region[i].y_low + lattice_size_ghost[1];
    local_sector_region[i].z_low = local_sector_region[i].z_low + lattice_size_ghost[2];
    local_sector_region[i].x_high = local_sector_region[i].x_high + lattice_size_ghost[0];
    local_sector_region[i].y_high = local_sector_region[i].y_high + lattice_size_ghost[1];
    local_sector_region[i].z_high = local_sector_region[i].z_high + lattice_size_ghost[2];
  }

  // extend local index region of 8 sectors with its ghost area.
  for (unsigned int i = (X_LOW | Y_LOW | Z_LOW); i <= (X_HIGH | Y_HIGH | Z_HIGH); i++) {
    local_ghost_ext_sector_region[i].x_low = local_sector_region[i].x_low - lattice_size_ghost[0];
    local_ghost_ext_sector_region[i].y_low = local_sector_region[i].y_low - lattice_size_ghost[1];
    local_ghost_ext_sector_region[i].z_low = local_sector_region[i].z_low - lattice_size_ghost[2];
    local_ghost_ext_sector_region[i].x_high = local_sector_region[i].x_high + lattice_size_ghost[0];
    local_ghost_ext_sector_region[i].y_high = local_sector_region[i].y_high + lattice_size_ghost[1];
    local_ghost_ext_sector_region[i].z_high = local_sector_region[i].z_high + lattice_size_ghost[2];
  }

  // set split coord (sector size of sector 0 plus ghost size.)
  local_split_coord[0] = sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0] + lattice_size_ghost[0];
  local_split_coord[1] = sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1] + lattice_size_ghost[1];
  local_split_coord[2] = sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2] + lattice_size_ghost[2];
}

comm::ColoredDomain::ColoredDomain(const std::array<uint64_t, DIMENSION_SIZE> _phase_space, const double _lattice_const,
                                   const double _cutoff_radius_factor)
    : Domain(_phase_space, _lattice_const, _cutoff_radius_factor) {}

comm::ColoredDomain::ColoredDomain(const comm::Domain &domain) : Domain(domain) { splitSector(domain); }

comm::ColoredDomain *comm::ColoredDomain::Builder::build() {
  ColoredDomain *p_domain = new ColoredDomain(_phase_space, _lattice_const, _cutoff_radius_factor);
  decomposition(*p_domain);
  createGlobalDomain(*p_domain);
  buildLatticeDomain(*p_domain);
  buildMeasuredDomain(*p_domain);
  p_domain->splitSector(*p_domain); // set sector size and sectors region.
  return p_domain;
}

comm::ColoredDomain *comm::ColoredDomain::Builder::localBuild(const int *_grid_size, const int *_grid_coord) {
  ColoredDomain *p_domain = new ColoredDomain(_phase_space, _lattice_const, _cutoff_radius_factor);
  for (int i = 0; i < 3; i++) {
    p_domain->_grid_size[i] = _grid_size[i];
    p_domain->_grid_coord[i] = _grid_coord[i];
  }
  createGlobalDomain(*p_domain);
  buildLatticeDomain(*p_domain);
  buildMeasuredDomain(*p_domain);
  p_domain->splitSector(*p_domain);
  return p_domain;
}
