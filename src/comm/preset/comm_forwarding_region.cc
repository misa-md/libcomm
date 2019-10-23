//
// Created by genshen on 2019-04-16.
//

#include "comm/domain/bcc_domain.h"
#include "comm/domain/region.hpp"
#include "comm/types_define.h"
#include "comm_forwarding_region.h"

comm::Region<comm::_type_lattice_size> comm::fwCommLocalRegion(
        const comm::BccDomain *p_domain, const unsigned int dimension, const unsigned int direction) {
    switch (dimension << 2 | direction) {
        case DIM_X << 2 | DIR_LOWER: { // x dimension, lower direction
            _type_lattice_size xstart = p_domain->dbx_lattice_size_ghost[0];
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop =
                    xstart + p_domain->dbx_lattice_size_ghost[0]; // note: this is ghost lattice size.
            _type_lattice_size ystop = ystart + p_domain->dbx_sub_box_lattice_size[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_sub_box_lattice_size[2];

            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_X << 2 | DIR_HIGHER: { // x dimension, higher direction
            _type_lattice_size xstart = p_domain->dbx_lattice_size_ghost[0] +
                                        p_domain->dbx_sub_box_lattice_size[0] -
                                        p_domain->dbx_lattice_size_ghost[0]; // fixme note: 2 was double for BCC lattice.
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_lattice_size_ghost[0] + p_domain->dbx_sub_box_lattice_size[0];
            _type_lattice_size ystop = ystart + p_domain->dbx_sub_box_lattice_size[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_sub_box_lattice_size[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Y << 2 | DIR_LOWER: { // y dimension, lower direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_ghost_extended_lattice_size[0];
            _type_lattice_size ystop = ystart + p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_sub_box_lattice_size[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Y << 2 | DIR_HIGHER: { // y dimension, higher direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1] +
                                        p_domain->dbx_sub_box_lattice_size[1] - p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_ghost_extended_lattice_size[0];
            _type_lattice_size ystop = p_domain->dbx_lattice_size_ghost[1] + p_domain->dbx_sub_box_lattice_size[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_sub_box_lattice_size[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Z << 2 | DIR_LOWER: { // z dimension, lower direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = 0;
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_ghost_extended_lattice_size[0];
            _type_lattice_size ystop = p_domain->dbx_ghost_extended_lattice_size[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_lattice_size_ghost[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Z << 2 | DIR_HIGHER: { // z dimension, higher direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = 0;
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_sub_box_lattice_size[2]
                                        - p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_ghost_extended_lattice_size[0];
            _type_lattice_size ystop = p_domain->dbx_ghost_extended_lattice_size[1];
            _type_lattice_size zstop = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_sub_box_lattice_size[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        default:
            // this case is not allowed.
            assert(false);
    }
}

comm::Region<double> comm::fwCommLocalMeaRegion(
        const comm::Domain *p_domain, const unsigned int dimension, const unsigned int direction) {
    switch (dimension << 2 | direction) {
        case DIM_X << 2 | DIR_LOWER: { // x dimension, lower direction
            return Region<double>(
                    p_domain->meas_sub_box_region.x_low,
                    p_domain->meas_sub_box_region.x_low + p_domain->meas_ghost_length[dimension],
                    p_domain->meas_sub_box_region.y_low,
                    p_domain->meas_sub_box_region.y_high,
                    p_domain->meas_sub_box_region.z_low,
                    p_domain->meas_sub_box_region.z_high
            );
        }
        case DIM_X << 2 | DIR_HIGHER: { // x dimension, higher direction
            return Region<double>(
                    p_domain->meas_sub_box_region.x_high - p_domain->meas_ghost_length[dimension],
                    p_domain->meas_sub_box_region.x_high,
                    p_domain->meas_sub_box_region.y_low,
                    p_domain->meas_sub_box_region.y_high,
                    p_domain->meas_sub_box_region.z_low,
                    p_domain->meas_sub_box_region.z_high
            );
        }
        case DIM_Y << 2 | DIR_LOWER: { // y dimension, lower direction
            return Region<double>(
                    p_domain->meas_ghost_ext_region.x_low,
                    p_domain->meas_sub_box_region.y_low,
                    p_domain->meas_sub_box_region.z_low,
                    p_domain->meas_ghost_ext_region.x_high,
                    p_domain->meas_sub_box_region.y_low + p_domain->meas_ghost_length[dimension],
                    p_domain->meas_sub_box_region.z_high
            );
        }
        case DIM_Y << 2 | DIR_HIGHER: { // y dimension, higher direction
            return Region<double>(
                    p_domain->meas_ghost_ext_region.x_low,
                    p_domain->meas_sub_box_region.y_high - p_domain->meas_ghost_length[dimension],
                    p_domain->meas_sub_box_region.z_low,
                    p_domain->meas_ghost_ext_region.x_high,
                    p_domain->meas_sub_box_region.y_high,
                    p_domain->meas_sub_box_region.z_high
            );
        }
        case DIM_Z << 2 | DIR_LOWER: { // z dimension, lower direction
            return Region<double>(
                    p_domain->meas_ghost_ext_region.x_low,
                    p_domain->meas_ghost_ext_region.y_low,
                    p_domain->meas_sub_box_region.z_low,
                    p_domain->meas_ghost_ext_region.x_high,
                    p_domain->meas_ghost_ext_region.y_high,
                    p_domain->meas_sub_box_region.z_low + p_domain->meas_ghost_length[dimension]
            );
        }
        case DIM_Z << 2 | DIR_HIGHER: { // z dimension, higher direction
            return Region<double>(
                    p_domain->meas_ghost_ext_region.x_low,
                    p_domain->meas_ghost_ext_region.y_low,
                    p_domain->meas_sub_box_region.z_high - p_domain->meas_ghost_length[dimension],
                    p_domain->meas_ghost_ext_region.x_high,
                    p_domain->meas_ghost_ext_region.y_high,
                    p_domain->meas_sub_box_region.z_high
            );
        }
        default:
            // this case is not allowed.
            assert(false);
    }

}