//
// Created by genshen on 2019-04-16.
//

#include "domain/region.hpp"
#include "domain/domain.h"
#include "types_define.h"
#include "comm_forwarding_region.h"

comm::Region<comm::_type_lattice_size> comm::fwCommLocalRegion(
        const comm::Domain *p_domain, const int dimension, const int direction) {
    switch (dimension << 2 | direction) {
        case DIM_X << 2 | DIR_LOWER: { // x dimension, lower direction
            _type_lattice_size xstart = p_domain->dbx_lattice_size_ghost[0];
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop =
                    xstart + p_domain->dbx_lattice_size_ghost[0]; // note: this is ghost lattice size.
            _type_lattice_size ystop = ystart + p_domain->dbx_lattice_size_sub_box[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];

            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_X << 2 | DIR_HIGHER: { // x dimension, higher direction
            _type_lattice_size xstart = p_domain->dbx_lattice_size_ghost[0] +
                                        p_domain->dbx_lattice_size_sub_box[0] -
                                        ((p_domain->cut_lattice) * 2); // fixme note: 2 was double for BCC lattice.
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_lattice_size_ghost[0] + p_domain->dbx_lattice_size_sub_box[0];
            _type_lattice_size ystop = ystart + p_domain->dbx_lattice_size_sub_box[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Y << 2 | DIR_LOWER: { // y dimension, lower direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1];
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_lattice_size_ghost_extended[0];
            _type_lattice_size ystop = ystart + p_domain->cut_lattice;
            _type_lattice_size zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Y << 2 | DIR_HIGHER: { // y dimension, higher direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = p_domain->dbx_lattice_size_ghost[1] +
                                        p_domain->dbx_lattice_size_sub_box[1] - p_domain->cut_lattice;
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_lattice_size_ghost_extended[0];
            _type_lattice_size ystop = p_domain->dbx_lattice_size_ghost[1] + p_domain->dbx_lattice_size_sub_box[1];
            _type_lattice_size zstop = zstart + p_domain->dbx_lattice_size_sub_box[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Z << 2 | DIR_LOWER: { // z dimension, lower direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = 0;
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2];
            _type_lattice_size xstop = p_domain->dbx_lattice_size_ghost_extended[0];
            _type_lattice_size ystop = p_domain->dbx_lattice_size_ghost_extended[1];
            _type_lattice_size zstop = zstart + p_domain->cut_lattice;
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        case DIM_Z << 2 | DIR_HIGHER: { // z dimension, higher direction
            _type_lattice_size xstart = 0;
            _type_lattice_size ystart = 0;
            _type_lattice_size zstart = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_lattice_size_sub_box[2]
                                        - (p_domain->cut_lattice);
            _type_lattice_size xstop = p_domain->dbx_lattice_size_ghost_extended[0];
            _type_lattice_size ystop = p_domain->dbx_lattice_size_ghost_extended[1];
            _type_lattice_size zstop = p_domain->dbx_lattice_size_ghost[2] + p_domain->dbx_lattice_size_sub_box[2];
            return Region<_type_lattice_size>(xstart, ystart, zstart, xstop, ystop, zstop);
        }
        default:
            // this case is not allowed.
            assert(false);
    }
}