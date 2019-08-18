//
// Created by genshen on 2019-08-18.
//

#include "sector_forwarding_region.h"

// table of sending regions of each sectors in each dimension.
comm::type_region_array comm::fwCommSectorLocalRegion(
        const comm::ColoredDomain *p_domain, const unsigned int sector_id,
        const unsigned int dimension) {
    type_region_array regions;
    switch (sector_id) {
        case X_LOW | Y_LOW | Z_LOW: { // sector 0
            const _type_lattice_size &x_c = p_domain->local_sector_region[sector_id].z_high;
            const _type_lattice_size &y_c = p_domain->local_sector_region[sector_id].y_high;
            const _type_lattice_size &z_c = p_domain->local_sector_region[sector_id].z_high;
            switch (dimension) {
                case DIM_X:
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_high -
                            p_domain->lattice_size_ghost[DIM_X], // x_e-g
                            p_domain->local_sub_box_lattice_region.y_low, // y_b
                            p_domain->local_sub_box_lattice_region.z_low, // z_b
                            p_domain->local_sub_box_lattice_region.x_high, // x_e
                            y_c + p_domain->lattice_size_ghost[DIM_Y], // y_c+g
                            z_c + p_domain->lattice_size_ghost[DIM_Z]// z_c+g
                    ));
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_high -
                            p_domain->lattice_size_ghost[DIM_X], // x_e-g
                            p_domain->local_sub_box_lattice_region.y_low, // y_b
                            p_domain->local_sub_box_lattice_region.z_high
                            - p_domain->lattice_size_ghost[DIM_Z], // z_e-g
                            p_domain->local_sub_box_lattice_region.x_high, // x_e
                            y_c + p_domain->lattice_size_ghost[DIM_Y], // y_c+g
                            p_domain->local_sub_box_lattice_region.z_high // z_e
                    ));
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_high -
                            p_domain->lattice_size_ghost[DIM_X], // x_e-g
                            p_domain->local_sub_box_lattice_region.y_high
                            - p_domain->lattice_size_ghost[DIM_Y], // y_e-g
                            p_domain->local_sub_box_lattice_region.z_low, // z_b
                            p_domain->local_sub_box_lattice_region.x_high, // x_e
                            p_domain->local_sub_box_lattice_region.y_high, // y_e
                            z_c + p_domain->lattice_size_ghost[DIM_Z] // z_c+g
                    ));
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_high
                            - p_domain->lattice_size_ghost[DIM_X], // x_e-g
                            p_domain->local_sub_box_lattice_region.y_high
                            - p_domain->lattice_size_ghost[DIM_Y], // y_e-g
                            p_domain->local_sub_box_lattice_region.z_high
                            - p_domain->lattice_size_ghost[DIM_Z], // z_e-g
                            p_domain->local_sub_box_lattice_region.x_high, // x_e
                            p_domain->local_sub_box_lattice_region.y_high, // y_e
                            p_domain->local_sub_box_lattice_region.z_high  // z_e
                    ));
                    break;
                case DIM_Y:
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_low
                            - p_domain->lattice_size_ghost[DIM_X], // x_b-g
                            p_domain->local_sub_box_lattice_region.y_high
                            - p_domain->lattice_size_ghost[DIM_Y], // y_e-g
                            p_domain->local_sub_box_lattice_region.z_low, // z_b
                            x_c + p_domain->lattice_size_ghost[DIM_X], // x_c+g
                            p_domain->local_sub_box_lattice_region.y_high, // y_e
                            z_c + p_domain->lattice_size_ghost[DIM_Z] // z_c+g
                    ));
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_low
                            - p_domain->lattice_size_ghost[DIM_X], // x_b-g
                            p_domain->local_sub_box_lattice_region.y_high
                            - p_domain->lattice_size_ghost[DIM_Y], // y_e-g
                            p_domain->local_sub_box_lattice_region.z_high
                            - p_domain->lattice_size_ghost[DIM_Z], // z_e-g
                            x_c + p_domain->lattice_size_ghost[DIM_X], // x_c+g
                            p_domain->local_sub_box_lattice_region.y_high, // y_e
                            p_domain->local_sub_box_lattice_region.z_high  // z_e
                    ));
                    break;
                case DIM_Z:
                    regions.push_back(comm::Region<comm::_type_lattice_size>(
                            p_domain->local_sub_box_lattice_region.x_low
                            - p_domain->lattice_size_ghost[DIM_X], // x_b-g
                            p_domain->local_sub_box_lattice_region.x_low
                            - p_domain->lattice_size_ghost[DIM_Y], // y_b-g
                            p_domain->local_sub_box_lattice_region.z_high
                            - p_domain->lattice_size_ghost[DIM_Z], // z_e-g
                            x_c + p_domain->lattice_size_ghost[DIM_X], // x_c+g
                            y_c + p_domain->lattice_size_ghost[DIM_Y], // y_c+g
                            p_domain->local_sub_box_lattice_region.z_high  // z_e
                    ));
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_HIGH | Y_LOW | Z_LOW: { // sector 1 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_LOW | Y_HIGH | Z_LOW: { // sector 2 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_HIGH | Y_HIGH | Z_LOW: { // sector 3 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_LOW | Y_LOW | Z_HIGH: { // sector 4 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_HIGH | Y_LOW | Z_HIGH: {  // sector 5 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_LOW | Y_HIGH | Z_HIGH: { // sector 6 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        case X_HIGH | Y_HIGH | Z_HIGH: { // sector 7 todo
            switch (dimension) {
                case DIM_X:
                    break;
                case DIM_Y:
                    break;
                case DIM_Z:
                    break;
                default:
                    assert(false); // this branch is not allow.
            }
        }
            break;
        default:
            assert(false); // this branch is not allow.
    }
    return regions;
}
