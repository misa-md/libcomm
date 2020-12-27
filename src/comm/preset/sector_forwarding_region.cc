//
// Created by genshen on 2019-08-18.
//

#include <cassert>
#include <stdexcept>
#include "sector_forwarding_region.h"

// table of sending regions of each sectors in each dimension.
comm::type_region_array comm::fwCommSectorSendRegion(const _type_sector_id sector_id, const unsigned int dim,
                                                     const _type_lattice_size ghost_size[DIMENSION_SIZE],
                                                     const _type_lattice_coord split_coord[DIMENSION_SIZE],
                                                     const Region<comm::_type_lattice_coord> local_box_region) {
    if (sector_id < 0 || sector_id >= 8 || dim < 0 || dim >= 3) {
        throw std::invalid_argument("error sector id or dimension id.");
    }
    switch (dim) {
        case DIM_X: // x dimension
            return regionsDimXSend(sector_id, ghost_size,
                                   local_box_region.x_low, local_box_region.y_low, local_box_region.z_low,
                                   split_coord[0], split_coord[1], split_coord[2],
                                   local_box_region.x_high, local_box_region.y_high, local_box_region.z_high);
        case DIM_Y: // y dimension
            return regionsDimYSend(sector_id, ghost_size,
                                   local_box_region.x_low, local_box_region.y_low, local_box_region.z_low,
                                   split_coord[0], split_coord[1], split_coord[2],
                                   local_box_region.x_high, local_box_region.y_high, local_box_region.z_high);
        case DIM_Z: // z dimension
            return regionsDimZSend(sector_id, ghost_size,
                                   local_box_region.x_low, local_box_region.y_low, local_box_region.z_low,
                                   split_coord[0], split_coord[1], split_coord[2],
                                   local_box_region.x_high, local_box_region.y_high, local_box_region.z_high);
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::fwCommSectorRecvRegion(const _type_sector_id sector_id, const unsigned int dim,
                                                     const _type_lattice_size ghost_size[DIMENSION_SIZE],
                                                     const _type_lattice_coord split_coord[DIMENSION_SIZE],
                                                     const Region<comm::_type_lattice_coord> local_box_region) {
    if (sector_id < 0 || sector_id >= 8 || dim < 0 || dim >= 3) {
        throw std::invalid_argument("error sector id or dimension id.");
    }
    switch (dim) {
        case DIM_X: // x dimension
            return regionsDimXRecv(sector_id, ghost_size,
                                   local_box_region.x_low, local_box_region.y_low, local_box_region.z_low,
                                   split_coord[0], split_coord[1], split_coord[2],
                                   local_box_region.x_high, local_box_region.y_high, local_box_region.z_high);
        case DIM_Y: // y dimension
            return regionsDimYRecv(sector_id, ghost_size,
                                   local_box_region.x_low, local_box_region.y_low, local_box_region.z_low,
                                   split_coord[0], split_coord[1], split_coord[2],
                                   local_box_region.x_high, local_box_region.y_high, local_box_region.z_high);
        case DIM_Z: // z dimension
            return regionsDimZRecv(sector_id, ghost_size,
                                   local_box_region.x_low, local_box_region.y_low, local_box_region.z_low,
                                   split_coord[0], split_coord[1], split_coord[2],
                                   local_box_region.x_high, local_box_region.y_high, local_box_region.z_high);
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::regionsDimXSend(const _type_sector_id sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                              const _type_lattice_coord xb, const _type_lattice_coord yb,
                                              const _type_lattice_coord zb, const _type_lattice_coord xc,
                                              const _type_lattice_coord yc, const _type_lattice_coord zc,
                                              const _type_lattice_coord xe, const _type_lattice_coord ye,
                                              const _type_lattice_coord ze) {
    switch (sector_id) {
        case 0: //sector 0
            return type_region_array{
                    {xe - g[0], yb,        zb,        xe, yc + g[1], zc + g[2]},
                    {xe - g[0], yb,        ze - g[2], xe, yc + g[1], ze},
                    {xe - g[0], ye - g[1], zb,        xe, ye,        zc + g[2]},
                    {xe - g[0], ye - g[1], ze - g[2], xe, ye,        ze},
            };
        case 1: //sector 1
            return type_region_array{
                    {xb, yb,        zb,        xb + g[0], yc + g[1], zc + g[2]},
                    {xb, yb,        ze - g[2], xb + g[0], yc + g[1], ze},
                    {xb, ye - g[1], zb,        xb + g[0], ye,        zc + g[2]},
                    {xb, ye - g[1], ze - g[2], xb + g[0], ye,        ze},
            };
        case 2: //sector 2
            return type_region_array{
                    {xe - g[0], yb,        zb,        xe, yb + g[1], zc + g[2]},
                    {xe - g[0], yb,        ze - g[2], xe, yb + g[1], ze},
                    {xe - g[0], yc - g[1], zb,        xe, ye,        zc + g[2]},
                    {xe - g[0], yc - g[1], ze - g[2], xe, ye,        ze},
            };
        case 3: //sector 3
            return type_region_array{
                    {xb, yb,        zb,        xb + g[0], yb + g[1], zc + g[2]},
                    {xb, yb,        ze - g[2], xb + g[0], yb + g[1], ze},
                    {xb, yc - g[1], zb,        xb + g[0], ye,        zc + g[2]},
                    {xb, yc - g[1], ze - g[2], xb + g[0], ye,        ze},
            };
        case 4: //sector 4
            return type_region_array{
                    {xe - g[0], yb,        zb,        xe, yc + g[1], zb + g[2]},
                    {xe - g[0], yb,        zc - g[2], xe, yc + g[1], ze},
                    {xe - g[0], ye - g[1], zb,        xe, ye,        zb + g[2]},
                    {xe - g[0], ye - g[1], zc - g[2], xe, ye,        ze},
            };
        case 5: //sector 5
            return type_region_array{
                    {xb, yb,        zb,        xb + g[0], yc + g[1], zb + g[2]},
                    {xb, yb,        zc - g[2], xb + g[0], yc + g[1], ze},
                    {xb, ye - g[1], zb,        xb + g[0], ye,        zb + g[2]},
                    {xb, ye - g[1], zc - g[2], xb + g[0], ye,        ze},
            };
        case 6: //sector 6
            return type_region_array{
                    {xe - g[0], yb,        zb,        xe, yb + g[1], zb + g[2]},
                    {xe - g[0], yb,        zc - g[2], xe, yb + g[1], ze},
                    {xe - g[0], yc - g[1], zb,        xe, ye,        zb + g[2]},
                    {xe - g[0], yc - g[1], zc - g[2], xe, ye,        ze},
            };
        case 7: //sector 7
            return type_region_array{
                    {xb, yb,        zb,        xb + g[0], yb + g[1], zb + g[2]},
                    {xb, yb,        zc - g[2], xb + g[0], yb + g[1], ze},
                    {xb, yc - g[1], zb,        xb + g[0], ye,        zb + g[2]},
                    {xb, yc - g[1], zc - g[2], xb + g[0], ye,        ze},
            };
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::regionsDimYSend(const _type_sector_id sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                              const _type_lattice_coord xb, const _type_lattice_coord yb,
                                              const _type_lattice_coord zb, const _type_lattice_coord xc,
                                              const _type_lattice_coord yc, const _type_lattice_coord zc,
                                              const _type_lattice_coord xe, const _type_lattice_coord ye,
                                              const _type_lattice_coord ze) {

    switch (sector_id) {
        case 0: //sector 0
            return type_region_array{
                    {xb - g[0], ye - g[1], zb,        xc + g[0], ye, zc + g[2]},
                    {xb - g[0], ye - g[1], ze - g[2], xc + g[0], ye, ze},
            };
        case 1: //sector 1
            return type_region_array{
                    {xc - g[0], ye - g[1], zb,        xe + g[0], ye, zc + g[2]},
                    {xc - g[0], ye - g[1], ze - g[2], xe + g[0], ye, ze},
            };
        case 2: //sector 2
            return type_region_array{
                    {xb - g[0], yb, zb,        xc + g[0], yb + g[1], zc + g[2]},
                    {xb - g[0], yb, ze - g[2], xc + g[0], yb + g[1], ze},
            };
        case 3: //sector 3
            return type_region_array{
                    {xc - g[0], yb, zb,        xe + g[0], yb + g[1], zc + g[2]},
                    {xc - g[0], yb, ze - g[2], xe + g[0], yb + g[1], ze},
            };
        case 4: //sector 4
            return type_region_array{
                    {xb - g[0], ye - g[1], zb,        xc + g[0], ye, zb + g[2]},
                    {xb - g[0], ye - g[1], zc - g[2], xc + g[0], ye, ze},
            };
        case 5: //sector 5
            return type_region_array{
                    {xc - g[0], ye - g[1], zb,        xe + g[0], ye, zb + g[2]},
                    {xc - g[0], ye - g[1], zc - g[2], xe + g[0], ye, ze},
            };
        case 6: //sector 6
            return type_region_array{
                    {xb - g[0], yb, zb,        xc + g[0], yb + g[1], zb + g[2]},
                    {xb - g[0], yb, zc - g[2], xc + g[0], yb + g[1], ze},
            };
        case 7: //sector 7
            return type_region_array{
                    {xc - g[0], yb, zb,        xe + g[0], yb + g[1], zb + g[2]},
                    {xc - g[0], yb, zc - g[2], xe + g[0], yb + g[1], ze},
            };
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::regionsDimZSend(const _type_sector_id sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                              const _type_lattice_coord xb, const _type_lattice_coord yb,
                                              const _type_lattice_coord zb, const _type_lattice_coord xc,
                                              const _type_lattice_coord yc, const _type_lattice_coord zc,
                                              const _type_lattice_coord xe, const _type_lattice_coord ye,
                                              const _type_lattice_coord ze) {
    switch (sector_id) {
        case 0: //sector 0
            return type_region_array{
                    {xb - g[0], yb - g[1], ze - g[2], xc + g[0], yc + g[1], ze},
            };
        case 1: //sector 1
            return type_region_array{
                    {xc - g[0], yb - g[1], ze - g[2], xe + g[0], yc + g[1], ze},
            };
        case 2: //sector 2
            return type_region_array{
                    {xb - g[0], yc - g[1], ze - g[2], xc + g[0], ye + g[1], ze},
            };
        case 3: //sector 3
            return type_region_array{
                    {xc - g[0], yc - g[1], ze - g[2], xe + g[0], ye + g[1], ze},
            };
        case 4: //sector 4
            return type_region_array{
                    {xb - g[0], yb - g[1], zb, xc + g[0], yc + g[1], zb + g[2]},
            };
        case 5: //sector 5
            return type_region_array{
                    {xc - g[0], yb - g[1], zb, xe + g[0], yc + g[1], zb + g[2]},
            };
        case 6: //sector 6
            return type_region_array{
                    {xb - g[0], yc - g[1], zb, xc + g[0], ye + g[1], zb + g[2]},
            };
        case 7: //sector 7
            return type_region_array{
                    {xc - g[0], yc - g[1], zb, xe + g[0], ye + g[1], zb + g[2]},
            };
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::regionsDimXRecv(const _type_sector_id sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                              const _type_lattice_coord xb, const _type_lattice_coord yb,
                                              const _type_lattice_coord zb, const _type_lattice_coord xc,
                                              const _type_lattice_coord yc, const _type_lattice_coord zc,
                                              const _type_lattice_coord xe, const _type_lattice_coord ye,
                                              const _type_lattice_coord ze) {
    switch (sector_id) {
        case 0: //sector 0
            return type_region_array{
                    {xb - g[0], yb,        zb,        xb, yc + g[1], zc + g[2]},
                    {xb - g[0], yb,        ze - g[2], xb, yc + g[1], ze},
                    {xb - g[0], ye - g[1], zb,        xb, ye,        zc + g[2]},
                    {xb - g[0], ye - g[1], ze - g[2], xb, ye,        ze},
            };
        case 1: //sector 1
            return type_region_array{
                    {xe, yb,        zb,        xe + g[0], yc + g[1], zc + g[2]},
                    {xe, yb,        ze - g[2], xe + g[0], yc + g[1], ze},
                    {xe, ye - g[1], zb,        xe + g[0], ye,        zc + g[2]},
                    {xe, ye - g[1], ze - g[2], xe + g[0], ye,        ze},
            };
        case 2: //sector 2
            return type_region_array{
                    {xb - g[0], yb,        zb,        xb, yb + g[1], zc + g[2]},
                    {xb - g[0], yb,        ze - g[2], xb, yb + g[1], ze},
                    {xb - g[0], yc - g[1], zb,        xb, ye,        zc + g[2]},
                    {xb - g[0], yc - g[1], ze - g[2], xb, ye,        ze},
            };
        case 3: //sector 3
            return type_region_array{
                    {xe, yb,        zb,        xe + g[0], yb + g[1], zc + g[2]},
                    {xe, yb,        ze - g[2], xe + g[0], yb + g[1], ze},
                    {xe, yc - g[1], zb,        xe + g[0], ye,        zc + g[2]},
                    {xe, yc - g[1], ze - g[2], xe + g[0], ye,        ze},
            };
        case 4: //sector 4
            return type_region_array{
                    {xb - g[0], yb,        zb,        xb, yc + g[1], zb + g[2]},
                    {xb - g[0], yb,        zc - g[2], xb, yc + g[1], ze},
                    {xb - g[0], ye - g[1], zb,        xb, ye,        zb + g[2]},
                    {xb - g[0], ye - g[1], zc - g[2], xb, ye,        ze},
            };
        case 5: //sector 5
            return type_region_array{
                    {xe, yb,        zb,        xe + g[0], yc + g[1], zb + g[2]},
                    {xe, yb,        zc - g[2], xe + g[0], yc + g[1], ze},
                    {xe, ye - g[1], zb,        xe + g[0], ye,        zb + g[2]},
                    {xe, ye - g[1], zc - g[2], xe + g[0], ye,        ze},
            };
        case 6: //sector 6
            return type_region_array{
                    {xb - g[0], yb,        zb,        xb, yb + g[1], zb + g[2]},
                    {xb - g[0], yb,        zc - g[2], xb, yb + g[1], ze},
                    {xb - g[0], yc - g[1], zb,        xb, ye,        zb + g[2]},
                    {xb - g[0], yc - g[1], zc - g[2], xb, ye,        ze},
            };
        case 7: //sector 7
            return type_region_array{
                    {xe, yb,        zb,        xe + g[0], yb + g[1], zb + g[2]},
                    {xe, yb,        zc - g[2], xe + g[0], yb + g[1], ze},
                    {xe, yc - g[1], zb,        xe + g[0], ye,        zb + g[2]},
                    {xe, yc - g[1], zc - g[2], xe + g[0], ye,        ze},
            };
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::regionsDimYRecv(const _type_sector_id sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                              const _type_lattice_coord xb, const _type_lattice_coord yb,
                                              const _type_lattice_coord zb, const _type_lattice_coord xc,
                                              const _type_lattice_coord yc, const _type_lattice_coord zc,
                                              const _type_lattice_coord xe, const _type_lattice_coord ye,
                                              const _type_lattice_coord ze) {
    switch (sector_id) {
        case 0: //sector 0
            return type_region_array{
                    {xb - g[0], yb - g[1], zb,        xc + g[0], yb, zc + g[2]},
                    {xb - g[0], yb - g[1], ze - g[2], xc + g[0], yb, ze},
            };
        case 1: //sector 1
            return type_region_array{
                    {xc - g[0], yb - g[1], zb,        xe + g[0], yb, zc + g[2]},
                    {xc - g[0], yb - g[1], ze - g[2], xe + g[0], yb, ze},
            };
        case 2: //sector 2
            return type_region_array{
                    {xb - g[0], ye, zb,        xc + g[0], ye + g[1], zc + g[2]},
                    {xb - g[0], ye, ze - g[2], xc + g[0], ye + g[1], ze},
            };
        case 3: //sector 3
            return type_region_array{
                    {xc - g[0], ye, zb,        xe + g[0], ye + g[1], zc + g[2]},
                    {xc - g[0], ye, ze - g[2], xe + g[0], ye + g[1], ze},
            };
        case 4: //sector 4
            return type_region_array{
                    {xb - g[0], yb - g[1], zb,        xc + g[0], yb, zb + g[2]},
                    {xb - g[0], yb - g[1], zc - g[2], xc + g[0], yb, ze},
            };
        case 5: //sector 5
            return type_region_array{
                    {xc - g[0], yb - g[1], zb,        xe + g[0], yb, zb + g[2]},
                    {xc - g[0], yb - g[1], zc - g[2], xe + g[0], yb, ze},
            };
        case 6: //sector 6
            return type_region_array{
                    {xb - g[0], ye, zb,        xc + g[0], ye + g[1], zb + g[2]},
                    {xb - g[0], ye, zc - g[2], xc + g[0], ye + g[1], ze},
            };
        case 7: //sector 7
            return type_region_array{
                    {xc - g[0], ye, zb,        xe + g[0], ye + g[1], zb + g[2]},
                    {xc - g[0], ye, zc - g[2], xe + g[0], ye + g[1], ze},
            };
        default:
            assert(false);
            return type_region_array{};
    }
}

comm::type_region_array comm::regionsDimZRecv(const _type_sector_id sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                              const _type_lattice_coord xb, const _type_lattice_coord yb,
                                              const _type_lattice_coord zb, const _type_lattice_coord xc,
                                              const _type_lattice_coord yc, const _type_lattice_coord zc,
                                              const _type_lattice_coord xe, const _type_lattice_coord ye,
                                              const _type_lattice_coord ze) {
    switch (sector_id) {
        case 0: //sector 0
            return type_region_array{
                    {xb - g[0], yb - g[1], zb - g[2], xc + g[0], yc + g[1], zb},
            };
        case 1: //sector 1
            return type_region_array{
                    {xc - g[0], yb - g[1], zb - g[2], xe + g[0], yc + g[1], zb},
            };
        case 2: //sector 2
            return type_region_array{
                    {xb - g[0], yc - g[1], zb - g[2], xc + g[0], ye + g[1], zb},
            };
        case 3: //sector 3
            return type_region_array{
                    {xc - g[0], yc - g[1], zb - g[2], xe + g[0], ye + g[1], zb},
            };
        case 4: //sector 4
            return type_region_array{
                    {xb - g[0], yb - g[1], ze, xc + g[0], yc + g[1], ze + g[2]},
            };
        case 5: //sector 5
            return type_region_array{
                    {xc - g[0], yb - g[1], ze, xe + g[0], yc + g[1], ze + g[2]},
            };
        case 6: //sector 6
            return type_region_array{
                    {xb - g[0], yc - g[1], ze, xc + g[0], ye + g[1], ze + g[2]},
            };
        case 7: //sector 7
            return type_region_array{
                    {xc - g[0], yc - g[1], ze, xe + g[0], ye + g[1], ze + g[2]},
            };
        default:
            assert(false);
            return type_region_array{};
    }
}
