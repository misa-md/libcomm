//
// Created by genshen on 2019-08-18.
//

#ifndef COMM_SECTOR_FORWARDING_REGION_H
#define COMM_SECTOR_FORWARDING_REGION_H


#include <types_define.h>
#include <domain/colored_domain.h>

namespace comm {

    typedef std::vector<comm::Region<comm::_type_lattice_coord>> type_region_array;

    /**
     * return table of sending regions of each sectors in each dimension.
     * \param sector_id sector id from 0 to 7.
     * \param dim dimension (can be DIM_X, DIM_Y, DIM_Z).
     * \param ghost_size ghost_size of each dimension.
     * \param split_coord The split coordinate to split 8 sectors.
     * \param local_box_region local region of current sub box.
     * \return regions to be sent.
     */
    type_region_array fwCommSectorSendRegion(const unsigned int sector_id, const unsigned int dim,
                                             const _type_lattice_size ghost_size[DIMENSION_SIZE],
                                             const _type_lattice_coord split_coord[DIMENSION_SIZE],
                                             const Region<comm::_type_lattice_coord> local_box_region);

    /**
     * similar as above, but it returns receive regions for specified sector at a dimension.
     * \return receive regions for specified sector at a dimension.
     */
    type_region_array fwCommSectorRecvRegion(const unsigned int sector_id, const unsigned int dim,
                                             const _type_lattice_size ghost_size[DIMENSION_SIZE],
                                             const _type_lattice_coord split_coord[DIMENSION_SIZE],
                                             const Region<comm::_type_lattice_coord> local_box_region);

    // private methods
    type_region_array regionsDimXSend(const unsigned int sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                      const _type_lattice_coord xb, const _type_lattice_coord yb,
                                      const _type_lattice_coord zb, const _type_lattice_coord xc,
                                      const _type_lattice_coord yc, const _type_lattice_coord zc,
                                      const _type_lattice_coord xe, const _type_lattice_coord ye,
                                      const _type_lattice_coord ze);

    type_region_array regionsDimYSend(const unsigned int sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                      const _type_lattice_coord xb, const _type_lattice_coord yb,
                                      const _type_lattice_coord zb, const _type_lattice_coord xc,
                                      const _type_lattice_coord yc, const _type_lattice_coord zc,
                                      const _type_lattice_coord xe, const _type_lattice_coord ye,
                                      const _type_lattice_coord ze);

    type_region_array regionsDimZSend(const unsigned int sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                      const _type_lattice_coord xb, const _type_lattice_coord yb,
                                      const _type_lattice_coord zb, const _type_lattice_coord xc,
                                      const _type_lattice_coord yc, const _type_lattice_coord zc,
                                      const _type_lattice_coord xe, const _type_lattice_coord ye,
                                      const _type_lattice_coord ze);

    type_region_array regionsDimXRecv(const unsigned int sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                      const _type_lattice_coord xb, const _type_lattice_coord yb,
                                      const _type_lattice_coord zb, const _type_lattice_coord xc,
                                      const _type_lattice_coord yc, const _type_lattice_coord zc,
                                      const _type_lattice_coord xe, const _type_lattice_coord ye,
                                      const _type_lattice_coord ze);

    type_region_array regionsDimYRecv(const unsigned int sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                      const _type_lattice_coord xb, const _type_lattice_coord yb,
                                      const _type_lattice_coord zb, const _type_lattice_coord xc,
                                      const _type_lattice_coord yc, const _type_lattice_coord zc,
                                      const _type_lattice_coord xe, const _type_lattice_coord ye,
                                      const _type_lattice_coord ze);

    type_region_array regionsDimZRecv(const unsigned int sector_id, const _type_lattice_size g[DIMENSION_SIZE],
                                      const _type_lattice_coord xb, const _type_lattice_coord yb,
                                      const _type_lattice_coord zb, const _type_lattice_coord xc,
                                      const _type_lattice_coord yc, const _type_lattice_coord zc,
                                      const _type_lattice_coord xe, const _type_lattice_coord ye,
                                      const _type_lattice_coord ze);

}

#endif //COMM_SECTOR_FORWARDING_REGION_H
