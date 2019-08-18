//
// Created by genshen on 2019-08-18.
//

#ifndef COMM_SECTOR_FORWARDING_REGION_H
#define COMM_SECTOR_FORWARDING_REGION_H


#include <types_define.h>
#include <domain/colored_domain.h>

namespace comm {

    typedef std::vector<comm::Region<comm::_type_lattice_size>> type_region_array;

    /**
     * return table of sending regions of each sectors in each dimension.
     * \param p_domain pointer of domain.
     * \param sector_id sector id from 0 to 7.
     * \param dimension dimension (can be DIM_X, DIM_Y, DIM_Z).
     * \return regions to be sent.
     */
    type_region_array fwCommSectorLocalRegion(
            const ColoredDomain *p_domain, const unsigned int sector_id,
            const unsigned int dimension);
}

#endif //COMM_SECTOR_FORWARDING_REGION_H
