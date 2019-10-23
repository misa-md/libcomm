//
// Created by genshen on 2019-04-15.
//

#ifndef COMM_COMM_FORWARDING_REGION_H
#define COMM_COMM_FORWARDING_REGION_H

#include <cassert>
#include "comm/domain/region.hpp"
#include "comm/domain/bcc_domain.h"

/**
 * this header file describes the communication region
 * when performing communication forwarding in x,y,z dimension.
 *
 */
namespace comm {
    /**
     * This function returns communication region when performing communication forwarding.
     * The region is the area belongs to current process,
     * but has contribution to its corresponding neighbour processes.
     *
     * \param p_domain pointer to the domain.
     * \param dimension dimension for communication, 0 for x dimension, 1 for y dimension, 2 for z dimension
     * \param direction direction for communication, values: DIR_LOWER or DIR_HIGHER.
     * \return region for communication forwarding, unit: lattice size.
     * \note: the region in x dimension in return value is double due to BCC lattice structure.
     */
    Region <_type_lattice_size> fwCommLocalRegion(const BccDomain *p_domain, const unsigned int dimension,
                                                  const unsigned int direction);

    /**
     * This function returns communication region when performing communication forwarding.
     * The unit is measure length, not lattice size as above function.
     * \deprecated do not use measured length to pack atoms for communication
     * \param p_domain pointer to the domain.
     * \param dimension dimension for communication, 0 for x dimension, 1 for y dimension, 2 for z dimension
     * \param direction direction for communication, values: DIR_LOWER or DIR_HIGHER.
     * \return region for communication forwarding, unit: measured length, which is double.
     */
    Region<double> fwCommLocalMeaRegion(const Domain *p_domain, const unsigned int dimension,
                                        unsigned const int direction);

}

#endif //COMM_COMM_FORWARDING_REGION_H
