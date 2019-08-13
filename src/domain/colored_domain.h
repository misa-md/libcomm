//
// Created by genshen on 2019-07-21.
//

#ifndef COMM_COLORED_DOMAIN_H
#define COMM_COLORED_DOMAIN_H

#include "domain.h"

namespace comm {
    const unsigned int X_LOW = 0;
    const unsigned int X_HIGH = 1;
    const unsigned int Y_LOW = 0;
    const unsigned int Y_HIGH = 2;
    const unsigned int Z_LOW = 0;
    const unsigned int Z_HIGH = 4;

    /**
     * Colored domain is usually used in parallel KMC sub-lattice or synchronous parallel KMC algorithm.
     */
    class ColoredDomain : public Domain {
    public:
        class Builder;

    public:
        /**
         * lattice size of each sector in this sub box.
         * each region is indexed by high/low flag, for example, binary (X_HIGH | Y_LOW | Z_HIGH) represents the region (1,0,1)
         */
        std::array<std::array<_type_lattice_size, 3>, 8> sector_lattice_size;

        /**
         * Regions of 8 sectors in sub box.
         *
         * @note the start index(zero index) is started from lower ghost boundary of local sub box.
         * index denotation of local sub-box as well as its ghost, this is in local box Coordinate System(not global box).
         *  | 0:ghost_lower                 | sub_box_lower                   | sub_box_upper               | ghost_upper
         *  |-------------------------------|---------...---------------------|-----------------------------|
         *  for example, a 8*8*8 sub-box with ghost size 3, a sector whose coordinate is (0,1,0) will
         *  have region x from 3 to 7 (3+8/2=7) and y from 7 to 11 (7+4 =11), and z from 3 to 7.
         */
        std::array<Region<_type_lattice_size>, 8> local_sector_region;

        /**
         * Regions of 8 sectors plus its ghost area in sub box.
         * For example, a 8*8*8 sub-box with ghost size 3, a sector whose coordinate is (0,1,0) will
         * have region x from 0 to 10 (3+8/2+3=10) and y from 4 (3+8/2-3=4) to 14 (3+8+3=14), and z from 0 to 10.
         */
        std::array<Region<_type_lattice_size>, 8> local_sector_ghost_region;

        /**
         * construct colored domain.
         * \param _phase_space the space size in each dimension.
         * \param _lattice_const the measured length of between two lattice(grid)
         * \param _cutoff_radius_factor cutoff for building ghost area:
         *  _cutoff_radius_factor = measured ghost length /_lattice_const
         */
        ColoredDomain(const std::array<uint64_t, DIMENSION_SIZE> _phase_space,
                      const double _lattice_const, const double _cutoff_radius_factor);

        /**
        * create colored domain from existed parent domain.
        * \param domain reference of parent domain.
        */
        explicit ColoredDomain(const Domain &domain);

        /**
         * set lattice size and region for each sector in sub-box.
         * \param domain reference of domain.
         */
        void splitSector(const comm::Domain &domain);
    };

    /**
     * builder for building colored domain
     */
    class ColoredDomain::Builder : public comm::Builder<ColoredDomain::Builder, ColoredDomain> {
    public:
        /**
         * builder domain using MPI_Cart_create.
         * \return
         */
        ColoredDomain *build() override;

        /**
         * build a local domain for a process.
         * \param _grid_size grid size in 3d
         * \param _grid_coord the grid coord of current process.
         * \return
         */
        ColoredDomain *localBuild(const int _grid_size[DIMENSION_SIZE],
                                  const int _grid_coord[DIMENSION_SIZE]) override;

    private:
        typedef comm::Builder<ColoredDomain::Builder, ColoredDomain> parent_domain;
    };

}

#endif //COMM_COLORED_DOMAIN_H
