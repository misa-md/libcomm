//
// Created by genshen on 2019-04-19.
//

#ifndef COMM_BCC_DOMAIN_H
#define COMM_BCC_DOMAIN_H

#include "domain.h"

namespace comm {
    /**
     * additional domain data for BCC structural lattice.
     * The main different between Domain and BccDomain is that:
     * 1. add double x lattice size.
     * 2. measured length is not the same as origin value in Domain class.
     */
    class BccDomain : public Domain {
    public:
        class Builder;

    public:
        const Region<_type_lattice_coord> &dbx_local_ghost_lattice_coord_region = _dbx_local_ghost_lattice_coord_region;

        const Region<_type_lattice_coord> &dbx_local_sub_box_lattice_coord_region = _dbx_local_sub_box_lattice_coord_region;

        const Region<_type_lattice_coord> &dbx_lattice_coord_ghost_region = _dbx_lattice_coord_ghost_region;

        const Region<_type_lattice_coord> &dbx_lattice_coord_sub_box_region = _dbx_lattice_coord_sub_box_region;

        const _type_lattice_size (&dbx_lattice_size_sub_box)[DIMENSION_SIZE] = _dbx_lattice_sub_box_size;

        const _type_lattice_size (&dbx_lattice_size_ghost_extended)[DIMENSION_SIZE] = _dbx_lattice_size_ghost_extended;

        const _type_lattice_size (&dbx_lattice_size_ghost)[DIMENSION_SIZE] = _dbx_lattice_size_ghost;

        /**
         * construct bcc domain.
         * \param _phase_space the space size in each dimension.
         * \param _lattice_const the measured length of between two lattice(grid)
         * \param _cutoff_radius_factor cutoff for building ghost area:
         *  _cutoff_radius_factor = measured ghost length /_lattice_const
         */
        BccDomain(const std::array<uint64_t, DIMENSION_SIZE> _phase_space,
                  const double _lattice_const, const double _cutoff_radius_factor);

        /**
         * create bcc domain from existed normal domain.
         * \param domain ref of normal domain.
         */
        explicit BccDomain(const Domain &domain);

        /**
         * renew the value in class BCCDomain.
         */
        void rescale(const comm::Domain &domain);

    private:
        // doubled x size for BCC lattice.
        _type_lattice_size _dbx_lattice_sub_box_size[DIMENSION_SIZE];
        _type_lattice_size _dbx_lattice_size_ghost_extended[DIMENSION_SIZE];
        _type_lattice_size _dbx_lattice_size_ghost[DIMENSION_SIZE];
        Region<_type_lattice_coord> _dbx_lattice_coord_sub_box_region;
        Region<_type_lattice_coord> _dbx_lattice_coord_ghost_region;

        Region<_type_lattice_coord> _dbx_local_sub_box_lattice_coord_region;
        Region<_type_lattice_coord> _dbx_local_ghost_lattice_coord_region;
    };

    /**
     * builder for building bcc domain
     */
    class BccDomain::Builder : public comm::Builder<BccDomain::Builder, BccDomain> {
    public:
        BccDomain *build() override;

        BccDomain *localBuild(const int _grid_size[DIMENSION_SIZE],
                              const int _grid_coord[DIMENSION_SIZE]) override;

    private:
        typedef comm::Builder<BccDomain::Builder, BccDomain> parent_domain;
    };

}

#endif //COMM_BCC_DOMAIN_H
