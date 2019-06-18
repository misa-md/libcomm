//
// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.
//

#ifndef COMM_DOMAIN_H
#define COMM_DOMAIN_H

#include <mpi.h>
#include <vector>
#include <array>

#include "types_define.h"
#include "region.hpp"
#include "builder.h"

/**
 * If N can be decomposed as N = N_x * N_y * N_z, where N, N_x, N_y, N_z are all integer bigger than or equal to 1,
 * then the whole simulation box will be divided into N sub-box with N_z levels in z axis,
 * and in each level, it has  N_x * N_y sub-boxes.
 * Then, we can bind each processor to a cartesian coordinate (x,y,z) due to the boxes partition,
 * where 0 <= x < N_x, 0 <= z < N_z, 0 <= z < N_z.
 * Last, based on the cartesian coordinate (x,y,z),
 * each processor can get the cartesian coordinate of its contiguous sub-boxes.
 *
 *
 * variable naming:
 * 1.variables for real length(measured length) and real boundary(measured boundary) of sub-box or global box
 * have a prefix "meas" or "_meas", with double type.
 *
 * 2.variables for cartesian coordinate of box decomposition have a prefix "grid_coord" or "_grid_coord", with int type;
 * and the grid count of box decomposition at each dimension have a prefix "grid_size" or "_grid_size", with int type;
 *
 * 3.variables for lattice coordinate have prefix "lattice_coord" or "_lattice_coord", with int type.
 *
 * 4.variables for lattice count in box or in global box have a prefix "lattice_size" or "_lattice_size", with int type.
 */

namespace comm {
    class Domain {
    public:
        class Builder;

        template<typename, typename>
        friend
        class comm::Builder;

    public:
        const double lattice_const;
        const double cutoff_radius_factor;
        // cut off lattice size.
        const _type_lattice_size cut_lattice;
        const std::array<u_int64_t, DIMENSION_SIZE> phase_space;
        /**
         * global measured length of the simulation box at each dimension.
         */
        const double (&meas_global_length)[DIMENSION_SIZE] = _meas_global_length;

        /**
         * the count of processors at each dimension.
         * Or we can say the decomposed grid size at each dimension.
         */
        const int (&grid_size)[DIMENSION_SIZE] = _grid_size;

        /**
         * the measured coordinate of lower and upper boundary of global simulation box.
         */
        const Region<double> &meas_global_box_coord_region = _meas_global_box_coord_region;

        /** local information for current simulation sub-box. **/
        /**
         * The cartesian coordinate of the sub-box bound to this processor after running grid decomposition.
         */
        const int (&grid_coord_sub_box)[DIMENSION_SIZE] = _grid_coord_sub_box;

        /**
         * the rank ids of contiguous processors in space.
         */
        const _MPI_Rank (&rank_id_neighbours)[DIMENSION_SIZE][2] = _rank_id_neighbours;

        /** boundary of local sub-box  **/
        /**
         * the measured lower and upper boundary of current sub-box at each dimension in global box Coordinate System.
         */
        const Region<double> &meas_sub_box_region = _meas_sub_box_region;

        /**boundary of ghost of local sub-box**/
        /**
         * measured ghost length at each dimension, which equals to the cutoff radius.
         */
        const double (&meas_ghost_length)[DIMENSION_SIZE] = _meas_ghost_length;
        /**
         * the measured ghost lower and upper bound of current sub-box.
         */
        const Region<double> &meas_ghost_region = _meas_ghost_region;

        /*lattice count in local sub-box*/
        /**
         * lattice count in local sub-box area at each dimension (upper boundary - lower boundary).
         */
        const _type_lattice_size (&lattice_size_sub_box)[DIMENSION_SIZE] = _lattice_sub_box_size;

        /**
         * lattice count in ghost area plus sub-box area at each dimension (upper boundary - lower boundary).
         * which  @var _lattice_size_ghost_extended[d] = @var _lattice_size_sub_box[d] + 2 * @var_ lattice_size_ghost[d];
         * and also equals to _lattice_coord_ghost_region.high[d] - _lattice_coord_ghost_region.low[d]
         */
        const _type_lattice_size (&lattice_size_ghost_extended)[DIMENSION_SIZE] = _lattice_size_ghost_extended;

        /**
         * purge ghost size, just lattice count in ghost area.
         */
        const _type_lattice_size (&lattice_size_ghost)[DIMENSION_SIZE] = _lattice_size_ghost;

        /*lattice boundary of local sub-box and ghost, but, the Coordinate System is still the global box.*/
        /**
         * lower and upper boundary(not included) of lattice coordinate of current local sub-box at each dimension
         * in global coordinate system(GCY).
         *
         */
        const Region<_type_lattice_coord> &lattice_coord_sub_box_region = _lattice_coord_sub_box_region;

        /**
         * lower and upper boundary(not included) of lattice coordinate in ghost area of current sub-box area
         * at each dimension in global coordinate system(GCY)
         */
        const Region<_type_lattice_coord> &lattice_coord_ghost_region = _lattice_coord_ghost_region;

        /*
         * lattice boundary of local sub-box and ghost, this is in local box Coordinate System(not global box.).
         * For convenience usage for atom index in local sub-box.
         *  | 0:ghost_lower                 | sub_box_lower                   | sub_box_upper               | ghost_upper
         *  |-------------------------------|---------...---------------------|-----------------------------|
         */
        /**
         * lower and upper boundary(not included) of lattice coordinate of local sub-box
         * at each dimension in local coordinate system(LCY).
         */
        const Region<_type_lattice_coord> &local_sub_box_lattice_coord_region = _local_sub_box_lattice_coord_region;

        // lower and upper boundary(not included) of lattice coordinate in ghost area of current sub-box area
        // at each dimension in local coordinate system(LCY).
        const Region<_type_lattice_coord> &local_ghost_lattice_coord_region = _local_ghost_lattice_coord_region;

    protected:
        Domain(const std::array<u_int64_t, DIMENSION_SIZE> _phase_space,
               const double _lattice_const, const double _cutoff_radius_factor);

        /** the private variables are referenced in preview public filed.*/
        double _meas_global_length[DIMENSION_SIZE];
        int _grid_size[DIMENSION_SIZE] = {0};

        Region<double> _meas_global_box_coord_region;

        int _grid_coord_sub_box[DIMENSION_SIZE];
        _MPI_Rank _rank_id_neighbours[DIMENSION_SIZE][2];

        /**
         * the measured lower bound of current sub-box at a dimension
         */
        Region<double> _meas_sub_box_region;
        double _meas_ghost_length[DIMENSION_SIZE];  // measured ghost length, which equals to the cutoff radius.
        Region<double> _meas_ghost_region; // measured sub box region plus measured ghost length.

        _type_lattice_size _lattice_sub_box_size[DIMENSION_SIZE];
        _type_lattice_size _lattice_size_ghost_extended[DIMENSION_SIZE];
        _type_lattice_size _lattice_size_ghost[DIMENSION_SIZE];
        Region<_type_lattice_coord> _lattice_coord_sub_box_region;
        Region<_type_lattice_coord> _lattice_coord_ghost_region;
        Region<_type_lattice_coord> _local_sub_box_lattice_coord_region;
        Region<_type_lattice_coord> _local_ghost_lattice_coord_region;
    };

    class Domain::Builder : public comm::Builder<Domain::Builder, Domain> {
    public:
        Domain *build() override;

        Domain *localBuild(const int _grid_size[DIMENSION_SIZE], const int _grid_coord[DIMENSION_SIZE]) override;
    };
}

#endif // COMM_DOMAIN_H
