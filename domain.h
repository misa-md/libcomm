//
// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.
//

#ifndef CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
#define CRYSTAL_MD_DOMAIN_DECOMPOSITION_H

#include <mpi.h>
#include <vector>
#include <array>
#include <utils/mpi_utils.h>

#include "types/pre_define.h"
#include "region.hpp"

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

#define LOWER  0
#define HIGHER 1


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

class Domain {
public:
    const double lattice_const;
    const double cutoff_radius_factor;
    // cut off lattice size.
    const _type_lattice_size cut_lattice;
    const std::array<u_int64_t, DIMENSION> phase_space;
    /**
     * global measured length of the simulation box at each dimension.
     */
    const double (&meas_global_length)[DIMENSION];

    /**
     * the count of processors at each dimension.
     * Or we can say the decomposed grid size at each dimension.
     */
    const int (&grid_size)[DIMENSION] = {0};

    /**
     * the measured coordinate of lower and upper boundary of global simulation box.
     */
    const Region<double> &meas_global_box_coord_region;

    /** local information for current simulation sub-box. **/
    /**
     * The cartesian coordinate of the sub-box bound to this processor after running grid decomposition.
     */
    const int (&grid_coord_sub_box)[DIMENSION];

    /**
     * the rank ids of contiguous processors in space.
     */
    const int (&rank_id_neighbours)[DIMENSION][2];

    /** boundary of local sub-box  **/
    /**
     * the measured lower and upper boundary of current sub-box at each dimension in global box Coordinate System.
     */
    const Region<double> &meas_sub_box_region;

    /**boundary of ghost of local sub-box**/
    /**
     * measured ghost length at each dimension, which equals to the cutoff radius.
     */
    const double (&meas_ghost_length)[DIMENSION];
    /**
     * the measured ghost lower and upper bound of current sub-box.
     */
    const Region<double> &meas_ghost_region;

    /*lattice count in local sub-box*/
    /**
     * lattice count in local sub-box area at each dimension (upper boundary - lower boundary).
     */
    const _type_lattice_size (&lattice_size_sub_box)[DIMENSION];
    const _type_lattice_size (&dbx_lattice_size_sub_box)[DIMENSION];
    /**
     * lattice count in ghost area plus sub-box area at each dimension (upper boundary - lower boundary).
     * which  @var _lattice_size_ghost_extended[d] = @var _lattice_size_sub_box[d] + 2 * @var_ lattice_size_ghost[d];
     * and also equals to _lattice_coord_ghost_region.high[d] - _lattice_coord_ghost_region.low[d]
     */
    const _type_lattice_size (&lattice_size_ghost_extended)[DIMENSION];
    const _type_lattice_size (&dbx_lattice_size_ghost_extended)[DIMENSION];
    /**
     * purge ghost size, just lattice count in ghost area.
     */
    const _type_lattice_size (&lattice_size_ghost)[DIMENSION];
    const _type_lattice_size (&dbx_lattice_size_ghost)[DIMENSION];

    /*lattice boundary of local sub-box and ghost, but, the Coordinate System is still the global box.*/
    /**
     * lower and upper boundary(not included) of lattice coordinate of current local sub-box at each dimension
     * in global coordinate system(GCY).
     *
     */
    const Region<_type_lattice_coord> &lattice_coord_sub_box_region;
    const Region<_type_lattice_coord> &dbx_lattice_coord_sub_box_region;

    /**
     * lower and upper boundary(not included) of lattice coordinate in ghost area of current sub-box area
     * at each dimension in global coordinate system(GCY)
     */
    const Region<_type_lattice_coord> &lattice_coord_ghost_region;
    const Region<_type_lattice_coord> &dbx_lattice_coord_ghost_region;

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
    const Region<_type_lattice_coord> &local_sub_box_lattice_coord_region;
    const Region<_type_lattice_coord> &dbx_local_sub_box_lattice_coord_region;

    // lower and upper boundary(not included) of lattice coordinate in ghost area of current sub-box area
    // at each dimension in local coordinate system(LCY).
    const Region<_type_lattice_coord> &local_ghost_lattice_coord_region;
    const Region<_type_lattice_coord> &dbx_local_ghost_lattice_coord_region;

private:

    /** the private variables are referenced in preview public filed.*/
    double _meas_global_length[DIMENSION];
    int _grid_size[DIMENSION] = {0};

    Region<double> _meas_global_box_coord_region;

    int _grid_coord_sub_box[DIMENSION];
    int _rank_id_neighbours[DIMENSION][2];

    /**
     * the measured lower bound of current sub-box at a dimension
     */
    Region<double> _meas_sub_box_region;
    double _meas_ghost_length[DIMENSION];  // measured ghost length, which equals to the cutoff radius.
    Region<double> _meas_ghost_region;

    _type_lattice_size _lattice_sub_box_size[DIMENSION];
    _type_lattice_size _lattice_size_ghost_extended[DIMENSION];
    _type_lattice_size _lattice_size_ghost[DIMENSION];
    Region<_type_lattice_coord> _lattice_coord_sub_box_region;
    Region<_type_lattice_coord> _lattice_coord_ghost_region;
    Region<_type_lattice_coord> _local_sub_box_lattice_coord_region;
    Region<_type_lattice_coord> _local_ghost_lattice_coord_region;

    // double x
    _type_lattice_size _dbx_lattice_sub_box_size[DIMENSION];
    _type_lattice_size _dbx_lattice_size_ghost_extended[DIMENSION];
    _type_lattice_size _dbx_lattice_size_ghost[DIMENSION];
    Region<_type_lattice_coord> _dbx_lattice_coord_sub_box_region;
    Region<_type_lattice_coord> _dbx_lattice_coord_ghost_region;

    Region<_type_lattice_coord> _dbx_local_sub_box_lattice_coord_region;
    Region<_type_lattice_coord> _dbx_local_ghost_lattice_coord_region;

    Domain(const std::array<u_int64_t, DIMENSION> _phase_space,
           const double _lattice_const, const double _cutoff_radius_factor);

public:
    class Builder {
    public:
        /**
         * set mpi rank and communications
         * @param mpi_process current MPI rank id, the ranks in current communicator and communicator
         * @param comm the new communicator after decomposition.
         * @return
         */
        Builder &setComm(kiwi::mpi_process mpi_process, MPI_Comm *comm);

        Builder &setPhaseSpace(const int64_t phaseSpace[DIMENSION]);

        Builder &setLatticeConst(const double latticeConst);

        Builder &setCutoffRadius(const double cutoffRadius);

        /**
         * remember to delete it when it is used
         * @return pointer to domain.
         */
        Domain *build();

    private:
        kiwi::mpi_process _mpi_pro;
        MPI_Comm *_p_comm;
        double _cutoff_radius_factor;
        double _lattice_const;
        std::array<u_int64_t, DIMENSION> _phase_space;

        /**
         * In this method, each processor will be bound to a cartesian coordinate.
         *
         * It first divide the simulation box into N pieces(sub-box) (N is the count of all processors).
         * And each processor will be bound to a sub-box, and tagged with a cartesian coordinate(x,y,z).
         */
        void decomposition(Domain &domain);

        /**
         * set length of global simulation box.
         * and set upper and lower bound of global simulation box.
         * @param domain reference to domain
         */
        void createGlobalDomain(Domain &domain);

        /**
         * set boundary for current sub-box.
         *
         * Note: why some variable in x dimension is multiplied by 2:
         * the basic unit in x dimension is half the lattice, not a lattice like in y,z dimension.
         * in x dimension:
         *    | :  |  :  |  :  |  :  |  :  |  ...
         * x: 0 1  2  3  4  5  6  7  8  9  10 ...
         * but, in y or z dimension:
         *    | :  |  :  |  :  |  :  |  :  | ...
         * y: 0    1     2     3     4     5  ...
         * in above figure, |  :  | represents a lattice length.
         *
         */
        void buildLatticeDomain(Domain &domain);

        /**
         * set lattice coordinate boundary of current sub-box in local coordinate system(LCY).
         */
        void buildMeasuredDomain(Domain &domain); // todo test.
    };
};

#endif // CRYSTAL_MD_DOMAIN_DECOMPOSITION_H
