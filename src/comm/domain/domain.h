//
// Created by baihe back to 2016-12-22.
// refactored by genshen on 2018-12-31.
//

#ifndef COMM_DOMAIN_H
#define COMM_DOMAIN_H

#include <array>
#include <vector>

#include "builder.h"
#include "cartesian_context.h"
#include "comm/types_define.h"
#include "region.hpp"

/**
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
  class MeasuredDomain {
  public:
    /**
     * global measured length of the simulation box at each dimension.
     */
    const double (&meas_global_length)[DIMENSION_SIZE] = _meas_global_length;

    /**
     * the measured coordinate of lower and upper boundary of global simulation box.
     */
    const Region<double> &meas_global_region = _meas_global_region;

    /** boundary of local sub-box  **/
    /**
     * the global measured lower and upper boundary of current sub-box at each dimension in global box Coordinate
     * System.
     */
    const Region<double> &meas_sub_box_region = _meas_sub_box_region;

    /**boundary of ghost of local sub-box**/
    /**
     * measured ghost length at each dimension, which equals to the cutoff radius.
     */
    const double (&meas_ghost_length)[DIMENSION_SIZE] = _meas_ghost_length;
    /**
     * the global measured ghost lower and upper bound of current sub-box.
     */
    const Region<double> &meas_ghost_ext_region = _meas_ghost_ext_region;

    inline double volume() const { return _meas_global_length[0] * _meas_global_length[1] * _meas_global_length[2]; }
    inline double Lx() const { return _meas_global_length[0]; }
    inline double Ly() const { return _meas_global_length[1]; }
    inline double Lz() const { return _meas_global_length[2]; }

  protected:
    /** the private variables are referenced in preview public filed.*/
    double _meas_global_length[DIMENSION_SIZE] = {0.0, 0.0, 0.0};
    Region<double> _meas_global_region; // the start and stop position of the global simulation box.

    /**
     * the measured lower bound of current sub-box at a dimension
     */
    Region<double> _meas_sub_box_region;
    /**
     * measured ghost length, which equals to the cutoff radius.
     */
    double _meas_ghost_length[DIMENSION_SIZE] = {0.0, 0.0, 0.0};
    Region<double> _meas_ghost_ext_region; // measured sub box region plus measured ghost length.

    void rescale_measured(const std::array<double, DIMENSION_SIZE>);
  };

  class LatticeDomain : public MeasuredDomain, public CartesianContext {
  public:
    template <typename, typename> friend class comm::Builder;

    class Builder;

  public:
    std::array<double, DIMENSION_SIZE> lattice_const; // the lattice constant

    const double cutoff_radius;
    /**
     * cut off lattice size.
     */
    const _type_lattice_size cut_lattice;
    const std::array<uint64_t, DIMENSION_SIZE> phase_space;

    inline double lattice_const_1d() const { return lattice_const[0]; }

    /**
     * @deprecated
     * _cutoff_radius_factor = measured ghost length /_lattice_const
     */
    inline double cutoff_radius_factor() const { return cutoff_radius / lattice_const_1d(); }

    /*lattice count in local sub-box*/
    /**
     * lattice count in local sub-box area at each dimension (upper boundary - lower boundary).
     */
    const _type_lattice_size (&sub_box_lattice_size)[DIMENSION_SIZE] = _sub_box_lattice_size;

    /**
     * lattice count in ghost area plus sub-box area at each dimension (upper boundary - lower boundary).
     * which  @var _lattice_size_ghost_extended[d] = @var _lattice_size_sub_box[d] + 2 * @var_ lattice_size_ghost[d];
     * and also equals to _lattice_coord_ghost_region.high[d] - _lattice_coord_ghost_region.low[d]
     */
    const _type_lattice_size (&ghost_extended_lattice_size)[DIMENSION_SIZE] = _ghost_extended_lattice_size;

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
    const Region<_type_lattice_coord> &sub_box_lattice_region = _sub_box_lattice_region;

    /**
     * lower and upper boundary(not included) of lattice coordinate with ghost area of current sub-box area
     * at each dimension in global coordinate system(GCY)
     */
    const Region<_type_lattice_coord> &ghost_ext_lattice_region = _ghost_ext_lattice_region;

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
    const Region<_type_lattice_coord> &local_sub_box_lattice_region = _local_sub_box_lattice_region;

    /**
     * lower and upper boundary(not included) of lattice coordinate in ghost area of current sub-box area
     * at each dimension in local coordinate system(LCY).
     */
    const Region<_type_lattice_coord> &local_ghost_ext_lattice_region = _local_ghost_ext_lattice_region;

    /**
     * rescale the domain with given scale factor.
     * make the coordinates and lengths multiplied by the scale factor,
     * which is useful for NPT ensemble in the MD simulation.
     * @param scale_factor the scale factor to rescale the domain.
     * @return the rescaled new domain.
     * @note rescale only rescale the measured length and regions.
     * the lattice size/region and ghost size/region (including the measured ghost size/region) will remain unchanged.
     */
    void rescale(const double scale_factor);

    void rescale(const std::array<double, DIMENSION_SIZE> scale_factors);

  protected:
    LatticeDomain(const std::array<uint64_t, DIMENSION_SIZE> _phase_space,
                  const std::array<double, DIMENSION_SIZE> _lattice_const, const double _cutoff_radius_factor);

    /** the private variables are referenced in preview public filed.*/
    _type_lattice_size _sub_box_lattice_size[DIMENSION_SIZE];
    _type_lattice_size _ghost_extended_lattice_size[DIMENSION_SIZE];
    _type_lattice_size _lattice_size_ghost[DIMENSION_SIZE];
    Region<_type_lattice_coord> _sub_box_lattice_region;
    Region<_type_lattice_coord> _ghost_ext_lattice_region;
    Region<_type_lattice_coord> _local_sub_box_lattice_region;
    Region<_type_lattice_coord> _local_ghost_ext_lattice_region;
  };

  using Domain = LatticeDomain;

  class Domain::Builder : public comm::Builder<Domain::Builder, Domain> {
  public:
    Domain *build() override;

    Domain *localBuild(const int _grid_size[DIMENSION_SIZE], const int _grid_coord[DIMENSION_SIZE]) override;
  };
} // namespace comm

#endif // COMM_DOMAIN_H
