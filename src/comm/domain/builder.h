//
// Created by genshen on 2019-04-19.
//

#ifndef COMM_BUILDER_H
#define COMM_BUILDER_H

#include <cstdint>
#include <array>
#include <cmath>
#include <mpi.h>
//#include <utils/mpi_utils.h>

#include "comm/types_define.h"

namespace comm {
  /**
   * \brief builder for building general domain.
   * we use [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) to implement this feature.
   * \tparam B Build type, derived class of com::Builder
   * \tparam D Domain type, comm::Domain or derived class of comm::Domain
   */
  template <typename B, typename D> class Builder {
  public:
    /**
     * set mpi rank and communications
     * @param mpi_process current MPI rank id, the ranks in current communicator and communicator
     * @param comm the new communicator after decomposition.
     * @return
     */
    B &setComm(mpi_process mpi_process, MPI_Comm *comm);

    B &setPhaseSpace(const int64_t phaseSpace[DIMENSION_SIZE]);

    B &setLatticeConst(const double latticeConst);

    /**
     * set ghost size
     * \param ghost_size the size of ghost size, unit: lattice, default cut_lattice
     * (should set cut_lattice before setting ghost size).
     * \return reference of Builder.
     */
    B &setGhostSize(const unsigned int ghost_size);

    B &setCutoffRadius(const double cutoff_radius_factor);

    B &setMPIMap3dSubDim(const int mpi_map_3d_sub_dim[DIMENSION_SIZE]);

    /**
     * remember to delete it when it is used
     * @return pointer to domain.
     */
    virtual D *build() = 0;

    /** local builder can build domain for one processor, without connection to other processors.
     * \note rank_id_neighbours is not set in local build.
     * @param _grid_size user defined grid size of decomposition.
     * @param _grid_coord the grid coordinate of current local domain/sub-box.
     * @return pointer to the Domain of sub-box.
     */
    virtual D *localBuild(const int _grid_size[DIMENSION_SIZE], const int _grid_coord[DIMENSION_SIZE]) = 0;

  protected:
    mpi_process _mpi_pro;
    MPI_Comm *_p_comm;
    double _cutoff_radius_factor;
    double _lattice_const;
    int _ghost_size;
    std::array<uint64_t, DIMENSION_SIZE> _phase_space;

    /**
     * In this method, each processor will be bound to a cartesian coordinate.
     *
     * It first divide the simulation box into N pieces(sub-box) (N is the count of all processors).
     * And each processor will be bound to a sub-box, and tagged with a cartesian coordinate(x,y,z).
     */
    virtual void decomposition(D &domain);
    void decomposition_imp_v1(D &domain);
    void decomposition_imp_v2(D &domain);

    /**
     * set length of global simulation box.
     * and set upper and lower bound of global simulation box.
     * @param domain reference to domain
     */
    virtual void createGlobalDomain(D &domain);

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
    virtual void buildLatticeDomain(D &domain);

    /**
     * set lattice coordinate boundary of current sub-box in local coordinate system(LCY).
     */
    virtual void buildMeasuredDomain(D &domain); // todo test.

    /*
     * the 3D mapping of MPI process in node, default: 1*1*1
     */
    int _mpi_map_3d_sub_dim[DIMENSION_SIZE] = {1, 1, 1};
  };
} // namespace comm

#include "builder.inl"

#endif // COMM_BUILDER_H
