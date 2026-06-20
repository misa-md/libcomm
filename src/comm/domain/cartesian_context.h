//
// Created by genshen on 2026/6/20.
//

#ifndef COMM_CARTESIAN_CONTEXT_H
#define COMM_CARTESIAN_CONTEXT_H

#include <mpi.h>

#include "comm/types_define.h"

namespace comm {

  /**
   * store MPI rank info here:
   *
   * If N can be decomposed as N = N_x * N_y * N_z, where N, N_x, N_y, N_z are all integer bigger than or equal to 1,
   * then the whole simulation box will be divided into N sub-box with N_z levels in z axis,
   * and in each level, it has  N_x * N_y sub-boxes.
   * Then, we can bind each processor to a cartesian coordinate (x,y,z) due to the boxes partition,
   * where 0 <= x < N_x, 0 <= z < N_z, 0 <= z < N_z.
   * Last, based on the cartesian coordinate (x,y,z),
   * each processor can get the cartesian coordinate of its contiguous sub-boxes.
   *
   */
  class CartesianContext {
  public:
    /**
     * the count of processors at each dimension.
     * Or we can say the decomposed grid size at each dimension.
     */
    const int (&grid_size)[DIMENSION_SIZE] = _grid_size;

    /** local information for current simulation sub-box. **/
    /**
     * The cartesian coordinate of the sub-box bound to this processor after running grid decomposition.
     */
    const int (&grid_coord)[DIMENSION_SIZE] = _grid_coord;

    /**
     * the rank ids of contiguous processors in space.
     */
    const _MPI_Rank (&rank_id_neighbours)[DIMENSION_SIZE][2] = _rank_id_neighbours;

  protected:
    void decomposition(const mpi_process _mpi_pro, MPI_Comm *_p_comm, const int _mpi_map_3d_sub_dim[DIMENSION_SIZE]);
    void decomposition_imp_v1(const mpi_process _mpi_pro, MPI_Comm *_p_comm);
    void decomposition_imp_v2(const mpi_process _mpi_pro, MPI_Comm *_p_comm,
                              const int _mpi_map_3d_sub_dim[DIMENSION_SIZE]);

  protected:
    int _grid_size[DIMENSION_SIZE] = {0};

    int _grid_coord[DIMENSION_SIZE] = {0, 0, 0};
    _MPI_Rank _rank_id_neighbours[DIMENSION_SIZE][2] = {{0, 0}, {0, 0}, {0, 0}};
  };
} // namespace comm

#endif // COMM_CARTESIAN_CONTEXT_H
