//
// Created by genshen on 2026/1/15.
//

#ifndef COMM_CART_3D_NODE_AFFINITY_H
#define COMM_CART_3D_NODE_AFFINITY_H

#include "comm/types_define.h"

#include <array>
#include <mpi.h>

namespace comm {
  namespace topology {
    struct Comm3dCart {
      Comm3dCart() : node_root_comm(MPI_COMM_NULL), node_comm(MPI_COMM_NULL) {};

      ~Comm3dCart() {
        if (node_root_comm != MPI_COMM_NULL) {
          MPI_Comm_free(&node_root_comm);
        }
        if (node_comm != MPI_COMM_NULL) {
          MPI_Comm_free(&node_comm);
        }
        node_root_comm = MPI_COMM_NULL;
        node_comm = MPI_COMM_NULL;
      }

      MPI_Comm
          node_root_comm; // this comm contains all rank 0 on each node. only rank0 in the node has valid comm value.
      MPI_Comm node_comm; // this comm contains the ranks on a node. A node is a comm.

      /**
       * map 3d coordination of mpi process to 1d mpi rank.
       */
      static inline int linear_rank(const std::array<int, comm::DIMENSION_SIZE> grid_size,
                                    const std::array<int, comm::DIMENSION_SIZE> grid_coord) {
        const int gx = grid_size[0];
        const int gy = grid_size[1];
        const int gz = grid_size[2];

        const int nx = grid_coord[0];
        const int ny = grid_coord[1];
        const int nz = grid_coord[2];
        return nx + gx * (ny + gy * nz);
      }
    };

    /**
     * just like: MPI_Cart_create.
     * LibComm_Cart_3d_create_shared is an optimization version of MPI_Cart_create.
     * Note:
     * period[i] = 1.
     * reorder = true.
     */
    void LibComm_shared_Cart_3d_create(MPI_Comm old_comm, const int _dimensions,
                                       const int node_dims[comm::DIMENSION_SIZE],
                                       const int inner_node_dim[comm::DIMENSION_SIZE],
                                       int _period[comm::DIMENSION_SIZE], bool _reorder, Comm3dCart *cart_comm);

    /**
     * just like: MPI_Cart_coords.
     * It returns the new coords of current rank.
     */
    void LibComm_shared_Cart_3d_coords(comm::topology::Comm3dCart *cart_comm, const int node_dims[comm::DIMENSION_SIZE],
                                       const int inner_node_dim[comm::DIMENSION_SIZE],
                                       const int grid_size[comm::DIMENSION_SIZE], int grid_coord[comm::DIMENSION_SIZE]);

    void LibComm_shared_Cart_3d_shift(const int grid_size[comm::DIMENSION_SIZE], const int dimension, const int disp,
                                      const int grid_coord[comm::DIMENSION_SIZE], int *rank_id_neighbours_low,
                                      int *rank_id_neighbours_high);

    void LibComm_shared_Cart_3d_commit(MPI_Comm old_comm, MPI_Comm *new_comm, const int grid_size[comm::DIMENSION_SIZE],
                                       const int grid_coord[comm::DIMENSION_SIZE]);
  } // namespace topology
} // namespace comm
#endif // COMM_CART_3D_NODE_AFFINITY_H
