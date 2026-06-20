//
// Created by genshen on 2026/6/20.
//

#include <stdexcept>
#include <string>

#include "cartesian_context.h"

#include "comm/topology/cart_3d_node_affinity.h"

void comm::CartesianContext::decomposition(const mpi_process _mpi_pro, MPI_Comm *_p_comm,
                                           const int _mpi_map_3d_sub_dim[DIMENSION_SIZE]) {
  if (_mpi_map_3d_sub_dim[0] == 1 && _mpi_map_3d_sub_dim[1] == 1 && _mpi_map_3d_sub_dim[2] == 1) {
    decomposition_imp_v1(_mpi_pro, _p_comm);
  } else {
    // note1: for inner node 1x1x1, it can also call imp_v2.
    // but for the given mpi process, it may produce different mpi rank when comparing with imp_v1.
    decomposition_imp_v2(_mpi_pro, _p_comm, _mpi_map_3d_sub_dim);
  }
}

void comm::CartesianContext::decomposition_imp_v1(const mpi_process _mpi_pro, MPI_Comm *_p_comm) {
  // Assume N can be decomposed as N = N_x * N_y * N_z,
  // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
  // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION_SIZE-1 equals N.
  MPI_Dims_create(_mpi_pro.all_ranks, DIMENSION_SIZE, this->_grid_size);

  int period[DIMENSION_SIZE];
  // 3D topology
  for (int &d : period) {
    d = 1;
  }
  // sort the processors to fit 3D cartesian topology.
  // the rank id may change.
  MPI_Cart_create(_mpi_pro.comm, DIMENSION_SIZE, this->_grid_size, period, true, _p_comm);

  comm::_MPI_Rank new_rank;
  MPI_Comm_rank(*_p_comm, &new_rank);
  // get cartesian coordinate of current processor.
  MPI_Cart_coords(*_p_comm, new_rank, DIMENSION_SIZE, this->_grid_coord);

  // get the rank ids of contiguous processors of current processor.
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    MPI_Cart_shift(*_p_comm, d, 1, &(this->_rank_id_neighbours[d][DIR_LOWER]),
                   &(this->_rank_id_neighbours[d][DIR_HIGHER]));
  }
}

void comm::CartesianContext::decomposition_imp_v2(const mpi_process _mpi_pro, MPI_Comm *_p_comm,
                                                  const int _mpi_map_3d_sub_dim[DIMENSION_SIZE]) {
  // Assume N can be decomposed as N = N_x * N_y * N_z,
  // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
  // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION_SIZE-1 equals N.
  MPI_Dims_create(_mpi_pro.all_ranks, DIMENSION_SIZE, this->_grid_size);

  int node_dims[DIMENSION_SIZE];
  for (int d = 0; d < 3; d++) {
    if ((this->_grid_size[d] % _mpi_map_3d_sub_dim[d] != 0) || (this->_grid_size[d] < _mpi_map_3d_sub_dim[d])) {
      throw std::runtime_error("error: MPI map 3D sub-dimension does not match grid size in dimension " +
                               std::to_string(d) + ".");
      MPI_Abort(MPI_COMM_WORLD, -1);
    } else {
      node_dims[d] = this->_grid_size[d] / _mpi_map_3d_sub_dim[d];
    }
  }

  int period[DIMENSION_SIZE];
  // 3d with period boundary
  for (int &d : period) {
    d = 1;
  }

  comm::topology::Comm3dCart node_shared_cart_comm;
  // map node and ranks in a node into 3d topology
  comm::topology::LibComm_shared_Cart_3d_create(_mpi_pro.comm, DIMENSION_SIZE, node_dims, _mpi_map_3d_sub_dim, period,
                                                true, &node_shared_cart_comm);

  // get cartesian coordinate of current processor.
  comm::topology::LibComm_shared_Cart_3d_coords(&node_shared_cart_comm, node_dims, _mpi_map_3d_sub_dim,
                                                this->_grid_size, this->_grid_coord);

  // get rank is of neighbor mpi rank
  for (int d = 0; d < DIMENSION_SIZE; d++) {
    comm::topology::LibComm_shared_Cart_3d_shift(this->_grid_size, d, 1, this->_grid_coord,
                                                 &(this->_rank_id_neighbours[d][DIR_LOWER]),
                                                 &(this->_rank_id_neighbours[d][DIR_HIGHER]));
  }

  comm::topology::LibComm_shared_Cart_3d_commit(_mpi_pro.comm, _p_comm, this->_grid_size, this->_grid_coord);

// #define LIBCOMM_DEBUG
#ifdef LIBCOMM_DEBUG
  int world_rank, world_size;
  MPI_Comm_rank(_mpi_pro.comm, &world_rank);
  MPI_Comm_size(_mpi_pro.comm, &world_size);

  // get process numbers and quantities within the node
  int node_id = 0, rank_in_node, node_size;
  MPI_Comm_rank(node_shared_cart_comm.node_comm, &rank_in_node);
  MPI_Comm_size(node_shared_cart_comm.node_comm, &node_size);
  if (rank_in_node == 0) {
    MPI_Comm_rank(node_shared_cart_comm.node_root_comm, &node_id);
  }
  MPI_Bcast(&node_id, 1, MPI_INT, 0, node_shared_cart_comm.node_comm);

  int new_rank;
  int new_size;
  MPI_Comm_rank(*_p_comm, &new_rank);
  MPI_Comm_size(*_p_comm, &new_size);
  MPI_Barrier(*_p_comm);
  for (int r = 0; r < world_size; r++) {
    if (new_rank == r) {
      char processor_name[MPI_MAX_PROCESSOR_NAME];
      int name_len;
      MPI_Get_processor_name(processor_name, &name_len);
      printf("[new-rank %2d | node %2d#%s] [total_size %d | inner_size %d] coord=(%d,%d,%d) "
             "neigh[x-%d x+%d y-%d y+%d z-%d z+%d] grid_size=(%d,%d,%d) inner_grid=(%d,%d,%d)\n",
             new_rank, node_id, processor_name, new_size, node_size, this->_grid_coord[0], this->_grid_coord[1],
             this->_grid_coord[2], this->_rank_id_neighbours[0][DIR_LOWER], this->_rank_id_neighbours[0][DIR_HIGHER],
             this->_rank_id_neighbours[1][DIR_LOWER], this->_rank_id_neighbours[1][DIR_HIGHER],
             this->_rank_id_neighbours[2][DIR_LOWER], this->_rank_id_neighbours[2][DIR_HIGHER], this->_grid_size[0],
             this->_grid_size[1], this->_grid_size[2], _mpi_map_3d_sub_dim[0], _mpi_map_3d_sub_dim[1],
             _mpi_map_3d_sub_dim[2]);
      fflush(stdout);
    }
    MPI_Barrier(*_p_comm);
  }
#endif
}
