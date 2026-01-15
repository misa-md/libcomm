//
// Created by genshen on 2026/1/15.
//

#include "cart_3d_node_affinity.h"
#include "comm/types_define.h"

#include <mpi.h>
#include <stdexcept>

void comm::topology::LibComm_shared_Cart_3d_create(MPI_Comm old_comm, const int _dimensions,
                                                   const int node_dims[comm::DIMENSION_SIZE],
                                                   const int inner_node_dim[comm::DIMENSION_SIZE],
                                                   int _period[comm::DIMENSION_SIZE], bool _reorder,
                                                   comm::topology::Comm3dCart *cart_comm) {
  int world_rank;
  MPI_Comm_rank(old_comm, &world_rank);

  // split the node communication domain
  MPI_Comm node_comm;
  MPI_Comm_split_type(old_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);

  // get MPI  numbers and rank in the node
  int rank_in_node, rank_size_in_node;
  MPI_Comm_rank(node_comm, &rank_in_node);
  MPI_Comm_size(node_comm, &rank_size_in_node);

  if (rank_size_in_node != inner_node_dim[0] * inner_node_dim[1] * inner_node_dim[2]) {
    throw std::runtime_error("error: rank number in node not match.");
    MPI_Abort(MPI_COMM_WORLD, 2);
  }

  // separate out the representatives of the nodes (rank 0) form a comm
  // get the number and id of the root node.
  MPI_Comm _node_root_comm;
  MPI_Comm_split(old_comm, (rank_in_node == 0) ? 0 : MPI_UNDEFINED, world_rank, &_node_root_comm);

  // map the nodes into 3d grid.
  if (rank_in_node == 0) {
    MPI_Comm new_node_3d_comm;
    const int periods[comm::DIMENSION_SIZE] = {1, 1, 1};
    MPI_Cart_create(_node_root_comm, comm::DIMENSION_SIZE, node_dims, periods, true, &new_node_3d_comm);
    cart_comm->node_root_comm = new_node_3d_comm;
  } else {
    cart_comm->node_root_comm = MPI_COMM_NULL;
  }

  cart_comm->node_comm = node_comm;
}

void comm::topology::LibComm_shared_Cart_3d_coords(comm::topology::Comm3dCart *cart_comm,
                                                   const int node_dims[comm::DIMENSION_SIZE],
                                                   const int inner_node_dim[comm::DIMENSION_SIZE],
                                                   const int grid_size[comm::DIMENSION_SIZE],
                                                   int grid_coord[comm::DIMENSION_SIZE]) {

  int rank_in_node = 0;
  MPI_Comm_rank(cart_comm->node_comm, &rank_in_node);

  int node_coord[3] = {0, 0, 0};
  if (rank_in_node == 0) {
    int node_root_rank; // , node_root_size;
    MPI_Comm_rank(cart_comm->node_root_comm, &node_root_rank);
    MPI_Cart_coords(cart_comm->node_root_comm, node_root_rank, comm::DIMENSION_SIZE, node_coord);
  }

  MPI_Bcast(node_coord, comm::DIMENSION_SIZE, MPI_INT, 0, cart_comm->node_comm);

  // local coordinate inside the node
  const int local_coord[comm::DIMENSION_SIZE] = {rank_in_node % inner_node_dim[0],
                                           (rank_in_node / inner_node_dim[0]) % inner_node_dim[1],
                                           rank_in_node / (inner_node_dim[0] * inner_node_dim[1])};

  // global coord of each MPI rank
  for (int d = 0; d < 3; d++) {
    grid_coord[d] = node_coord[d] * inner_node_dim[d] + local_coord[d];
  }
}

void comm::topology::LibComm_shared_Cart_3d_shift(const int grid_size[comm::DIMENSION_SIZE], const int dimension,
                                                  const int disp, const int grid_coord[comm::DIMENSION_SIZE],
                                                  int *rank_id_neighbours_low, int *rank_id_neighbours_high) {
  const int gx = grid_size[0];
  const int gy = grid_size[1];
  const int gz = grid_size[2];

  const int nx = grid_coord[0];
  const int ny = grid_coord[1];
  const int nz = grid_coord[2];

  const int xlo = (nx - 1 + gx) % gx;
  const int xhi = (nx + 1 + gx) % gx;
  const int ylo = (ny - 1 + gy) % gy;
  const int yhi = (ny + 1 + gy) % gy;
  const int zlo = (nz - 1 + gz) % gz;
  const int zhi = (nz + 1 + gz) % gz;

  if (dimension == 0) {
    *rank_id_neighbours_low = comm::topology::Comm3dCart::linear_rank({gx, gy, gz}, {xlo, ny, nz});
    *rank_id_neighbours_high = comm::topology::Comm3dCart::linear_rank({gx, gy, gz}, {xhi, ny, nz});
  }
  if (dimension == 1) {
    *rank_id_neighbours_low = comm::topology::Comm3dCart::linear_rank({gx, gy, gz}, {nx, ylo, nz});
    *rank_id_neighbours_high = comm::topology::Comm3dCart::linear_rank({gx, gy, gz}, {nx, yhi, nz});
  }
  if (dimension == 2) {
    *rank_id_neighbours_low = comm::topology::Comm3dCart::linear_rank({gx, gy, gz}, {nx, ny, zlo});
    *rank_id_neighbours_high = comm::topology::Comm3dCart::linear_rank({gx, gy, gz}, {nx, ny, zhi});
  }
}

void comm::topology::LibComm_shared_Cart_3d_commit(MPI_Comm old_comm, MPI_Comm *new_comm,
                                                   const int grid_size[comm::DIMENSION_SIZE],
                                                   const int grid_coord[comm::DIMENSION_SIZE]) {
  if (grid_coord == nullptr || grid_size == nullptr || new_comm == nullptr) {
    throw std::runtime_error("error: MPI communicator or meta data is null.");
  }

  const int linear_rank = comm::topology::Comm3dCart::linear_rank({grid_size[0], grid_size[1], grid_size[2]},
                                                                  {grid_coord[0], grid_coord[1], grid_coord[2]});
  MPI_Comm_split(old_comm, 0, linear_rank, new_comm);
}
