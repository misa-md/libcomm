//
// Created by genshen on 2019-03-15.
//

#ifndef COMM_COMM_H
#define COMM_COMM_H

#include "packer.h"
#include "types_define.h"

namespace comm {
    /**
     *
     * \brief communicate with its neighbour processes to exchange data.
     * This function communicate at x dimension, then y dimension, at last z dimension.
     * \tparam T type used in \param packer; the data type to be packed/unpacked.
     * \param packer packer pointer which does pack and unpack data.
     * \param processes the MPI rank and communicator in process.
     * \param data_type the MPI DataType used in communication.
     * \param neighbours_rank the rank id of all neighbour processes.
     * \param reversed the order of dimension loop will be z,y,x if \param reversed is true.
     * otherwise the loop order will be x,y,z. Default value is false (order: x,y,z)
     */
    template<typename T>
    void neiSendReceive(Packer<T> *packer,
                        const mpi_process processes,
                        const MPI_Datatype data_type,
                        const _MPI_Rank (&neighbours_rank)[DIMENSION_SIZE][2],
                        const bool reversed = false);

    /**
     *
     * \brief the similar function as above, but the communication DataType is in template param.
     * For some simple type like MPI_DOUBLE (constance, and can be determined at compiling time),
     * can use this function.
     * \tparam data_type the MPI DataType used in communication.
     */
    template<typename T, MPI_Datatype DT>
    inline void neiSendReceive(Packer<T> *packer, const mpi_process processes,
                               const _MPI_Rank (&neighbours_rank)[DIMENSION_SIZE][2],
                               const bool reversed = false) {
        neiSendReceive(packer, processes, DT, neighbours_rank, reversed);
    }

    /**
     * single side forwarding communication.
     * \tparam T type of data packer.
     * \param packer data packer.
     * \param send_dirs send directions of each dimension.
     * \param recv_dirs receive directions of each dimension.
     * \param processes communication domain.
     * \param data_type mpi data type to be exchanged.
     * \param neighbours_rank rank ids of neighbour ranks in each dimension.
     */
    template<typename T>
    void singleSideForwardComm(Packer<T> *packer,
                               const mpi_process processes,
                               const MPI_Datatype data_type,
                               const unsigned int send_dirs[DIMENSION_SIZE],
                               const unsigned int recv_dirs[DIMENSION_SIZE],
                               const _MPI_Rank neighbours_rank[DIMENSION_SIZE][2]);

};

#include "comm.inl"

#endif //COMM_COMM_H
