//
// Created by genshen on 2019-03-15.
//

#ifndef COMM_COMM_H
#define COMM_COMM_H

#include <array>
#include <mpi.h>
#include "packer.h"
#include "types_define.h"
#include "region_packer.h"

namespace comm {
    /**
     *
     * \brief communicate with its neighbour processes to exchange data.
     * This function communicate at x dimension, then y dimension, at last z dimension.
     * \tparam T type used in \param packer; the data type to be packed/unpacked.
     * \tparam F reversed dimension order to z->y->x (not x->y->z),
     *  the order of dimension loop will be z,y,x if \tparam F is true.
     *  otherwise the loop order will be x,y,z. Default value is false (order: x,y,z)
     * \param packer packer pointer which does pack and unpack data.
     * \param processes the MPI rank and communicator in process.
     * \param data_type the MPI DataType used in communication.
     * \param neighbours_rank the rank id of all neighbour processes.
     */
    template<typename T, bool F = false>
    void neiSendReceive(Packer<T> *packer,
                        const mpi_process processes,
                        const MPI_Datatype data_type,
                        const _MPI_Rank (&neighbours_rank)[DIMENSION_SIZE][2]);

    /**
     *
     * \brief the similar function as above, but the communication DataType is in template param.
     * For some simple type like MPI_DOUBLE (constance, and can be determined at compiling time),
     * can use this function.
     * \tparam F reversed dimension order to z->y->x (not x->y->z)
     * \tparam data_type the MPI DataType used in communication.
     * \param packer the data packer
     * \param processes mpi communication domain
     * \param neighbours_rank neighbor processes' rank id in each dimension and direction
     *
     */
    template<typename T, MPI_Datatype DT, bool F = false>
    inline void neiSendReceive(Packer<T> *packer, const mpi_process processes,
                               const _MPI_Rank (&neighbours_rank)[DIMENSION_SIZE][2]) {
        neiSendReceive<T, F>(packer, processes, DT, neighbours_rank);
    }

    /**
     * single side forwarding communication.
     * \tparam T type of data packer.
     * \tparam RT type of region type.
     * \tparam F reversed dimension order
     * \param packer data packer.
     * \param processes communication domain.
     * \param send_regions send regions used in communication in each dimension.
     * \param recv_regions receive regions used in communication in each dimension.
     * \param data_type mpi data type to be exchanged.
     * \param send_ranks MPI rand id in each dimension for sending communication.
     * \param recv_ranks MPI rand id in each dimension for receiving communication.
     */
    template<typename T, typename RT, bool F = false>
    void singleSideForwardComm(RegionPacker<T, RT> *packer,
                               const mpi_process processes,
                               const MPI_Datatype data_type, // todo return type in packer
                               const std::array<std::vector<comm::Region<RT>>, DIMENSION_SIZE> send_regions,
                               const std::array<std::vector<comm::Region<RT>>, DIMENSION_SIZE> recv_regions,
                               const std::array<unsigned int, DIMENSION_SIZE> send_ranks,
                               const std::array<unsigned int, DIMENSION_SIZE> recv_ranks);

}

#include "comm.inl"

#endif //COMM_COMM_H
