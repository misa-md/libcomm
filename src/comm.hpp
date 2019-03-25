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
     * \tparam T type used in \param packer; the data type to be packed/unpacked.
     * \param packer packer pointer which does pack and unpack data.
     * \param processes the MPI rank and communicator in process.
     * \param data_type the MPI DataType used in communication.
     * \param neighbours_rank the rank id of all neighbour processes.
     */
    template<typename T>
    void neiSendReceive(Packer<T> *packer,
                        const mpi_process processes,
                        const MPI_Datatype data_type,
                        const _MPI_Rank (&neighbours_rank)[DIMENSION][2]);

    /**
     *
     * \brief the similar function as above, but the communication DataType is in template param.
     * For some simple type like MPI_DOUBLE (constance, and can be determined at compiling time),
     * can use this function.
     * \tparam data_type the MPI DataType used in communication.
     */
    template<typename T, MPI_Datatype DT>
    inline void neiSendReceive(Packer<T> *packer, const mpi_process processes,
                               const _MPI_Rank (&neighbours_rank)[DIMENSION][2]) {
        neiSendReceive(packer, processes, DT, neighbours_rank);
    }
};

#include "comm.inl"

#endif //COMM_COMM_H
