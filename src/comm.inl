//
// Created by genshen on 2019-03-15.
//

#include <cassert>
#include "comm.hpp"

template<typename T>
void comm::neiSendReceive(Packer<T> *packer,
                          const mpi_process processes,
                          const MPI_Datatype data_type,
                          const _MPI_Rank (&neighbours_rank)[DIMENSION][2],
                          const bool reversed) {
    unsigned int num_send[DIMENSION][2];
    unsigned int num_receive[DIMENSION][2];
    T *send_buff[2];
    T *receive_buff[2];

    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    // if dimension loop is reversed, the loop order will be z,y,x.
    // otherwise the loop order will be x,y,z.
    for (int d = 0; d < DIMENSION; d++) {
        int dimension;
        if (reversed) {
            dimension = DIMENSION - 1 - d;
        } else {
            dimension = d;
        }

        for (int direction = LOWER; direction <= HIGHER; direction++) {
            num_send[dimension][direction] = packer->sendLength(dimension, direction);
            assert(num_send[dimension][direction] >= 0); // todo remove
            send_buff[direction] = new T[num_send[dimension][direction]];
            packer->onSend(send_buff[direction], num_send[dimension][direction], dimension, direction);
        }
        for (int direction = LOWER; direction <= HIGHER; direction++) {
            unsigned int &numsend = num_send[dimension][direction];
            int numrecv = 0;

            MPI_Isend(send_buff[direction], numsend, data_type,
                      neighbours_rank[dimension][direction], 99,
                      processes.comm, &send_requests[dimension][direction]);
            // test the status of neighbor process.
            MPI_Probe(neighbours_rank[dimension][(direction + 1) % 2], 99, processes.comm, &status);
            // test the data length to be received.
            MPI_Get_count(&status, data_type, &numrecv);
            // initialize receive buffer via receiving size.
            // the receiving length is get bt MPI_Probe from its neighbour process.
            assert(numrecv >= 0);  // todo remove
            receive_buff[direction] = new T[numrecv];
            num_receive[dimension][direction] = numrecv;
            MPI_Irecv(receive_buff[direction], numrecv, data_type,
                      neighbours_rank[dimension][(direction + 1) % 2], 99,
                      processes.comm, &recv_requests[dimension][direction]);
        }
        for (int direction = LOWER; direction <= HIGHER; direction++) {
            MPI_Wait(&send_requests[dimension][direction], &send_statuses[dimension][direction]);
            MPI_Wait(&recv_requests[dimension][direction], &recv_statuses[dimension][direction]);
            packer->onReceive(receive_buff[direction], num_receive[dimension][direction], dimension, direction);
            // release buffer
            delete[] send_buff[direction];
            delete[] receive_buff[direction];
        }
    }
    // all finished
    packer->onFinish();
}
