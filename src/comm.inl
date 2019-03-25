//
// Created by genshen on 2019-03-15.
//

#include <cassert>

template<typename T>
void comm::neiSendReceive(Packer<T> *packer,
                          const mpi_process processes,
                          const MPI_Datatype data_type,
                          const _MPI_Rank (&neighbours_rank)[DIMENSION][2]) {
    unsigned int num_send[DIMENSION][2];
    unsigned int num_receive[DIMENSION][2];
    T *send_buff[2];
    T *receive_buff[2];

    MPI_Status status;
    MPI_Status send_statuses[DIMENSION][2];
    MPI_Status recv_statuses[DIMENSION][2];
    MPI_Request send_requests[DIMENSION][2];
    MPI_Request recv_requests[DIMENSION][2];

    int direction;
    for (int d = 0; d < DIMENSION; d++) {
        for (direction = LOWER; direction <= HIGHER; direction++) {
            num_send[d][direction] = packer->sendLength(d, direction);
            assert(num_send[d][direction] >= 0); // todo remove
            send_buff[direction] = new T[num_send[d][direction]];
            packer->onSend(send_buff[direction], num_send[d][direction], d, direction);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            unsigned int &numsend = num_send[d][direction];
            int numrecv = 0;

            MPI_Isend(send_buff[direction], numsend, data_type,
                      neighbours_rank[d][direction], 99,
                      processes.comm, &send_requests[d][direction]);
            // test the status of neighbor process.
            MPI_Probe(neighbours_rank[d][(direction + 1) % 2], 99, processes.comm, &status);
            // test the data length to be received.
            MPI_Get_count(&status, data_type, &numrecv);
            // initialize receive buffer via receiving size.
            // the receiving length is get bt MPI_Probe from its neighbour process.
            assert(numrecv >= 0);  // todo remove
            receive_buff[direction] = new T[numrecv];
            num_receive[d][direction] = numrecv;
            MPI_Irecv(receive_buff[direction], numrecv, data_type,
                      neighbours_rank[d][(direction + 1) % 2], 99,
                      processes.comm, &recv_requests[d][direction]);
        }
        for (direction = LOWER; direction <= HIGHER; direction++) {
            MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]);
            MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]);
            packer->onReceive(receive_buff[direction], num_receive[d][direction], d, direction);
            // release buffer
            delete[] send_buff[direction];
            delete[] receive_buff[direction];
        }
    }
    // all finished
    packer->onFinish();
}
