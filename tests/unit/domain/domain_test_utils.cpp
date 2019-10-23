//
// Created by genshen on 5/12/18.
//

#include <mpi.h>
#include <comm/domain/bcc_domain.h>
#include "domain_test_utils.h"

comm::Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    MPI_Comm mpi_comm; // new mpi comm
    comm::mpi_process m_process{};
    m_process.comm = MPI_COMM_WORLD;
    MPI_Comm_size(m_process.comm, &m_process.all_ranks);
    MPI_Comm_rank(m_process.comm, &m_process.own_rank);

    comm::Domain *p_domain = comm::Domain::Builder()
            .setComm(m_process, &mpi_comm)
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
    return p_domain;
}

comm::Domain *getDomainInstance(int64_t space[3], const int ghost_size,
                                double lattice_const, double cutoff_radius) {
    MPI_Comm mpi_comm; // new mpi comm
    comm::mpi_process m_process{};
    m_process.comm = MPI_COMM_WORLD;
    MPI_Comm_size(m_process.comm, &m_process.all_ranks);
    MPI_Comm_rank(m_process.comm, &m_process.own_rank);

    comm::Domain *p_domain = comm::Domain::Builder()
            .setComm(m_process, &mpi_comm)
            .setPhaseSpace(space)
            .setGhostSize(ghost_size)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
    return p_domain;
}

// creat bcc domain.
comm::BccDomain *getBccDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    MPI_Comm mpi_comm; // new mpi comm
    comm::mpi_process m_process{};
    m_process.comm = MPI_COMM_WORLD;
    MPI_Comm_size(m_process.comm, &m_process.all_ranks);
    MPI_Comm_rank(m_process.comm, &m_process.own_rank);

    comm::BccDomain *p_domain = comm::BccDomain::Builder()
            .setComm(m_process, &mpi_comm)
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
    return p_domain;
}
