//
// Created by genshen on 2019-03-15.
//

#ifndef COMM_TYPES_DEFINE_H
#define COMM_TYPES_DEFINE_H

#include <mpi.h>

// some value definition here.
#define DIMENSION 3

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

#define LOWER  0
#define HIGHER 1

namespace comm {
    const static int COMM_MASTER = 0;

    typedef int _type_lattice_size;
    typedef _type_lattice_size _type_lattice_coord;

    typedef double _type_atom_mass;
    typedef double _type_atom_location;
    typedef double _type_atom_velocity;
    typedef double _type_atom_force;
    typedef double _type_atom_rho;
    typedef double _type_atom_df;

    typedef int _MPI_Rank;
    struct mpi_process {
        _MPI_Rank own_rank, all_ranks;
        MPI_Comm comm;
    };
}

#endif //COMM_TYPES_DEFINE_H
