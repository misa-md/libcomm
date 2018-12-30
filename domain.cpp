#include <logs/logs.h>
#include "domain.h"
#include "utils/mpi_domain.h"
#include "utils/mpi_data_types.h"
#include "pack/pack.h"


Domain::Domain(const std::array<u_int64_t, DIMENSION> _phase_space,
               const double _lattice_const, const double _cutoff_radius_factor)
        : lattice_const(_lattice_const), cutoff_radius_factor(_cutoff_radius_factor),
          phase_space(_phase_space),
//      todo _grid_size(0),
//      todo _meas_global_length(0.0),
        /**initialize following references */
          meas_global_length(_meas_global_length),
          grid_size(_grid_size),
          meas_global_box_coord_region(_meas_global_box_coord_region),
          grid_coord_sub_box(_grid_coord_sub_box),
          rank_id_neighbours(_rank_id_neighbours),
          meas_sub_box_region(_meas_sub_box_region),
          meas_ghost_length(_meas_ghost_length),
          meas_ghost_region(_meas_ghost_region),
          lattice_size_sub_box(_lattice_sub_box_size),
          lattice_size_ghost_extended(_lattice_size_ghost_extended),
          lattice_size_ghost(_lattice_size_ghost),
          lattice_coord_sub_box_region(_lattice_coord_sub_box_region),
          lattice_coord_ghost_region(_lattice_coord_ghost_region),
          local_sub_box_lattice_coord_region(_local_sub_box_lattice_coord_region),
          local_ghost_lattice_coord_region(_local_ghost_lattice_coord_region) {}

Domain::Builder &Domain::Builder::setPhaseSpace(const int64_t phaseSpace[DIMENSION]) {
    for (int i = 0; i < DIMENSION; i++) {
        _phase_space[i] = phaseSpace[i]; // todo type not match
    }
    return *this;
}

Domain::Builder &Domain::Builder::setLatticeConst(const double latticeConst) {
    _lattice_const = latticeConst;
    return *this;
}

Domain::Builder &Domain::Builder::setCutoffRadius(const double cutoff_radius_factor) {
    _cutoff_radius_factor = cutoff_radius_factor;
    return *this;
}

void Domain::Builder::decomposition(Domain &domain) {
    // Assume N can be decomposed as N = N_x * N_y * N_z,
    // then we have: _grid_size[0] = N_x, _grid_size[1] = N_y, _grid_size[1] = N_z.
    // Fill in the _grid_size array such that the product of _grid_size[i] for i=0 to DIMENSION-1 equals N.
    MPI_Dims_create(MPIDomain::sim_processor.all_ranks, DIMENSION,
                    domain._grid_size); // fixme origin code: (int *) &domain._grid_size
    kiwi::logs::i(MASTER_PROCESSOR, "decomposition", "MPI grid dimensions: {0},{1},{2}\n",
                  domain._grid_size[0], domain._grid_size[1], domain._grid_size[2]);

    int period[DIMENSION];
    // 3维拓扑
    for (int d = 0; d < DIMENSION; d++) {
        period[d] = 1;
    }
    // sort the processors to fit 3D cartesian topology.
    // the rank id may change.
    MPI_Comm _comm;
    int _debug_old_rank = MPIDomain::sim_processor.own_rank;
    MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, domain._grid_size, period, true, &_comm);
    kiwi::mpiUtils::onGlobalCommChanged(_comm);
    MPIDomain::sim_processor = kiwi::mpiUtils::global_process; // set new domain.

    // get cartesian coordinate of current processor.
    MPI_Cart_coords(MPIDomain::sim_processor.comm, MPIDomain::sim_processor.own_rank, DIMENSION,
                    domain._grid_coord_sub_box);
    kiwi::logs::d("decomposition", "old_rank_id: {0}, MPI coordinate of current process: x:{1},y{2},z{3}\n",
                  _debug_old_rank, domain._grid_coord_sub_box[0], domain._grid_coord_sub_box[1],
                  domain._grid_coord_sub_box[2]);

    // get the rank ids of contiguous processors of current processor.
    for (int d = 0; d < DIMENSION; d++) {
        MPI_Cart_shift(_comm, d, 1, &domain._rank_id_neighbours[d][LOWER], &domain._rank_id_neighbours[d][HIGHER]);
    }
}

void Domain::Builder::createGlobalDomain(Domain &domain) {
    for (int d = 0; d < DIMENSION; d++) {
        //phaseSpace个单位长度(单位长度即latticeconst)
        domain._meas_global_length[d] = _phase_space[d] * _lattice_const;
        domain._meas_global_box_coord_region.low[d] = 0; // lower bounding is set to 0 by default.
        domain._meas_global_box_coord_region.high[d] = domain._meas_global_length[d];
    }
}

void Domain::Builder::createSubBoxDomain(Domain &domain) {
    // calculate measured length in each dimension.
    for (int d = 0; d < DIMENSION; d++) {
        // the lower and upper bounding of current sub-box.
        domain._meas_sub_box_region.low[d] = domain._meas_global_box_coord_region.low[d] +
                                             domain._grid_coord_sub_box[d] *
                                             (domain._meas_global_length[d] / domain._grid_size[d]);
        domain._meas_sub_box_region.high[d] = domain._meas_global_box_coord_region.low[d] +
                                              (domain._grid_coord_sub_box[d] + 1) *
                                              (domain._meas_global_length[d] / domain._grid_size[d]);

        domain._meas_ghost_length[d] = _cutoff_radius_factor * _lattice_const; // ghost length todo

        domain._meas_ghost_region.low[d] = domain._meas_sub_box_region.low[d] - domain._meas_ghost_length[d];
        domain._meas_ghost_region.high[d] = domain._meas_sub_box_region.high[d] + domain._meas_ghost_length[d];
    }

    // set lattice size of local sub-box.
    for (int d = 0; d < DIMENSION; d++) {
        domain._lattice_sub_box_size[d] = (domain._grid_coord_sub_box[d] + 1) * _phase_space[d] / domain._grid_size[d] -
                                          (domain._grid_coord_sub_box[d]) * _phase_space[d] / domain._grid_size[d];
    }
    domain._lattice_sub_box_size[0] *= 2;

    /*
    nghostx = p_domain->getSubBoxLatticeSize(0) + 2 * 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    nghosty = p_domain->getSubBoxLatticeSize(1) + 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    nghostz = p_domain->getSubBoxLatticeSize(2) + 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    */
    // set ghost lattice size.
    for (int d = 0; d < DIMENSION; d++) {
        // i * ceil(x) >= ceil(i*x) for all x ∈ R and i ∈ Z
        domain._lattice_size_ghost[d] = (d == 0) ? 2 * ceil(_cutoff_radius_factor) : ceil(_cutoff_radius_factor);
        domain._lattice_size_ghost_extended[d] = domain._lattice_sub_box_size[d] + 2 * domain._lattice_size_ghost[d];
    }

    // set lattice coordinate boundary in global and local coordinate system(GCS and LCS).
    buildSubBoxDomainGCS(domain);
    buildSubBoxDomainLCS(domain);
}

void Domain::Builder::buildSubBoxDomainGCS(Domain &domain) {
    // set lattice coordinate boundary of sub-box.
    for (int d = 0; d < DIMENSION; d++) {
        // floor equals to "/" if all operation number >=0.
        domain._lattice_coord_sub_box_region.low[d] = domain._grid_coord_sub_box[d] * _phase_space[d] /
                                                      domain._grid_size[d]; // todo set measure coord = lower*lattice_const.
        domain._lattice_coord_sub_box_region.high[d] =
                (domain._grid_coord_sub_box[d] + 1) * _phase_space[d] / domain._grid_size[d];
    }
    domain._lattice_coord_sub_box_region.x_low *= 2;
    domain._lattice_coord_sub_box_region.x_high *= 2;

    /*
   loghostx = p_domain->getSubBoxLatticeCoordLower(0) - 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
   loghosty = p_domain->getSubBoxLatticeCoordLower(1) - ( ceil( cutoffRadius / _latticeconst ) + 1 );
   loghostz = p_domain->getGlobalSubBoxLatticeCoordLower(2) - ( ceil( cutoffRadius / _latticeconst ) + 1 );
   */
    // set lattice coordinate boundary for ghost.
    for (int d = 0; d < DIMENSION; d++) {
        // todo too integer minus, cut too many??
        domain._lattice_coord_ghost_region.low[d] = domain._lattice_coord_sub_box_region.low[d] -
                                                    domain._lattice_size_ghost[d];
        domain._lattice_coord_ghost_region.high[d] = domain._lattice_coord_sub_box_region.high[d] +
                                                     domain._lattice_size_ghost[d];
    }
}

void Domain::Builder::buildSubBoxDomainLCS(Domain &domain) {
    for (int d = 0; d < DIMENSION; d++) {
        domain._local_ghost_lattice_coord_region.low[d] = 0;
        domain._local_ghost_lattice_coord_region.high[d] = domain._lattice_size_ghost_extended[d];
        domain._local_sub_box_lattice_coord_region.low[d] = domain._lattice_size_ghost[d];
        domain._local_sub_box_lattice_coord_region.high[d] =
                domain._lattice_size_ghost[d] + domain._lattice_sub_box_size[d];
    }
}

Domain *Domain::Builder::build() {
    Domain *p_domain = new Domain(_phase_space, _lattice_const, _cutoff_radius_factor); // todo
    decomposition(*p_domain);
    createGlobalDomain(*p_domain);
    createSubBoxDomain(*p_domain);
    return p_domain;
}
