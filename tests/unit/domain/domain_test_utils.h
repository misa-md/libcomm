//
// Created by genshen on 5/12/18.
//

#ifndef CRYSTAL_MD_DOMAIN_TEST_UTILS_H
#define CRYSTAL_MD_DOMAIN_TEST_UTILS_H

#include <comm/domain/bcc_domain.h>
#include <comm/domain/domain.h>
#include <cstdint>

// create normal domain
comm::Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius);

comm::Domain *getDomainInstance(int64_t space[3], const int ghost_size, double lattice_const, double cutoff_radius);

// create bcc domain
comm::BccDomain *getBccDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius);

#endif // CRYSTAL_MD_DOMAIN_TEST_UTILS_H
