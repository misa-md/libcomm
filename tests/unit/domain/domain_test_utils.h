//
// Created by genshen on 5/12/18.
//

#ifndef CRYSTAL_MD_DOMAIN_TEST_UTILS_H
#define CRYSTAL_MD_DOMAIN_TEST_UTILS_H

#include <cstdint>
#include <domain/domain.h>

comm::Domain *getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius);

#endif //CRYSTAL_MD_DOMAIN_TEST_UTILS_H
