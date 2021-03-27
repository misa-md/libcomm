//
// Created by genshen on 2019/10/13.
//

#ifndef COMM_SECTOR_FORWARDING_DIRECTION_H
#define COMM_SECTOR_FORWARDING_DIRECTION_H

#include <array>

#include "comm/types_define.h"

namespace comm {
  /**
   * \brief the send directions (can only be comm::DIR_LOW and comm::DIR_HIGH)
   * of single side forward communication in x,y,z dimensions for syncing ghost regions.
   * \param sector_id sector id.
   * \return send directions of each dimension.
   */
  inline std::array<unsigned int, comm::DIMENSION_SIZE> ssfdCommSendDirs(const _type_sector_id sector_id) {
    // for example, for sector id = Z_HIGH & Y_LOW & X_LOW = 0b100
    // it will return {1, 1, 0} = {X_HIGH,  Y_HIGH, Z_LOW}
    return {
        (0x7u - sector_id) & 0x1u,
        ((0x7u - sector_id) >> 1u) & 0x1u,
        ((0x7u - sector_id) >> 2u) & 0x1u,
    };
  }

  /**
   * \brief the receive directions (can only be comm::DIR_LOW and comm::DIR_HIGH)
   * of single side forward communication in x,y,z dimensions for syncing ghost regions.
   * \param sector_id sector id.
   * \return receive directions of each dimension.
   */
  inline std::array<unsigned int, comm::DIMENSION_SIZE> ssfdCommRecvDirs(const _type_sector_id sector_id) {
    // for example, for sector id = Z_HIGH & Y_LOW & X_LOW = 0b100
    // it will return {0, 0, 1} = {X_LOW,  Y_LOW, Z_HIGH}
    return {
        sector_id & 0x1u,
        (sector_id >> 1u) & 0x1u,
        (sector_id >> 2u) & 0x1u,
    };
  }
}; // namespace comm

#endif // COMM_SECTOR_FORWARDING_DIRECTION_H
