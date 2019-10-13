//
// Created by genshen on 2019/10/13.
//

#include <array>
#include <types_define.h>
#include <gtest/gtest.h>
#include <preset/sector_forwarding_direction.h>

// another implementation of get sending directions
std::array<unsigned int, comm::DIMENSION_SIZE> test_SSFDCommSendDirs_expected(const comm::_type_sector_id sector_id) {
    // for example, for sector id = Z_HIGH & Y_LOW & X_LOW = 0b100
    // it will return {0, 1, 1} = {X_HIGH,  Y_HIGH, Z_LOW}
    switch (sector_id) {
        case 0: // 0b000
            return {comm::DIR_HIGHER, comm::DIR_HIGHER, comm::DIR_HIGHER};
        case 1: // 0b001
            return {comm::DIR_LOWER, comm::DIR_HIGHER, comm::DIR_HIGHER};
        case 2: // 0b010
            return {comm::DIR_HIGHER, comm::DIR_LOWER, comm::DIR_HIGHER};
        case 3: // 0b011
            return {comm::DIR_LOWER, comm::DIR_LOWER, comm::DIR_HIGHER};
        case 4: // 0b100
            return {comm::DIR_HIGHER, comm::DIR_HIGHER, comm::DIR_LOWER};
        case 5: // 0b101
            return {comm::DIR_LOWER, comm::DIR_HIGHER, comm::DIR_LOWER};
        case 6: // 0b110
            return {comm::DIR_HIGHER, comm::DIR_LOWER, comm::DIR_LOWER};
        case 7: // 0x111
            return {comm::DIR_LOWER, comm::DIR_LOWER, comm::DIR_LOWER};
        default:
            assert(false);
            return {};
    }
}

// another implementation of get receiving directions
std::array<unsigned int, comm::DIMENSION_SIZE> test_SSFDCommRecvDirs_expected(const comm::_type_sector_id sector_id) {
    switch (sector_id) {
        case 0:
            return {comm::DIR_LOWER, comm::DIR_LOWER, comm::DIR_LOWER};
        case 1:
            return {comm::DIR_HIGHER, comm::DIR_LOWER, comm::DIR_LOWER};
        case 2:
            return {comm::DIR_LOWER, comm::DIR_HIGHER, comm::DIR_LOWER};
        case 3:
            return {comm::DIR_HIGHER, comm::DIR_HIGHER, comm::DIR_LOWER};
        case 4:
            return {comm::DIR_LOWER, comm::DIR_LOWER, comm::DIR_HIGHER};
        case 5:
            return {comm::DIR_HIGHER, comm::DIR_LOWER, comm::DIR_HIGHER};
        case 6:
            return {comm::DIR_LOWER, comm::DIR_HIGHER, comm::DIR_HIGHER};
        case 7:
            return {comm::DIR_HIGHER, comm::DIR_HIGHER, comm::DIR_HIGHER};
        default:
            assert(false);
            return {};
    }
}

TEST(sector_fd_dir_test_send_dirs, sector_fd_dir_test) {
    for (comm::_type_sector_id s = 0; s < 8; s++) {
        const std::array<unsigned int, comm::DIMENSION_SIZE> send_dirs = comm::ssfdCommSendDirs(s);
        const std::array<unsigned int, comm::DIMENSION_SIZE> send_dirs_expected = test_SSFDCommSendDirs_expected(s);
        for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
            EXPECT_EQ(send_dirs[d], send_dirs_expected[d]);
        }
    }
}

TEST(sector_fd_dir_test_recv_dirs, sector_fd_dir_test) {
    for (comm::_type_sector_id s = 0; s < 8; s++) {
        const std::array<unsigned int, comm::DIMENSION_SIZE> recv_dirs = comm::ssfdCommRecvDirs(s);
        const std::array<unsigned int, comm::DIMENSION_SIZE> recv_dirs__expected = test_SSFDCommRecvDirs_expected(s);
        for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
            EXPECT_EQ(recv_dirs[d], recv_dirs__expected[d]);
        }
    }
}
