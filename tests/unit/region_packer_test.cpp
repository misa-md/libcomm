//
// Created by genshen on 2019/10/12.
//

#include <comm/comm.hpp>
#include <comm/region_packer.h>
#include <gtest/gtest.h>

class TestRegionPacker : public comm::RegionPacker<int, int> {
public:
  const unsigned long sendLength(const std::vector<comm::Region<int>> send_regions, const int dimension,
                                 const int direction) override {
    return 1;
  }

  void onSend(int buffer[], const std::vector<comm::Region<int>> send_regions, const unsigned long send_len,
              const int dimension, const int direction) override {
    buffer[0] = 0x100 + dimension;
  }

  void onReceive(int buffer[], const std::vector<comm::Region<int>> recv_regions, const unsigned long receive_len,
                 const int dimension, const int direction) override {
    data[dimension] = buffer[0];
  }

private:
  int data[comm::DIMENSION_SIZE] = {0, 0, 0};

  FRIEND_TEST(region_packer_comm_test, region_packer_test);
};

// @MPI
TEST(region_packer_comm_test, region_packer_test) {
  int ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  // Get the rank of the process
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  unsigned int my_rank_uint = my_rank;

  TestRegionPacker packer;
  comm::mpi_process pro = comm::mpi_process{my_rank, ranks, MPI_COMM_WORLD};
  comm::singleSideForwardComm(&packer, pro, MPI_INT, {}, {},
                              {my_rank_uint, my_rank_uint, my_rank_uint},  // send to itself
                              {my_rank_uint, my_rank_uint, my_rank_uint}); // receive by itself
  EXPECT_EQ(packer.data[0], 0x100);
  EXPECT_EQ(packer.data[1], 0x101);
  EXPECT_EQ(packer.data[2], 0x102);
}
