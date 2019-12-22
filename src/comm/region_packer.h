//
// Created by genshen on 2019/10/12.
//

#ifndef COMM_REGION_PACKER_H
#define COMM_REGION_PACKER_H

#include <vector>
#include "domain/region.hpp"

namespace comm {
    /**
     * \tparam T type of data
     * \tparam RT type of region type
     * \brief abstract class for packing data by regions.
     */
    template<typename T, typename RT>
    class RegionPacker {
    public:
        typedef T pack_date_type;
        typedef RT pack_region_type;

        /**
         * count the length to be send to neighbour process.
         * \param send_regions regions to be send in communication.
         * \param dimension 0,1,2. the id of dimensions.
         * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
         * \return
         */
        virtual const unsigned long sendLength(const std::vector<comm::Region<RT>> send_regions,
                                               const int dimension, const int direction) = 0;

        /**
         * This function will be called before sending data.
         * \param buffer the empty buffer to be send to its neighbor process, you can fill the buffer here.
         * \param send_regions regions to be send in communication.
         * \param send_size the length to be send.
         * \param dimension 0,1,2. the id of dimensions.
         * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
         */
        virtual void onSend(T buffer[], const std::vector<comm::Region<RT>> send_regions,
                            const unsigned long send_len, const int dimension, const int direction) = 0;

        /**
         * This function will be called after the receive finished.
         * We can unpack data here.
         * \param buffer the buffer of received data.
         * \param recv_regions regions to be received in communication.
         * \param receive_len the length of received data.
         * \param dimension 0,1,2. the id of dimensions.
         * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
         */
        virtual void onReceive(T buffer[], const std::vector<comm::Region<RT>> recv_regions,
                               const unsigned long receive_len, const int dimension, const int direction) = 0;

        /**
         * \brief this function will be called after all communication finished.
         */
        virtual void onFinish() {};
    };
}

#endif //COMM_REGION_PACKER_H
