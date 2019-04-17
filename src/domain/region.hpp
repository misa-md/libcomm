//
// Created by genshen back to 2018-12-30.
//

#ifndef COMM_DOMAIN_REGION_HPP
#define COMM_DOMAIN_REGION_HPP

#include "../types_define.h"

namespace comm {
    template<typename T>
    struct Region {
        T low[DIMENSION_SIZE];
        T high[DIMENSION_SIZE];

        T &x_low = low[0];
        T &y_low = low[1];
        T &z_low = low[2];
        T &x_high = high[0];
        T &y_high = high[1];
        T &z_high = high[2];

        explicit Region();

        /**
         * initialize the region by given values.
         * \param x_start start position at x dimension.
         * \param y_start start position at y dimension.
         * \param z_start start position at z dimension.
         * \param x_end ending position at x dimension.
         * \param y_end ending position at y dimension.
         * \param z_end ending position at z dimension.
         */
        Region(const T x_start, const T y_start, const T z_start,
               const T x_end, const T y_end, const T z_end);

        /**
         * \param x,y,z a coordinate point.
         * \return return true if the coordinate \param x,y,z is in this region
         */
        inline bool isIn(const T x, const T y, const T z) const {
            return x > x_low && x < x_high && y > y_low && y < y_high && z > z_low && z < z_high;;
        }
    };

#include "region.inl"

}

#endif // COMM_DOMAIN_REGION_HPP