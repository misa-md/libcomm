//
// Created by genshen back to 2018-12-30.
//

#ifndef COMM_DOMAIN_REGION_HPP
#define COMM_DOMAIN_REGION_HPP

#include "comm/types_define.h"

namespace comm {
  template <typename T> struct Region {
    T x_low;
    T y_low;
    T z_low;
    T x_high;
    T y_high;
    T z_high;

    T *low = &x_low;

    T *high = &x_high;

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
    Region(const T x_start, const T y_start, const T z_start, const T x_end, const T y_end, const T z_end);

    Region(const Region &region); // copy ctor

    Region &operator=(const Region &arr);

    inline T *data() { return &x_low; }

    bool operator==(const Region &x) const {
      return x_low == x.x_low && y_low == x.y_low && z_low == x.z_low && x_high == x.x_high && y_high == x.y_high &&
             z_high == x.z_high;
    }

    bool operator!=(const Region &x) const { return !(*this == x); }

    /**
     * get the volume of region
     * \return
     */
    inline T volume() { return (x_high - x_low) * (y_high - y_low) * (z_high - z_low); }

    /**
     * Given a point (x,y,z), if x belongs to set [x_low, x_high), y belongs to set [y_low, y_high),
     * and z belongs to set [z_low, z_high),
     * this function will return true, otherwise, false will be returned.
     *
     * \param x,y,z a coordinate point.
     * \return return true if the coordinate \param x,y,z is in this region
     */
    inline bool isIn(const T x, const T y, const T z) const {
      return x >= x_low && x < x_high && y >= y_low && y < y_high && z >= z_low && z < z_high;
      ;
    }
  };
} // namespace comm

#include "region.inl"

#endif // COMM_DOMAIN_REGION_HPP