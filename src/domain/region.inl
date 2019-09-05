
#include "region.hpp"

//
// Created by genshen on 2019-04-15.
//

template<typename T>
comm::Region<T>::Region() {}

template<typename T>
comm::Region<T>::Region(const T x_start, const T y_start, const T z_start,
                        const T x_end, const T y_end, const T z_end) {
    low[0] = x_start;
    low[1] = y_start;
    low[2] = z_start;
    high[0] = x_end;
    high[1] = y_end;
    high[2] = z_end;
}

template<typename T>
comm::Region<T> &comm::Region<T>::operator=(const comm::Region<T> &r) {
    memcpy(this->low, r.low, DIMENSION_SIZE * sizeof(T));
    memcpy(this->high, r.high, DIMENSION_SIZE * sizeof(T));
    return *this;
}
