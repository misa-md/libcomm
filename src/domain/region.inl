//
// Created by genshen on 2019-04-15.
//

#include <string.h>

template<typename T>
comm::Region<T>::Region() {}

template<typename T>
comm::Region<T>::Region(const T x_start, const T y_start, const T z_start,
                        const T x_end, const T y_end, const T z_end)
        :low(&x_low), high(&x_high) {
    x_low = x_start;
    y_low = y_start;
    z_low = z_start;
    x_high = x_end;
    y_high = y_end;
    z_high = z_end;
}

template<typename T>
comm::Region<T>::Region(const comm::Region<T> &r):low(&x_low), high(&x_high) {
    memcpy(this->low, r.low, DIMENSION_SIZE * sizeof(T));
    memcpy(this->high, r.high, DIMENSION_SIZE * sizeof(T));
}

template<typename T>
comm::Region<T> &comm::Region<T>::operator=(const comm::Region<T> &r) {
    memcpy(this->low, r.low, DIMENSION_SIZE * sizeof(T));
    memcpy(this->high, r.high, DIMENSION_SIZE * sizeof(T));
    return *this;
}
