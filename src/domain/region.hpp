//
// Created by genshen back to 2018-12-30.
//

#ifndef COMM_DOMAIN_REGION_H
#define COMM_DOMAIN_REGION_H

template<typename T>
struct Region {
    T low[DIMENSION];
    T high[DIMENSION];

    T &x_low = low[0];
    T &y_low = low[1];
    T &z_low = low[2];
    T &x_high = high[0];
    T &y_high = high[1];
    T &z_high = high[2];
};

#endif // COMM_DOMAIN_REGION_H