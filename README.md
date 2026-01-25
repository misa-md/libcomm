# libcomm

> communication and domain lib for particle simulation.

## Install dependencies
Install [pkg](https://github.com/genshen/pkg) tool from [release](https://github.com/genshen/pkg/release) page, and add its path to `PATH` environment, then run:
```bash
$ pkg fetch
$ pkg install
```

## Build lib
```bash
$ cmake -B./build -H./ \
 -DCOMM_OpenMP_ENABLE_FLAG=ON \
 -DCOMM_MPI_ENABLE_FLAG=ON \
 -DCOMM_TEST_ENABLE_FLAG=ON \
 -DCOMM_TEST_MPI_ENABLE_FLAG=ON
$ cmake --build ./build -j 4
```
