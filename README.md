# libcomm

> communication and domain lib for particle simulation.

## Build
```bash
$ cmake -B./build -H./ \
 -DOpenMP_ENABLE_FLAG=ON \
 -DMPI_ENABLE_FLAG=ON \
 -DTEST_ENABLE_FLAG=ON \
 -DTEST_MPI_ENABLE_FLAG=ON
$ cmake --build ./build -j 4
```