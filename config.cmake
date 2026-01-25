set(CMAKE_CXX_STANDARD 11)
set(COMM_VERSION "0.4.0")
set(ProjectName comm)

#############
## options ##
#############
# change to mpicc and mpicxx
#set(CMAKE_C_COMPILER mpicc -cc=gcc -cxx=g++)
#set(CMAKE_CXX_COMPILER mpicxx -cc=gcc -cxx=g++)

option(COMM_OpenMP_ENABLE_FLAG "Use OpenMP" OFF) #change this flag to OFF to disable OpenMP
option(COMM_MPI_ENABLE_FLAG "Use MPI library" ON) #change this flag to false to disable mpi
option(COMM_TEST_BUILD_ENABLE_FLAG "Enable building test" ON) # enable test
option(COMM_TEST_MPI_ENABLE_FLAG "Enable MPI in test" ON) # enable mpi in test, its value depends on option COMM_MPI_ENABLE_FLAG.

## architecture ralated values.
# option(ARCH_SW "Enable sunway athread" OFF) # enable sunway athread if its running on sunway system.

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

if (CMAKE_BUILD_TYPE MATCHES "^(Debug|DEBUG|debug)$")
    set(COMM_DEBUG_ENABLE_FLAG ON)
endif ()

#############
## const ##
#############
set(COMM_LIB_NAME ${ProjectName}) # todo use PARENT_SCOPE to modify globle variable.

# test
set(COMM_UINT_TEST_NAME "unit-test")