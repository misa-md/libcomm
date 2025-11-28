if (NOT (COMM_TEST_BUILD_ENABLE_FLAG))
    set(COMM_TEST_MPI_ENABLE_FLAG OFF)  # if we disabled KMC_MPI_ENABLE_FLAG, we cann't also use KMC_TEST_MPI_ENABLE_FLAG.
endif ()
MESSAGE(STATUS "MPI in test: ${COMM_TEST_MPI_ENABLE_FLAG}")

add_subdirectory(unit)
