# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/test_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_config.h"
)

set(TEST_HEADER_FILES
        test_config.h
        domain/domain_test_utils.h
        )

set(TEST_SOURCE_FILES
        main.cpp
        domain/domain_test.cpp
        domain/domain_test_utils.cpp
        domain/region_test.cpp
        domain/bcc_domain_test.cpp
        domain/colored_domain_test.cpp
        preset/comm_forwarding_region_test.cpp
        preset/sector_forwarding_direction_test.cpp
        preset/sector_forwarding_region_test.cpp
        region_packer_test.cpp
        )

# set binary path and lib storage path wile compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_executable(${COMM_UINT_TEST_NAME} ${TEST_SOURCE_FILES} ${TEST_HEADER_FILES})

# Extra linking for the project.
# Extra linking for the project.
target_link_libraries(${COMM_UINT_TEST_NAME}
        PRIVATE GTest::gtest GTest::gtest_main ${COMM_LIB_NAME} ${EXTRA_LIBS}
        )
# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME unit-test COMMAND ${COMM_UINT_TEST_NAME})


##########################
##### install files
##########################
install(TARGETS ${COMM_UINT_TEST_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
