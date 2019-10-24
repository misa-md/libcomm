set(HEADER_FILES
        comm/domain/domain.h
        comm/domain/bcc_domain.h
        comm/domain/colored_domain.h
        comm/domain/builder.h
        comm/domain/region.hpp
        comm/preset/comm_forwarding_region.h
        comm/preset/sector_forwarding_direction.h
        comm/preset/sector_forwarding_region.h
        comm/comm.hpp
        comm/packer.h
        comm/region_packer.h
        comm/types_define.h
        )

set(SOURCE_FILES
        comm/domain/domain.cc
        comm/domain/bcc_domain.cc
        comm/domain/colored_domain.cc
        comm/domain/builder.inl
        comm/domain/region.inl
        comm/preset/comm_forwarding_region.cc
        comm/preset/sector_forwarding_region.cc
        comm/comm.inl
        )

# set(MAIN_FILE main.cpp)

# set binary path and lib storage path for compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# make lib
add_library(${COMM_LIB_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(
        ${COMM_LIB_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

##########################
##### install files
##########################
set(TARGETS_EXPORT_NAME ${COMM_LIB_NAME}-targets)
set(COMM_LIB_CMAKE_DIR lib/cmake/${COMM_LIB_NAME})
set(COMM_LIB_CMAKE_NAMESPACE ${COMM_LIB_NAME}::)
set(COMM_LIB_CMAKE_CONFIG_IN "${COMM_LIB_NAME}-config.cmake.in")
set(COMM_LIB_CMAKE_CONFIG_OUT "${CMAKE_CURRENT_BINARY_DIR}/${COMM_LIB_NAME}-config.cmake")

# install bin and libs.
install(TARGETS ${COMM_EXECUTE_BIN_NAME} ${COMM_LIB_NAME}
        EXPORT ${TARGETS_EXPORT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

# install header.
install(DIRECTORY "./"
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
        PATTERN "*.inl"
        )

# install export target and config for find_package
install(EXPORT ${TARGETS_EXPORT_NAME} DESTINATION ${COMM_LIB_CMAKE_DIR} NAMESPACE ${COMM_LIB_CMAKE_NAMESPACE}) # todo
include(CMakePackageConfigHelpers)
configure_package_config_file(
        ${COMM_LIB_CMAKE_CONFIG_IN} ${COMM_LIB_CMAKE_CONFIG_OUT}
        INSTALL_DESTINATION ${COMM_LIB_CMAKE_DIR}
)
install(FILES ${COMM_LIB_CMAKE_CONFIG_OUT} DESTINATION ${COMM_LIB_CMAKE_DIR})
