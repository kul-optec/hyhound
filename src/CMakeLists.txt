include(cmake/Library.cmake)
find_package(guanaqo REQUIRED)

# Configuration options
# ------------------------------------------------------------------------------
add_library(config INTERFACE)
list(GET HYHOUND_DENSE_REAL_TYPE 0 HYHOUND_DENSE_REAL_TYPE_0)
set(HYHOUND_WITH_DOUBLE Off)
set(HYHOUND_WITH_FLOAT Off)
if ("double" IN_LIST HYHOUND_DENSE_REAL_TYPE)
    set(HYHOUND_WITH_DOUBLE On)
endif()
if ("float" IN_LIST HYHOUND_DENSE_REAL_TYPE)
    set(HYHOUND_WITH_FLOAT On)
endif()
configure_file("config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/config/include/hyhound/config.hpp" @ONLY)
target_sources(config INTERFACE FILE_SET headers TYPE HEADERS
    BASE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/config/include"
    FILES "${CMAKE_CURRENT_BINARY_DIR}/config/include/hyhound/config.hpp"
)
target_compile_definitions(config INTERFACE
    $<$<BOOL:${HYHOUND_VERIFY_ASSUMPTIONS}>:HYHOUND_VERIFY_ASSUMPTIONS>
    $<$<BOOL:${HYHOUND_WITH_OCP}>:HYHOUND_WITH_OCP>
)
hyhound_configure_interface_library(config)

# Utility functions
# ------------------------------------------------------------------------------
add_library(util INTERFACE)
target_sources(util INTERFACE FILE_SET headers TYPE HEADERS
    BASE_DIRS "util/include"
    FILES "util/include/hyhound/assume.hpp"
          "util/include/hyhound/cneg.hpp"
          "util/include/hyhound/loop.hpp"
          "util/include/hyhound/lut.hpp"
          "util/include/hyhound/unroll.h"
)
target_link_libraries(util INTERFACE hyhound::config guanaqo::guanaqo)
hyhound_configure_interface_library(util)

# Main library for Cholesky factorization up/downdates
# ------------------------------------------------------------------------------
add_library(hyhound)
target_sources(hyhound PUBLIC FILE_SET headers TYPE HEADERS
    BASE_DIRS "hyhound/include"
    FILES "hyhound/include/hyhound/updown.hpp"
          "hyhound/include/hyhound/householder-updowndate.hpp"
          "hyhound/include/hyhound/householder-updowndate-serial.tpp"
          "hyhound/include/hyhound/micro-kernels/common.hpp"
          "hyhound/include/hyhound/micro-kernels/householder-updowndate-diag.tpp"
          "hyhound/include/hyhound/micro-kernels/householder-updowndate-full.tpp"
          "hyhound/include/hyhound/micro-kernels/householder-updowndate-tail.tpp"
          "hyhound/include/hyhound/micro-kernels/householder-updowndate.hpp"
          "hyhound/include/hyhound/micro-kernels/matrix-accessor.hpp"
          "hyhound/include/hyhound/householder-updowndate-micro-kernels.tpp"
)
target_link_libraries(hyhound PUBLIC hyhound::config hyhound::util)
target_link_libraries(hyhound PRIVATE hyhound::warnings)
hyhound_configure_library(hyhound)

include("hyhound/src/micro-kernels/Codegen.cmake")
codegen_hyhound_microkernels(hyhound)

# Library for solving OCPs using Riccati recursion or Schur complement method
# ------------------------------------------------------------------------------
if (HYHOUND_WITH_OCP)
    find_package(guanaqo REQUIRED COMPONENTS BLAS)
    find_package(Eigen3 REQUIRED)
    add_library(ocp
        "ocp/src/riccati/factor.cpp"
        "ocp/src/riccati/solve.cpp"
        "ocp/src/riccati/update.cpp"
        "ocp/src/schur/factor.cpp"
        "ocp/src/schur/solve.cpp"
        "ocp/src/schur/update.cpp"
    )
    target_sources(ocp PUBLIC FILE_SET headers TYPE HEADERS
        BASE_DIRS "ocp/include"
        FILES "ocp/include/hyhound/ocp/riccati.hpp"
              "ocp/include/hyhound/ocp/schur.hpp"
    )
    target_link_libraries(ocp PUBLIC hyhound::hyhound PRIVATE hyhound::warnings)
    target_link_libraries(ocp PUBLIC Eigen3::Eigen guanaqo::blas)
    hyhound_configure_library(ocp EXPORT_PREFIX hyhound)
endif()


# ==============================================================================

# Save the project version to a generated file
configure_file("cmake/hyhound-version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/hyhound-version.h" @ONLY)
target_sources(hyhound PUBLIC FILE_SET headers TYPE HEADERS
    BASE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/include"
    FILES "${CMAKE_CURRENT_BINARY_DIR}/include/hyhound-version.h"
)

# Build time
if (HYHOUND_WITH_ACCURATE_BUILD_TIME)
    set(HYHOUND_BUILD_TIME_CPP "${CMAKE_CURRENT_BINARY_DIR}/hyhound-build-time.cpp")
    add_custom_target(hyhound-build-time-generate
        BYPRODUCTS ${HYHOUND_BUILD_TIME_CPP}
        COMMAND ${CMAKE_COMMAND}
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildTime.cmake"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    target_sources(hyhound PRIVATE ${HYHOUND_BUILD_TIME_CPP})
    add_dependencies(hyhound hyhound-build-time-generate)
else()
    set(HYHOUND_BUILD_TIME_CPP "${CMAKE_CURRENT_BINARY_DIR}/hyhound-build-time.cpp")
    if (NOT EXISTS ${HYHOUND_BUILD_TIME_CPP})
        include(cmake/BuildTime.cmake)
    endif()
    target_sources(hyhound PRIVATE ${HYHOUND_BUILD_TIME_CPP})
endif()

# Installation
include(cmake/Install.cmake)
