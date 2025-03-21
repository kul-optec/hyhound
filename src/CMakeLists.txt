include(cmake/Library.cmake)
include(cmake/BLAS.cmake)
find_package(guanaqo REQUIRED)
if (HYHOUND_WITH_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Configuration options
# ------------------------------------------------------------------------------
add_library(config INTERFACE)
configure_file("config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/config/include/hyhound/config.hpp" @ONLY)
target_sources(config INTERFACE FILE_SET headers TYPE HEADERS
    BASE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/config/include"
    FILES "${CMAKE_CURRENT_BINARY_DIR}/config/include/hyhound/config.hpp"
)
target_compile_definitions(config INTERFACE
    $<$<BOOL:${HYHOUND_VERIFY_ASSUMPTIONS}>:HYHOUND_VERIFY_ASSUMPTIONS>
    $<$<BOOL:${HYHOUND_WITH_MKL}>:HYHOUND_WITH_MKL>
    $<$<BOOL:${HYHOUND_WITH_OPENMP}>:HYHOUND_WITH_OPENMP>
)
hyhound_configure_interface_library(config EXPORT_PREFIX hyhound)

# Utility functions
# ------------------------------------------------------------------------------
add_library(util INTERFACE)
target_sources(util INTERFACE FILE_SET headers TYPE HEADERS
    BASE_DIRS "util/include"
    FILES "util/include/hyhound/assume.hpp"
          "util/include/hyhound/cneg.hpp"
          "util/include/hyhound/loop.hpp"
          "util/include/hyhound/lut.hpp"
          "util/include/hyhound/matrix-view.hpp"
          "util/include/hyhound/openmp.h"
          "util/include/hyhound/unroll.h"
)
target_link_libraries(util INTERFACE hyhound::config guanaqo::guanaqo)
hyhound_configure_interface_library(util)
if (HYHOUND_WITH_OPENMP)
    target_link_libraries(util INTERFACE OpenMP::OpenMP_CXX)
endif()

# Linear algebra routines and BLAS interface
# ------------------------------------------------------------------------------
add_library(linalg "linalg/src/blas-interface.cpp")
target_sources(linalg PUBLIC FILE_SET headers TYPE HEADERS
    BASE_DIRS "linalg/include"
    FILES "linalg/include/hyhound/linalg/blas-interface.hpp"
          "linalg/include/hyhound/linalg/blas.hpp"
          "linalg/include/hyhound/linalg/lapack.hpp"
)
target_link_libraries(linalg PUBLIC hyhound::util hyhound::blas-lapack-lib)
target_compile_features(linalg PRIVATE cxx_std_20)
hyhound_configure_library(linalg EXPORT_PREFIX hyhound)

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
hyhound_configure_library(hyhound)

include("hyhound/src/micro-kernels/Codegen.cmake")
codegen_hyhound_microkernels(hyhound)


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
