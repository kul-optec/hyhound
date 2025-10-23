cmake_minimum_required(VERSION 3.24...4.1)
set(CMAKE_CXX_SCAN_FOR_MODULES Off)
project(hyhound
    VERSION 1.0.2
    DESCRIPTION "Hyperbolic Householder transformations for Cholesky factorization up- and downdates"
    HOMEPAGE_URL "https://github.com/kul-optec/hyhound"
    LANGUAGES CXX
)
include(CTest)
set(PY_VERSION_SUFFIX "a2")

# Make sure that the Python and CMake versions match
set(PY_FULL_VERSION ${PROJECT_VERSION}${PY_VERSION_SUFFIX})
if (DEFINED PY_BUILD_CMAKE_PROJECT_VERSION)
    if (NOT "${PY_BUILD_CMAKE_PROJECT_VERSION}" MATCHES "^${PY_FULL_VERSION}$")
        message(FATAL_ERROR "Version number does not match "
                             "(${PY_BUILD_CMAKE_PROJECT_VERSION} - ${PY_FULL_VERSION}).")
    endif()
endif()

# Options
include(CMakeDependentOption)
set(HYHOUND_DENSE_REAL_TYPE "double" "float" CACHE STRING
    "The floating point types that the functions are instantiated for")
set(HYHOUND_DENSE_INDEX_TYPE "long long" CACHE STRING
    "The main integer type for indices and sizes")
# Target options
option(HYHOUND_HAVE_TWO_512_FMA_UNITS
    "Does the CPU have 0.5 CPI 512-bit FMA" Off)
set(HYHOUND_MAX_HYH_KERNEL_WIDTH "8" CACHE STRING
    "Kernel size limit for the 'diagonal' kernels")
set(HYHOUND_MAX_HYH_KERNEL_HEIGHT "32" CACHE STRING
    "Kernel size limit for the 'below-diagonal' kernels")
# Enable/disable optional components
option(HYHOUND_WITH_TESTS
    "Build the tests" ${BUILD_TESTING})
option(HYHOUND_WITH_OCP
    "Build the OCP solvers" Off)
option(HYHOUND_WITH_BENCHMARKS
    "Build the benchmarks" Off)
option(HYHOUND_WITH_PYTHON
    "Build the Python interface" Off)
option(HYHOUND_WITH_PYTHON_DISPATCH
    "Build the AVX2/AVX-512 dispatch logic for the Python interface" Off)
# Developer options
option(HYHOUND_WARNINGS_AS_ERRORS
    "Enable -Werror or /WX" Off)
option(HYHOUND_WITH_ACCURATE_BUILD_TIME
    "Update the build time on every build" On)
option(HYHOUND_VERIFY_ASSUMPTIONS
    "Check assumptions at run time instead of blindly assuming" Off)
option(HYHOUND_FORCE_TEST_DISCOVERY
    "Query the test executables even when cross-compiling" Off)
option(HYHOUND_DEVELOPER_MODE
    "Enable developer options such as ccache and colored compiler output" Off)

# Installation paths
include(GNUInstallDirs)
set(HYHOUND_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}"
    CACHE PATH "Installation directory for archives and libraries")
set(HYHOUND_INSTALL_CMAKEDIR "${HYHOUND_INSTALL_LIBDIR}/cmake/hyhound"
    CACHE PATH "Installation directory for CMake configuration files")
set(HYHOUND_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}"
    CACHE PATH "Installation directory for binaries and DLLs")
set(HYHOUND_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}"
    CACHE PATH "Installation directory for headers")
option(HYHOUND_STANDALONE
    "Install with relative RPATH to locate its own shared libraries" On)

# Development
if (PROJECT_IS_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS On)
endif()
if (HYHOUND_DEVELOPER_MODE)
    include(cmake/Develop.cmake)
endif()

# Compiler warnings
include(cmake/Warnings.cmake)
add_warnings_target(warnings ${HYHOUND_WARNINGS_AS_ERRORS})
add_library(hyhound::warnings ALIAS warnings)

# Other compiler and CMake options
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "_d")
set(CMAKE_RELWITHDEBINFO_POSTFIX "_rd")
set(CMAKE_MINSIZEREL_POSTFIX "_rs")

# Locating dependencies
if (DEFINED ENV{MKLROOT})
    list(PREPEND CMAKE_FIND_ROOT_PATH "$ENV{MKLROOT}")
endif()

# Libraries
add_subdirectory(src)

# Tests
if (HYHOUND_WITH_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# Benchmarks
if (HYHOUND_WITH_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Python interface
if (HYHOUND_WITH_PYTHON)
    add_subdirectory(python)
endif()
