cmake_minimum_required(VERSION 4.1)

set(HYHOUND_PYTHON_POSTFIX "_generic" CACHE STRING
    "Additional postfix in Python module file name")

# Find Python and nanobind
if (CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    find_package(Python REQUIRED COMPONENTS Development.Module
                        OPTIONAL_COMPONENTS Development.SABIModule)
else()
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module
                        OPTIONAL_COMPONENTS Development.SABIModule)
endif()
find_package(nanobind REQUIRED CONFIG)

# Determine whether we're doing a stable or a free-threaded build
set(USE_SABI)
set(NANOBIND_TARGET "nanobind-hyhound")
if (PY_BUILD_CMAKE_PACKAGE_ABI_TAG MATCHES "(^|;)cp([0-9dmu]+)t($|;)")
    message(STATUS "Free-threaded build")
    string(APPEND NANOBIND_TARGET "-ft")
elseif (DEFINED PY_BUILD_CMAKE_PACKAGE_USE_SABI)
    message(STATUS "Stable ABI build (${PY_BUILD_CMAKE_PACKAGE_USE_SABI})")
    string(APPEND NANOBIND_TARGET "-abi3")
    set(USE_SABI USE_SABI ${PY_BUILD_CMAKE_PACKAGE_USE_SABI})
    set(USE_LIMITED_API ${PY_BUILD_CMAKE_PACKAGE_LIMITED_API})
endif()

# Build the nanobind library
nanobind_build_library(${NANOBIND_TARGET})
if (USE_LIMITED_API)
    target_compile_definitions(${NANOBIND_TARGET} PRIVATE
        Py_LIMITED_API=${USE_LIMITED_API})
endif()
install(TARGETS ${NANOBIND_TARGET}
        EXCLUDE_FROM_ALL
        COMPONENT python_nanobind
        DESTINATION ${PY_BUILD_CMAKE_IMPORT_NAME})

# Configure the dynamic loading paths to locate the nanobind libraries.
function(hyhound_rpath_origin tgt)
    if(APPLE)
        set_target_properties(${tgt} PROPERTIES INSTALL_RPATH "@loader_path")
    else()
        set_target_properties(${tgt} PROPERTIES INSTALL_RPATH "$ORIGIN")
    endif()
endfunction()

# Configure common options for Python extension modules
function(hyhound_configure_python_ext tgt)
    target_link_libraries(${tgt} PRIVATE ${NANOBIND_TARGET})
    hyhound_rpath_origin(${tgt})
    nanobind_compile_options(${tgt})
    nanobind_link_options(${tgt})
    hyhound_configure_visibility(${tgt})
    target_compile_definitions(${tgt} PRIVATE
        MODULE_NAME=$<TARGET_FILE_BASE_NAME:${tgt}>
        VERSION_INFO="${PY_FULL_VERSION}")
endfunction()

# Build the Python extension
Python_add_library(_hyhound MODULE WITH_SOABI ${USE_SABI} "hyhound.py.cpp")
hyhound_configure_python_ext(_hyhound)
target_link_libraries(_hyhound PRIVATE hyhound::hyhound)
set_target_properties(_hyhound PROPERTIES
    RELEASE_POSTFIX "${HYHOUND_PYTHON_POSTFIX}"
    DEBUG_POSTFIX "${HYHOUND_PYTHON_POSTFIX}"
    RELWITHDEBINFO_POSTFIX "${HYHOUND_PYTHON_POSTFIX}"
    MINSIZEREL_POSTFIX "${HYHOUND_PYTHON_POSTFIX}")

# Install the Python module
install(TARGETS _hyhound
        EXCLUDE_FROM_ALL
        COMPONENT python_modules
        DESTINATION ${PY_BUILD_CMAKE_IMPORT_NAME})

# Build the dispatch library
if (HYHOUND_WITH_PYTHON_DISPATCH)
    Python_add_library(_dispatch MODULE WITH_SOABI ${USE_SABI} "dispatch.py.cpp")
    hyhound_configure_python_ext(_dispatch)
    find_package(CpuFeatures REQUIRED)
    target_link_libraries(_dispatch PRIVATE CpuFeatures::cpu_features)
    hyhound_configure_visibility(_dispatch)
    install(TARGETS _dispatch
        EXCLUDE_FROM_ALL
        COMPONENT python_modules
        DESTINATION ${PY_BUILD_CMAKE_IMPORT_NAME})
endif()

# Generate stubs for the Python module
set(WITH_PY_STUBS_DEFAULT On)
if (CMAKE_CROSSCOMPILING)
    set(WITH_PY_STUBS_DEFAULT Off)
endif()
option(WITH_PY_STUBS
    "Generate Python stub files (.pyi) for the Python module."
    ${WITH_PY_STUBS_DEFAULT})
if (WITH_PY_STUBS)
    nanobind_add_stub(
        $<TARGET_FILE_BASE_NAME:_hyhound>_stub
        INSTALL_TIME
        COMPONENT python_stubs
        MODULE $<TARGET_FILE_BASE_NAME:_hyhound>
        OUTPUT "${PY_BUILD_CMAKE_IMPORT_NAME}/_hyhound.pyi"
        PYTHON_PATH "${PY_BUILD_CMAKE_IMPORT_NAME}"
        EXCLUDE_FROM_ALL
        VERBOSE)
endif()
