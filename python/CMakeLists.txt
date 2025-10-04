cmake_minimum_required(VERSION 4.1)

# Find Python and nanobind
if (CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    find_package(Python REQUIRED COMPONENTS Development.Module
                        OPTIONAL_COMPONENTS Development.SABIModule)
else()
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module
                        OPTIONAL_COMPONENTS Development.SABIModule)
endif()
find_package(nanobind REQUIRED CONFIG)

# Determine whether we're doing a free-threaded build
set(FREE_THREADED)
if (PY_BUILD_CMAKE_PACKAGE_ABI_TAG MATCHES "(^|;)cp([0-9dmu]+)t($|;)")
    message(STATUS "Enabling nanobind FREE_THREADED mode")
    set(FREE_THREADED "FREE_THREADED")
endif()

nanobind_add_module(_hyhound STABLE_ABI ${FREE_THREADED} "hyhound.py.cpp")
target_compile_definitions(_hyhound PRIVATE
    MODULE_NAME=$<TARGET_FILE_BASE_NAME:_hyhound>
    VERSION_INFO="${PY_FULL_VERSION}"
)
target_link_libraries(_hyhound PRIVATE hyhound::hyhound)
set_target_properties(_hyhound PROPERTIES
    RELEASE_POSTFIX ""
    DEBUG_POSTFIX ""
    RELWITHDEBINFO_POSTFIX ""
    MINSIZEREL_POSTFIX ""
)

# Hide all symbols by default (including external libraries on Linux)
set_target_properties(_hyhound PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN true)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_options(_hyhound PRIVATE "LINKER:--exclude-libs,ALL")
endif()

# Install the Python module
install(TARGETS _hyhound
        EXCLUDE_FROM_ALL
        COMPONENT python_modules
        DESTINATION ${PY_BUILD_CMAKE_IMPORT_NAME})

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
        MODULE ${PY_BUILD_CMAKE_IMPORT_NAME}.$<TARGET_FILE_BASE_NAME:_hyhound>
        OUTPUT ${PY_BUILD_CMAKE_IMPORT_NAME}/$<TARGET_FILE_BASE_NAME:_hyhound>.pyi
        PYTHON_PATH "."
    )
endif()
