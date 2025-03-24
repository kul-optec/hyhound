include(${PROJECT_SOURCE_DIR}/cmake/Debug.cmake)

# Set the runtime linker/loader search paths to make hyhound stand-alone
cmake_path(RELATIVE_PATH HYHOUND_INSTALL_LIBDIR
           BASE_DIRECTORY HYHOUND_INSTALL_BINDIR
           OUTPUT_VARIABLE HYHOUND_INSTALL_LIBRELBINDIR)

function(hyhound_add_if_target_exists OUT)
    foreach(TGT IN LISTS ARGN)
        if (TARGET ${TGT})
            list(APPEND ${OUT} ${TGT})
        endif()
    endforeach()
    set(${OUT} ${${OUT}} PARENT_SCOPE)
endfunction()

include(CMakePackageConfigHelpers)

set(HYHOUND_INSTALLED_COMPONENTS)
macro(hyhound_install_config PKG COMP)
    # Install the target CMake definitions
    install(EXPORT hyhound${PKG}Targets
        FILE hyhound${PKG}Targets.cmake
        DESTINATION "${HYHOUND_INSTALL_CMAKEDIR}"
            COMPONENT ${COMP}
        NAMESPACE hyhound::)
    # Add all targets to the build tree export set
    export(EXPORT hyhound${PKG}Targets
        FILE "${PROJECT_BINARY_DIR}/hyhound${PKG}Targets.cmake"
        NAMESPACE hyhound::)
    # Generate the config file that includes the exports
    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PKG}Config.cmake.in"
        "${PROJECT_BINARY_DIR}/hyhound${PKG}Config.cmake"
        INSTALL_DESTINATION "${HYHOUND_INSTALL_CMAKEDIR}"
        NO_SET_AND_CHECK_MACRO)
    write_basic_package_version_file(
        "${PROJECT_BINARY_DIR}/hyhound${PKG}ConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}"
        COMPATIBILITY SameMinorVersion)
    # Install the hyhoundConfig.cmake and hyhoundConfigVersion.cmake
    install(FILES
        "${PROJECT_BINARY_DIR}/hyhound${PKG}Config.cmake"
        "${PROJECT_BINARY_DIR}/hyhound${PKG}ConfigVersion.cmake"
        DESTINATION "${HYHOUND_INSTALL_CMAKEDIR}"
            COMPONENT ${COMP})
    list(APPEND HYHOUND_OPTIONAL_COMPONENTS ${PKG})
endmacro()

macro(hyhound_install_cmake FILES COMP)
    # Install a CMake script
    install(FILES ${FILES}
        DESTINATION "${HYHOUND_INSTALL_CMAKEDIR}"
            COMPONENT ${COMP})
endmacro()

set(HYHOUND_INSTALLED_TARGETS_MSG "\nSummary of hyhound components and targets to install:\n\n")

# Install the hyhound core libraries
set(HYHOUND_CORE_HIDDEN_TARGETS warnings common_options)
set(HYHOUND_CORE_TARGETS config util hyhound)
if (HYHOUND_CORE_TARGETS)
    install(TARGETS ${HYHOUND_CORE_HIDDEN_TARGETS} ${HYHOUND_CORE_TARGETS}
        EXPORT hyhoundCoreTargets
        RUNTIME DESTINATION "${HYHOUND_INSTALL_BINDIR}"
            COMPONENT lib
        LIBRARY DESTINATION "${HYHOUND_INSTALL_LIBDIR}"
            COMPONENT lib
            NAMELINK_COMPONENT dev
        ARCHIVE DESTINATION "${HYHOUND_INSTALL_LIBDIR}"
            COMPONENT dev
        FILE_SET headers DESTINATION "${HYHOUND_INSTALL_INCLUDEDIR}"
            COMPONENT dev)
    hyhound_install_config(Core dev)
    list(JOIN HYHOUND_CORE_TARGETS ", " TGTS)
    string(APPEND HYHOUND_INSTALLED_TARGETS_MSG " * Core:  ${TGTS}\n")
    list(APPEND HYHOUND_INSTALL_TARGETS ${HYHOUND_CORE_TARGETS})
endif()

# Install the hyhound OCP libraries
set(HYHOUND_OCP_TARGETS)
hyhound_add_if_target_exists(HYHOUND_OCP_TARGETS ocp)
if (HYHOUND_OCP_TARGETS)
    install(TARGETS ${HYHOUND_OCP_TARGETS}
        EXPORT hyhoundOCPTargets
        RUNTIME DESTINATION "${HYHOUND_INSTALL_BINDIR}"
            COMPONENT lib
        LIBRARY DESTINATION "${HYHOUND_INSTALL_LIBDIR}"
            COMPONENT lib
            NAMELINK_COMPONENT dev
        ARCHIVE DESTINATION "${HYHOUND_INSTALL_LIBDIR}"
            COMPONENT dev
        FILE_SET headers DESTINATION "${HYHOUND_INSTALL_INCLUDEDIR}"
            COMPONENT dev)
    hyhound_install_config(OCP dev)
    list(JOIN HYHOUND_OCP_TARGETS ", " TGTS)
    string(APPEND HYHOUND_INSTALLED_TARGETS_MSG " * OCP:   ${TGTS}\n")
    list(APPEND HYHOUND_INSTALL_TARGETS ${HYHOUND_OCP_TARGETS})
endif()

# Install the debug files
foreach(target IN LISTS HYHOUND_CORE_TARGETS HYHOUND_OCP_TARGETS)
    get_target_property(target_type ${target} TYPE)
    if (${target_type} STREQUAL "SHARED_LIBRARY")
        hyhound_install_debug_syms(${target} debug
                                  ${HYHOUND_INSTALL_LIBDIR}
                                  ${HYHOUND_INSTALL_BINDIR})
    elseif (${target_type} STREQUAL "EXECUTABLE")
        hyhound_install_debug_syms(${target} debug
                                  ${HYHOUND_INSTALL_BINDIR}
                                  ${HYHOUND_INSTALL_BINDIR})
    endif()
endforeach()

# Make stand-alone
if (HYHOUND_STANDALONE)
    foreach(target IN LISTS HYHOUND_CORE_TARGETS HYHOUND_OCP_TARGETS)
        set_target_properties(${TGT} PROPERTIES
            INSTALL_RPATH "$ORIGIN;$ORIGIN/${HYHOUND_INSTALL_LIBRELBINDIR}")
    endforeach()
endif()

# Generate the main config file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in"
    "${PROJECT_BINARY_DIR}/hyhoundConfig.cmake"
    INSTALL_DESTINATION "${HYHOUND_INSTALL_CMAKEDIR}"
    NO_SET_AND_CHECK_MACRO)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/hyhoundConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMinorVersion)
# Install the main hyhoundConfig.cmake and hyhoundConfigVersion.cmake files
install(FILES
    "${PROJECT_BINARY_DIR}/hyhoundConfig.cmake"
    "${PROJECT_BINARY_DIR}/hyhoundConfigVersion.cmake"
    DESTINATION "${HYHOUND_INSTALL_CMAKEDIR}"
        COMPONENT dev)

# Print the components and targets we're going to install
message(${HYHOUND_INSTALLED_TARGETS_MSG})
