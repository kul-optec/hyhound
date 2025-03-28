set(CMAKE_EXPORT_COMPILE_COMMANDS On)
set(CMAKE_C_COMPILER_LAUNCHER "sccache")
set(CMAKE_CXX_COMPILER_LAUNCHER "sccache")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fdiagnostics-color=always)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang$")
    add_compile_options(-fcolor-diagnostics)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-ffunction-sections)
    add_link_options(LINKER:--gc-sections)
endif()
