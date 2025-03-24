cmake_minimum_required(VERSION 3.21)
project(hyhound-test CXX)
include(GoogleTest)

# Dependencies
find_package(GTest REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(guanaqo REQUIRED COMPONENTS BLAS)
if (NOT TARGET hyhound::hyhound)
    find_package(hyhound REQUIRED OPTIONAL_COMPONENTS OCP)
endif()

add_executable(tests
    "test-downdate.cpp"
    "test-updowndate.cpp"
)
target_link_libraries(tests PRIVATE
    guanaqo::blas hyhound::hyhound hyhound::warnings
    GTest::gtest_main GTest::gmock
    Eigen3::Eigen)

if (TARGET hyhound::ocp)
    target_sources(tests PRIVATE "test-ocp.cpp")
    target_link_libraries(tests PRIVATE hyhound::ocp)
endif()

if (NOT CMAKE_CROSSCOMPILING OR HYHOUND_FORCE_TEST_DISCOVERY)
    gtest_discover_tests(tests DISCOVERY_TIMEOUT 60)
endif()
add_executable(hyhound::tests ALIAS tests)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.21) # Copy DLLs for test executable
    add_custom_command(TARGET tests PRE_LINK COMMAND
        COMMAND ${CMAKE_COMMAND} -E
            $<IF:$<BOOL:$<TARGET_RUNTIME_DLLS:tests>>,copy,true>
            $<TARGET_RUNTIME_DLLS:tests> $<TARGET_FILE_DIR:tests>
        COMMAND_EXPAND_LISTS)
endif()
