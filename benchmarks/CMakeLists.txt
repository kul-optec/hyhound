find_package(benchmark REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)
find_package(guanaqo REQUIRED COMPONENTS BLAS)
enable_testing()
include(GoogleTest)

add_executable(benchmark-hyh hyh.cpp)
target_link_libraries(benchmark-hyh PRIVATE
    Eigen3::Eigen benchmark::benchmark guanaqo::blas
    hyhound::hyhound hyhound::warnings)

add_library(ocp
    ocp/riccati/factor.cpp
    ocp/riccati/solve.cpp
    ocp/riccati/update.cpp
    ocp/schur/factor.cpp
    ocp/schur/solve.cpp
    ocp/schur/update.cpp
)
target_link_libraries(ocp
    PUBLIC Eigen3::Eigen hyhound::config guanaqo::blas
    PRIVATE hyhound::hyhound hyhound::warnings)
target_include_directories(ocp
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(benchmark-ocp ocp.cpp)
target_link_libraries(benchmark-ocp PRIVATE
    ocp benchmark::benchmark hyhound::hyhound hyhound::warnings)

add_executable(benchmark-ocp-test ocp-test.cpp)
target_link_libraries(benchmark-ocp-test PRIVATE
    ocp benchmark::benchmark GTest::gtest_main
    hyhound::hyhound hyhound::warnings)

if (NOT CMAKE_CROSSCOMPILING OR HYHOUND_FORCE_TEST_DISCOVERY)
    gtest_discover_tests(benchmark-ocp-test DISCOVERY_TIMEOUT 60)
endif()
