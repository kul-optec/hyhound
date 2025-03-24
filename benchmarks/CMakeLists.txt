find_package(benchmark REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)
find_package(guanaqo REQUIRED COMPONENTS BLAS)

add_executable(benchmark-hyh hyh.cpp)
target_link_libraries(benchmark-hyh PRIVATE
    Eigen3::Eigen benchmark::benchmark guanaqo::blas
    hyhound::hyhound hyhound::warnings)

if (HYHOUND_WITH_OCP)
    add_executable(benchmark-ocp ocp.cpp)
    target_link_libraries(benchmark-ocp PRIVATE
        benchmark::benchmark
        hyhound::ocp hyhound::warnings)
endif()
