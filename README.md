# hyhound

**Hy**perbolic **Ho**useholder transformations for **U**p- a**n**d **D**owndating Cholesky factorizations.

## Instructions for reproducing the benchmark results (Linux)

Requirements: [CMake](https://cmake.org/), [Ninja](https://ninja-build.org/),
[Conan](https://conan.io/), [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html).

```sh
git clone https://github.com/kul-optec/hyhound
cd hyhound
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes"
conan profile detect ||:
conan install . --build=missing -pr profiles/desktop \
    -s build_type=Release \
    -c tools.build:skip_test=True \
    -o \&:with_openmp=True -o \&:with_openblas=False -o \&:with_mkl=True \
    -o \&:with_benchmarks=True
cmake --preset conan-default
cmake --build --preset conan-release -j
```
```sh
OMP_NUM_THREADS=1 ./build/benchmarks/Release/benchmark-hyh \
    --benchmark_out=hyh.json --benchmark_repetitions=5 --benchmark_min_time=0.02s \
    --benchmark_enable_random_interleaving
```
```sh
OMP_NUM_THREADS=1 ./build/benchmarks/Release/benchmark-ocp \
    --benchmark_out=ocp.json --benchmark_repetitions=5 --benchmark_min_time=1000x
```

The `desktop` profile enables AVX-512. If this is not supported by your hardware,
you can use the `laptop` profile, which uses AVX2 only (for Intel Skylake and
newer).

OpenBLAS can be used instead of the Intel MKL by passing the options
`-o \&:with_openblas=True -o \&:with_mkl=False` to Conan. Be sure to use CMake's
`--fresh` flag to reconfigure the project after making this change.

Only libstdc++ version 14 is currently supported (GCC 14 or Clang 18).
