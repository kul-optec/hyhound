#pragma once

#include <hyhound/config.hpp>
#include <hyhound/stringify.h>

#if HYHOUND_WITH_OPENMP
#include <omp.h>
#define HYHOUND_OMP(X) _Pragma(HYHOUND_STRINGIFY(omp X))
#define HYHOUND_OMP_IF_ELSE(X, Y) X
#define HYHOUND_OMP_IF(X) X
#else
#define HYHOUND_OMP(X)
#define HYHOUND_OMP_IF_ELSE(X, Y) Y
#define HYHOUND_OMP_IF(X)                                                     \
    do {                                                                       \
    } while (0)
#endif
