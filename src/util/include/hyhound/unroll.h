#pragma once

#include <guanaqo/stringify.h>

#ifdef NDEBUG

#ifdef __clang__
#define HYHOUND_FULLY_UNROLL_LOOP _Pragma("clang loop unroll(full)")
#define HYHOUND_FULLY_UNROLLED_FOR(...)                                        \
    HYHOUND_FULLY_UNROLL_LOOP for (__VA_ARGS__)
#define HYHOUND_FULLY_UNROLLED_IVDEP_FOR(...)                                  \
    HYHOUND_FULLY_UNROLL_LOOP _Pragma(                                         \
        "clang loop vectorize(enable) interleave(enable)") for (__VA_ARGS__)
#define HYHOUND_UNROLLED_IVDEP_FOR(N, ...)                                     \
    _Pragma(GUANAQO_STRINGIFY(clang loop unroll_count(N))) _Pragma(            \
        "clang loop vectorize(enable) interleave(enable)") for (__VA_ARGS__)
#else
#define HYHOUND_FULLY_UNROLL_LOOP _Pragma("GCC unroll 99")
#define HYHOUND_FULLY_UNROLLED_FOR(...)                                        \
    HYHOUND_FULLY_UNROLL_LOOP for (__VA_ARGS__)
#define HYHOUND_FULLY_UNROLLED_IVDEP_FOR(...)                                  \
    HYHOUND_FULLY_UNROLL_LOOP _Pragma("GCC ivdep") for (__VA_ARGS__)
#define HYHOUND_UNROLLED_IVDEP_FOR(N, ...)                                     \
    _Pragma(GUANAQO_STRINGIFY(GCC unroll N))                                   \
        _Pragma("GCC ivdep") for (__VA_ARGS__)
#endif

#else

#define HYHOUND_FULLY_UNROLL_LOOP
#define HYHOUND_FULLY_UNROLLED_FOR(...) for (__VA_ARGS__)
#define HYHOUND_FULLY_UNROLLED_IVDEP_FOR(...) for (__VA_ARGS__)
#define HYHOUND_UNROLLED_IVDEP_FOR(N, ...) for (__VA_ARGS__)

#endif
