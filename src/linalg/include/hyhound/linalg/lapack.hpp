/**
 * @file
 * LAPACK error handling.
 */

#pragma once

#include <hyhound/config.hpp>

#if HYHOUND_WITH_MKL
#include <mkl.h>
#else
#include <lapack.h>
#endif

#include <stdexcept>
#include <type_traits>

namespace hyhound::linalg {

struct lapack_error : std::runtime_error {
    lapack_error(const std::string &name, index_t info)
        : std::runtime_error("LAPACK error: " + name +
                             " (info=" + std::to_string(info) + ")"),
          info {info} {}
    index_t info;
};

template <class Name>
inline void lapack_throw_on_err(Name &&name, index_t info) {
    if (info)
        throw lapack_error(std::forward<Name>(name), info);
}

static_assert(std::is_same_v<index_t, lapack_int>, "Unsupported index type");
static_assert(std::is_signed_v<index_t>);

} // namespace hyhound::linalg

#if !HYHOUND_WITH_MKL
#define dpotrf(...) LAPACK_dpotrf(__VA_ARGS__)
#define spotrf(...) LAPACK_spotrf(__VA_ARGS__)
#define dtrtri(...) LAPACK_dtrtri(__VA_ARGS__)
#define strtri(...) LAPACK_strtri(__VA_ARGS__)
#define dlauum(...) LAPACK_dlauum(__VA_ARGS__)
#define slauum(...) LAPACK_slauum(__VA_ARGS__)
#define dtrtrs(...) LAPACK_dtrtrs(__VA_ARGS__)
#define strtrs(...) LAPACK_strtrs(__VA_ARGS__)
#define dsytrf_rk(...) LAPACK_dsytrf_rk(__VA_ARGS__)
#define ssytrf_rk(...) LAPACK_ssytrf_rk(__VA_ARGS__)
#endif
