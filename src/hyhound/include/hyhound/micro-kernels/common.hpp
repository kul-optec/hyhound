#pragma once

#include <hyhound/assume.hpp>
#include <hyhound/cneg.hpp>
#include <hyhound/config.hpp>
#include <hyhound/matrix-view.hpp>
#include <hyhound/unroll.h>
#include <hyhound/updown.hpp>

#include <experimental/simd>
#include <type_traits>

#include "matrix-accessor.hpp"

namespace hyhound::micro_kernels {

namespace stdx                     = std::experimental;
using native_abi                   = stdx::simd_abi::native<real_t>;
constexpr index_t native_simd_size = stdx::simd_size_v<real_t, native_abi>;

template <index_t Bs, index_t MaxVecLen = 0>
struct optimal_simd_type {
    static constexpr index_t max_vec_len =
        MaxVecLen > 0 ? MaxVecLen : native_simd_size;
    static constexpr index_t vec_len =
        ((Bs > max_vec_len) && (Bs % max_vec_len == 0)) ? max_vec_len : Bs;
    using simd_abi = stdx::simd_abi::deduce_t<real_t, vec_len>;
    using type     = stdx::simd<real_t, simd_abi>;
};
template <index_t Bs, index_t MaxVecLen = 0>
using optimal_simd_type_t = typename optimal_simd_type<Bs, MaxVecLen>::type;

template <class UpDown>
struct UpDownArg;

template <>
struct UpDownArg<Update> {
    UpDownArg(Update) {}
    auto operator()(auto x, index_t) const { return x; }
    static constexpr bool negate = false;
};

template <>
struct UpDownArg<Downdate> {
    UpDownArg(Downdate) {}
    auto operator()(auto x, index_t) const { return x; }
    static constexpr bool negate = true;
};

template <>
struct UpDownArg<UpDowndate> {
    UpDownArg(UpDowndate ud) : signs{ud.signs.data()} {}
    const real_t *__restrict signs;
    template <class T>
    auto operator()(T x, index_t j) const {
        return hyhound::cneg(x, T{signs[j]});
    }
    static constexpr bool negate = false;
};

template <>
struct UpDownArg<DownUpdate> {
    UpDownArg(DownUpdate du) : signs{du.signs.data()} {}
    const real_t *__restrict signs;
    template <class T>
    auto operator()(T x, index_t j) const {
        return hyhound::cneg(x, T{signs[j]});
    }
    static constexpr bool negate = true;
};

template <>
struct UpDownArg<DiagonalUpDowndate> {
    UpDownArg(DiagonalUpDowndate dg) : diag{dg.diag.data()} {}
    const real_t *__restrict diag;
    template <class T>
    auto operator()(T x, index_t j) const {
        return x * T{diag[j]};
    }
    static constexpr bool negate = false;
};

} // namespace hyhound::micro_kernels
