#pragma once

#include <hyhound/assume.hpp>
#include <hyhound/cneg.hpp>
#include <hyhound/config.hpp>
#include <hyhound/unroll.h>
#include <hyhound/updown.hpp>

#include <experimental/simd>
#include <type_traits>

#include "matrix-accessor.hpp"

namespace hyhound::micro_kernels {

namespace stdx = std::experimental;
template <class T>
using native_abi = stdx::simd_abi::native<T>;
template <class T>
constexpr index_t native_simd_size = stdx::simd_size_v<T, native_abi<T>>;

template <class T, index_t Bs, index_t MaxVecLen = 0>
struct optimal_simd_type {
    static constexpr index_t max_vec_len =
        MaxVecLen > 0 ? MaxVecLen : native_simd_size<T>;
    static constexpr index_t vec_len =
        ((Bs > max_vec_len) && (Bs % max_vec_len == 0)) ? max_vec_len : Bs;
    using simd_abi = stdx::simd_abi::deduce_t<T, vec_len>;
    using type     = stdx::simd<T, simd_abi>;
};
template <class T, index_t Bs, index_t MaxVecLen = 0>
using optimal_simd_type_t = typename optimal_simd_type<T, Bs, MaxVecLen>::type;

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

template <class T>
struct UpDownArg<UpDowndate<T>> {
    UpDownArg(UpDowndate<T> ud) : signs{ud.signs.data()} {}
    const T *__restrict signs;
    template <class U>
    auto operator()(U x, index_t j) const {
        return hyhound::cneg(x, U{signs[j]});
    }
    static constexpr bool negate = false;
};

template <class T>
struct UpDownArg<DownUpdate<T>> {
    UpDownArg(DownUpdate<T> du) : signs{du.signs.data()} {}
    const T *__restrict signs;
    template <class U>
    auto operator()(U x, index_t j) const {
        return hyhound::cneg(x, U{signs[j]});
    }
    static constexpr bool negate = true;
};

template <class T>
struct UpDownArg<DiagonalUpDowndate<T>> {
    UpDownArg(DiagonalUpDowndate<T> dg) : diag{dg.diag.data()} {}
    const T *__restrict diag;
    template <class U>
    auto operator()(U x, index_t j) const {
        return x * U{diag[j]};
    }
    static constexpr bool negate = false;
};

} // namespace hyhound::micro_kernels
