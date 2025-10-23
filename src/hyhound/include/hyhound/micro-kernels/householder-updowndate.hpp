#pragma once

#include "common.hpp"

namespace hyhound::micro_kernels::householder {

struct Config {
    /// Block size of the block column of L to process in the micro-kernels.
    index_t block_size_r;
    /// Block size of the block row of L to process in the micro-kernels.
    index_t block_size_s;
    /// Column prefetch distance for the matrix A.
    index_t prefetch_dist_col_a = 4;
};

constexpr index_t MaxSizeR = 32;

// Since we're dealing with triangular matrices, it pays off to use a smaller
// vector length, because that means we're doing less work: for example, for an
// 8×8 triangular matrix with vector length 8, all 64 elements need to be
// processed, whereas only 48 elements (three quarters) would be processed for
// a vector length of 4.
// Since I'm using an Icelake Intel Core, AVX2 and AVX512 FMA instructions
// have the same throughput (0.5 CPI for 4 elements or 1 CPI for 8). For
// more modern/expensive systems, we may want to always use 8 elements,
// since these CPUs also have two 512-bit FMA units (0.5 CPI), resulting
// in a higher throughput for 8-element vectors.
#if __AVX512F__
#if HYHOUND_HAVE_TWO_512_FMA_UNITS
template <class T, index_t R>
using diag_simd_t = optimal_simd_type_t<T, R, native_simd_size<T>>;
#else
template <class T, index_t R>
using diag_simd_t = optimal_simd_type_t<T, R, native_simd_size<T> / 2>;
#endif
#else
template <class T, index_t R>
using diag_simd_t = optimal_simd_type_t<T, R>;
#endif
template <class T, Config Conf>
using tail_simd_L_t = optimal_simd_type_t<T, Conf.block_size_r>;
template <class T, Config Conf>
using tail_simd_A_t = optimal_simd_type_t<T, Conf.block_size_s>;

/// Ensures that the matrix W is aligned for SIMD.
template <class T, index_t R = MaxSizeR>
constexpr size_t W_align =
    stdx::memory_alignment_v<stdx::simd<T, stdx::simd_abi::deduce_t<T, R>>>;

/// Ensures that the first element of every column of W is aligned for SIMD.
template <class T, index_t R = MaxSizeR>
constexpr size_t W_stride = (R * sizeof(T) + W_align<T, R> - 1) /
                            W_align<T, R> * W_align<T, R> / sizeof(T);
/// Size of the matrix W.
template <class T, index_t R = MaxSizeR>
constexpr size_t W_size = (W_stride<T, R> * R * sizeof(T) + W_align<T, R> - 1) /
                          W_align<T, R> * W_align<T, R> / sizeof(T);

template <class T, index_t R = MaxSizeR>
using mut_W_accessor =
    mat_access_impl<T, std::integral_constant<index_t, W_stride<T, R>>>;

template <class T, index_t R = MaxSizeR>
struct matrix_W_storage {
    alignas(W_align<T, R>) T W[W_stride<T, R> * R]{};
    constexpr operator mut_W_accessor<T, R>() { return {W}; }
    T &operator()(index_t r, index_t c) {
        return mut_W_accessor<T, R>{W}(r, c);
    }
};

template <index_t R, class T, class UpDown>
void updowndate_diag(index_t colsA, mut_W_accessor<T> W, T *Ld, index_t ldL,
                     T *Ad, index_t ldA, UpDownArg<UpDown> signs) noexcept;

template <index_t R, class T, class UpDown>
void updowndate_full(index_t colsA, T *Ld, index_t ldL, T *Ad, index_t ldA,
                     UpDownArg<UpDown> signs) noexcept;

template <Config Conf, class T, class UpDown>
void updowndate_tail(index_t colsA0, index_t colsA, mut_W_accessor<T> W, T *Lp,
                     index_t ldL, const T *Bp, index_t ldB, T *Ap, index_t ldA,
                     UpDownArg<UpDown> signs) noexcept;

} // namespace hyhound::micro_kernels::householder
