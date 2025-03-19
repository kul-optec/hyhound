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
template <index_t R>
using diag_simd_t = optimal_simd_type_t<R, native_simd_size>;
#else
template <index_t R>
using diag_simd_t = optimal_simd_type_t<R, native_simd_size / 2>;
#endif
#else
template <index_t R>
using diag_simd_t = optimal_simd_type_t<R>;
#endif
template <Config Conf>
using tail_simd_L_t = optimal_simd_type_t<Conf.block_size_r>;
template <Config Conf>
using tail_simd_A_t = optimal_simd_type_t<Conf.block_size_s>;

/// Ensures that the matrix W is aligned for SIMD.
template <index_t R = MaxSizeR>
constexpr size_t W_align = stdx::memory_alignment_v<
    stdx::simd<real_t, stdx::simd_abi::deduce_t<real_t, R>>>;

/// Ensures that the first element of every column of W is aligned for SIMD.
template <index_t R = MaxSizeR>
constexpr size_t W_stride = (R * sizeof(real_t) + W_align<R> - 1) / W_align<R> *
                            W_align<R> / sizeof(real_t);
/// Size of the matrix W.
template <index_t R = MaxSizeR>
constexpr size_t W_size = (W_stride<R> * R * sizeof(real_t) + W_align<R> - 1) /
                          W_align<R> * W_align<R> / sizeof(real_t);

template <index_t R = MaxSizeR>
using mut_W_accessor =
    mat_access_impl<real_t, std::integral_constant<index_t, W_stride<R>>>;

template <index_t R = MaxSizeR>
struct matrix_W_storage {
    alignas(W_align<R>) real_t W[W_stride<R> * R]{};
    constexpr operator mut_W_accessor<R>() { return {W}; }
};

template <index_t R, class UpDown>
void updowndate_diag(index_t colsA, mut_W_accessor<> W, real_t *Ld, index_t ldL,
                     real_t *Ad, index_t ldA, UpDownArg<UpDown> signs) noexcept;

template <index_t R, class UpDown>
void updowndate_full(index_t colsA, real_t *Ld, index_t ldL, real_t *Ad,
                     index_t ldA, UpDownArg<UpDown> signs) noexcept;

template <Config Conf, class UpDown>
void updowndate_tail(index_t colsA0, index_t colsA, mut_W_accessor<> W,
                     real_t *Lp, index_t ldL, const real_t *Bp, index_t ldB,
                     real_t *Ap, index_t ldA, UpDownArg<UpDown> signs) noexcept;

} // namespace hyhound::micro_kernels::householder
