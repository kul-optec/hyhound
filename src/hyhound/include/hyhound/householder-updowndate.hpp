#pragma once

#include <hyhound/config.hpp>
#include <hyhound/updown.hpp>

#include <guanaqo/mat-view.hpp>

#include <type_traits>

namespace hyhound {

namespace detail {
template <class T>
struct DefaultMicroKernelSizes;
template <>
struct DefaultMicroKernelSizes<float> {
#ifdef __AVX512F__
    // AVX512 has 32 vector registers:
    static constexpr index_t DefaultSizeR = 8;
    static constexpr index_t DefaultSizeS = 32;
#elif defined(__ARM_NEON)
    // NEON has 32 vector registers:
    static constexpr index_t DefaultSizeR = 4; // TODO: tune
    static constexpr index_t DefaultSizeS = 32;
#else
    // AVX2 has 16 vector registers:
    static constexpr index_t DefaultSizeR = 4;
    static constexpr index_t DefaultSizeS = 32;
#endif
};
template <>
struct DefaultMicroKernelSizes<double> {
#ifdef __AVX512F__
    // AVX512 has 32 vector registers:
    static constexpr index_t DefaultSizeR = 8;
    static constexpr index_t DefaultSizeS = 24;
#elif defined(__ARM_NEON)
    // NEON has 32 vector registers:
    static constexpr index_t DefaultSizeR = 4; // TODO: tune
    static constexpr index_t DefaultSizeS = 12;
#else
    // AVX2 has 16 vector registers:
    static constexpr index_t DefaultSizeR = 4;
    static constexpr index_t DefaultSizeS = 12;
#endif
};
} // namespace detail

template <class T>
struct Config {
    /// Block size of the block column of L to process in the micro-kernels.
    index_t block_size_r = detail::DefaultMicroKernelSizes<T>::DefaultSizeR;
    /// Block size of the block row of L to process in the micro-kernels.
    index_t block_size_s = detail::DefaultMicroKernelSizes<T>::DefaultSizeS;
    /// Number of block columns per cache block.
    index_t num_blocks_r = 1;
    /// Column prefetch distance for the matrix A.
    index_t prefetch_dist_col_a = 4;
    /// Enable cache blocking by copying the current block row of A to a
    /// temporary buffer.
    bool enable_packing = true;
};

template <class T = real_t>
using MatrixView = guanaqo::MatrixView<T, index_t>;

inline namespace serial {
template <class T, Config<T> Conf = {}, class UpDown>
void update_cholesky(MatrixView<T> L, MatrixView<T> A, UpDown signs,
                     MatrixView<T> Ws = MatrixView<T>{{.rows = 0}});
template <class T, Config<T> Conf = {}, class UpDown>
void apply_householder(MatrixView<T> L, MatrixView<T> A, UpDown signs,
                       std::type_identity_t<MatrixView<const T>> Ws,
                       std::type_identity_t<MatrixView<const T>> B);
} // namespace serial

} // namespace hyhound
