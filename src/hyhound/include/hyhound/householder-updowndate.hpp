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

/// Blocking and packing parameters for the @ref update_cholesky and
/// @ref apply_householder functions.
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

/// Non-owning view of a dense matrix
/// (column-major storage order with unit inner stride).
/// @tparam T   The matrix element type.
template <class T = real_t>
using MatrixView = guanaqo::MatrixView<T, index_t>;

inline namespace serial {
/**
 * @brief Update the Cholesky factorization @f$ LL^\top = H @f$ such that the
 * updated factorization satisfies
 * @f$ \tilde L \tilde L^\top = H + A S A^\top @f$.
 *
 * @f$ S @f$ is a diagonal matrix with the diagonal elements determined by the
 * @p signs parameter:
 * - If @p signs is of type @ref Update, then @f$ S = \mathrm{I} @f$;
 * - If @p signs is of type @ref Downdate, then @f$ S = -\mathrm{I} @f$;
 * - If @p signs is of type @ref UpDowndate, then @f$ S = \mathrm{diag}(s) @f$,
 *   where @f$ s = \mathrm{copysign}(1, \texttt{signs.signs}) @f$, and where
 *   `signs.signs` contains only the values `+0.0` and `-0.0`;
 * - If @p signs is of type @ref DownUpdate, then @f$ S = -\mathrm{diag}(s) @f$,
 *   where @f$ s = \mathrm{copysign}(1, \texttt{signs.signs}) @f$, and where
 *   `signs.signs` contains only the values `+0.0` and `-0.0`;
 * - If @p signs is of type @ref DiagonalUpDowndate, then
 *   @f$ S = \mathrm{diag}(\texttt{signs.diag}) @f$.
 *
 * If @f$ L @f$ is a tall matrix, the function also returns a matrix
 * @f$ \tilde A @f$ such that @f$ \tilde L \tilde L^\top +
 * \begin{pmatrix} 0 \\ \tilde A \end{pmatrix} S
 * \begin{pmatrix} 0 & \tilde A^\top \end{pmatrix} = L L^\top + A S A^\top @f$.
 *
 * @param[in,out] L (k × n, lower trapezoidal)
 *      On input, the lower-triangular Cholesky factor @f$ L @f$ of @f$ H @f$.
 *      On output, the lower-triangular Cholesky factor @f$ \tilde L @f$ of
 *      @f$ H + A S A^\top @f$.
 * @param[in,out] A (k × m)
 *      On input, the matrix @f$ A @f$.
 *      On output, the top n rows are overwritten, and the bottom k-n rows
 *      contain the matrix @f$ \tilde A @f$.
 * @param[in] signs
 *      The update/downdate specification. See above for details.
 * @param[out] Ws (r × n)
 *      If specified, this matrix will contain part of the block Householder
 *      representation that was applied to perform the update/downdate, and the
 *      top n rows of A will contain the corresponding  Householder reflectors.
 *      The number of rows r depends on the block size specified by
 *      @ref Config::block_size_r.
 */
template <class T, Config<T> Conf = {}, class UpDown>
void update_cholesky(MatrixView<T> L, MatrixView<T> A, UpDown signs,
                     MatrixView<T> Ws = {{}});

/**
 * @brief Apply a block Householder transformation @f$ \breve Q @f$,
 * returning @f$ \begin{pmatrix} \tilde L & \tilde A \end{pmatrix} = 
 * \begin{pmatrix} L & A \end{pmatrix} \breve Q @f$.
 *
 * @param[in,out] L (l × n)
 *      On input, the matrix @f$ L @f$.
 *      On output, the matrix @f$ \tilde L @f$.
 * @param[in,out] A (l × m)
 *      On input, the matrix @f$ A @f$.
 *      On output, the matrix @f$ \tilde A @f$.
 * @param[in] signs
 *      The update/downdate specification. See @ref update_cholesky for details.
 * @param[in] B (n × m)
 *      The Householder reflectors generated by @ref update_cholesky.
 * @param[in] Ws (r × n)
 *      The upper triangular block Householder representation generated by
 *      @ref update_cholesky.
 */
template <class T, Config<T> Conf = {}, class UpDown>
void apply_householder(MatrixView<T> L, MatrixView<T> A, UpDown signs,
                       std::type_identity_t<MatrixView<const T>> B,
                       std::type_identity_t<MatrixView<const T>> Ws);
} // namespace serial

} // namespace hyhound
