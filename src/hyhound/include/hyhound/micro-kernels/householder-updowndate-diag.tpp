#pragma once

#include "householder-updowndate.hpp"

#define UNROLL_FOR(...) HYHOUND_FULLY_UNROLLED_IVDEP_FOR (__VA_ARGS__)
#define UNROLL_FOR_A_COLS(...) for (__VA_ARGS__)

namespace hyhound::micro_kernels::householder {

template <index_t R, class T, class UpDown>
[[gnu::hot]] void updowndate_diag(index_t colsA, mut_W_accessor<T> W,
                                  T *__restrict Ld, index_t ldL,
                                  T *__restrict Ad, index_t ldA,
                                  UpDownArg<UpDown> signs) noexcept {
    using std::copysign;
    using std::sqrt;
    using simd = diag_simd_t<T, R>;
    const mut_matrix_accessor<T> L{Ld, ldL}, A{Ad, ldA};
    static constexpr index_t N = simd::size();
    static_assert(R % N == 0);
    static_assert(R > 0);
    static_assert(W.outer_stride >= R);
    HYHOUND_ASSUME(colsA > 0);
    [[maybe_unused]] const auto W_addr = reinterpret_cast<uintptr_t>(W.data);
    HYHOUND_ASSUME((W_addr % W_align<T, R>) == 0);
    static_assert(W_align<T> % W_align<T, R> == 0);
    static_assert(R <= MaxSizeR);

    UNROLL_FOR (index_t k = 0; k < R; ++k) {
        // Compute all inner products between A and a
        simd bb[R / N]{};
        UNROLL_FOR_A_COLS (index_t j = 0; j < colsA; ++j) {
            T Akj = signs(A(k, j), j);
            UNROLL_FOR (index_t i = 0; i < R; i += N) {
                auto Aij = A.template load<simd>(i, j);
                if constexpr (signs.negate)
                    bb[i / N] -= Aij * Akj;
                else
                    bb[i / N] += Aij * Akj;
            }
        }
        // Energy condition and Householder coefficients
        const T α2 = bb[k / N][k % N], Lkk = L(k, k);
        const T L̃kk = copysign(sqrt(Lkk * Lkk + α2), Lkk), β = Lkk + L̃kk;
        const T γoβ = 2 * β / (β * β + α2), γ = β * γoβ, inv_β = 1 / β;
        L(k, k) = L̃kk;
        // Compute L̃
        const index_t kp1N = (k + 1 + N - 1) / N;
        UNROLL_FOR (index_t i = k + 1; i < kp1N * N; ++i) { // scalar part
            auto Lik         = L(i, k);
            bb[i / N][i % N] = γ * Lik + bb[i / N][i % N] * γoβ;
            L(i, k)          = bb[i / N][i % N] - Lik;
        }
        UNROLL_FOR (index_t i = kp1N; i < R / N; ++i) { // vectorized
            index_t ii = i * N;
            auto Lik   = L.template load<simd>(ii, k);
            bb[i]      = γ * Lik + bb[i] * γoβ;
            Lik        = bb[i] - Lik;
            L.store(Lik, ii, k);
        }
        // Update A
        UNROLL_FOR_A_COLS (index_t j = 0; j < colsA; ++j) {
            T Akj = A(k, j) *= inv_β; // Scale Householder vector
            UNROLL_FOR (index_t i = k + 1; i < kp1N * N; ++i) // scalar part
                A(i, j) -= bb[i / N][i % N] * Akj;
            UNROLL_FOR (index_t i = kp1N; i < R / N; ++i) { // vectorized
                index_t ii = i * N;
                auto Aij   = A.template load<simd>(ii, j);
                Aij -= bb[i] * Akj;
                A.store(Aij, ii, j);
            }
        }
        // Save block Householder matrix W
        static constexpr auto alignW = stdx::vector_aligned;
        UNROLL_FOR (index_t i = 0; i < k + 1; i += N)
            bb[i / N] *= inv_β;
        bb[k / N][k % N] = γ; // inverse of diagonal
        UNROLL_FOR (index_t i = 0; i < k + 1; i += N)
            W.store(bb[i / N], i, k, alignW);
    }
}

} // namespace hyhound::micro_kernels::householder

#undef UNROLL_FOR
#undef UNROLL_FOR_A_COLS
