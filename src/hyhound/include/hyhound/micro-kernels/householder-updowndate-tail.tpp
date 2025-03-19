#pragma once

#include "householder-updowndate.hpp"

#define UNROLL_FOR(...) HYHOUND_FULLY_UNROLLED_IVDEP_FOR (__VA_ARGS__)
#define UNROLL_FOR_A_COLS(...) for (__VA_ARGS__)

namespace hyhound::micro_kernels::householder {

template <Config Conf, class UpDown>
[[gnu::hot]] void
updowndate_tail(index_t colsA0, index_t colsA, mut_W_accessor<> W,
                real_t *__restrict Lp, index_t ldL, const real_t *__restrict Bp,
                index_t ldB, real_t *__restrict Ap, index_t ldA,
                UpDownArg<UpDown> signs) noexcept {
    using simdA = tail_simd_A_t<Conf>;
    using simdL = tail_simd_L_t<Conf>;
    const mut_matrix_accessor L{Lp, ldL}, A{Ap, ldA};
    const matrix_accessor B{Bp, ldB};
    static constexpr index_t NA = simdA::size(), NL = simdL::size();
    static constexpr index_t S = Conf.block_size_s, R = Conf.block_size_r;
    static_assert(S % NA == 0);
    static_assert(R % NL == 0);
    static_assert(R > 0);
    static_assert(W.outer_stride >= R);
    HYHOUND_ASSUME(colsA0 >= 0);
    HYHOUND_ASSUME(colsA >= colsA0);
    [[maybe_unused]] const auto W_addr = reinterpret_cast<uintptr_t>(W.data);
    HYHOUND_ASSUME(W_addr % W_align<Conf.block_size_r> == 0);

    // Compute product U = A B
    simdA V[S / NA][R]{};
    UNROLL_FOR_A_COLS (index_t j = colsA0; j < colsA; ++j) {
        if (Conf.prefetch_dist_col_a > 0)
            _mm_prefetch(&A(0, j + Conf.prefetch_dist_col_a), _MM_HINT_T0);
        UNROLL_FOR (index_t kk = 0; kk < R; kk += NL) {
            auto Akj = signs(B.load<simdL>(kk, j), j);
            UNROLL_FOR (index_t k = 0; k < NL; ++k)
                UNROLL_FOR (index_t i = 0; i < S; i += NA)
                    if constexpr (signs.negate)
                        V[i / NA][kk + k] -= A.load<simdA>(i, j) * Akj[k];
                    else
                        V[i / NA][kk + k] += A.load<simdA>(i, j) * Akj[k];
        }
    }
    // Solve system V = (L+U)W⁻¹ (in-place)
    UNROLL_FOR (index_t k = 0; k < R; ++k) {
        simdL Wk[R / NL];
        UNROLL_FOR (index_t l = 0; l < k; l += NL)
            Wk[l / NL] = W.template load<simdL>(l, k, stdx::vector_aligned);
        UNROLL_FOR (index_t i = 0; i < S; i += NA) {
            auto Lik = L.load<simdA>(i, k);
            V[i / NA][k] += Lik;
            UNROLL_FOR (index_t l = 0; l < k; ++l)
                V[i / NA][k] -= V[i / NA][l] * Wk[l / NL][l % NL];
            V[i / NA][k] *= W(k, k); // diagonal already inverted
            Lik = V[i / NA][k] - Lik;
            L.store(Lik, i, k);
        }
    }
    // Update A = -V Bᵀ
    UNROLL_FOR_A_COLS (index_t j = 0; j < colsA0; ++j) {
        simdL Akj[R / NL];
        UNROLL_FOR (index_t kk = 0; kk < R; kk += NL)
            Akj[kk / NL] = B.load<simdL>(kk, j);
        UNROLL_FOR (index_t i = 0; i < S; i += NA) {
            simdA Aij{0};
            UNROLL_FOR (index_t kk = 0; kk < R; kk += NL) {
                UNROLL_FOR (index_t k = 0; k < NL; ++k)
                    Aij -= V[i / NA][kk + k] * Akj[kk / NL][k];
            }
            A.store(Aij, i, j);
        }
    }
    // Update A -= V Bᵀ
    UNROLL_FOR_A_COLS (index_t j = colsA0; j < colsA; ++j) {
        if (Conf.prefetch_dist_col_a > 0)
            _mm_prefetch(&A(0, j + Conf.prefetch_dist_col_a), _MM_HINT_T0);
        simdL Akj[R / NL];
        UNROLL_FOR (index_t kk = 0; kk < R; kk += NL)
            Akj[kk / NL] = B.load<simdL>(kk, j);
        UNROLL_FOR (index_t i = 0; i < S; i += NA) {
            auto Aij = A.load<simdA>(i, j);
            UNROLL_FOR (index_t kk = 0; kk < R; kk += NL) {
                UNROLL_FOR (index_t k = 0; k < NL; ++k)
                    Aij -= V[i / NA][kk + k] * Akj[kk / NL][k];
            }
            A.store(Aij, i, j);
        }
    }
}

} // namespace hyhound::micro_kernels::householder

#undef UNROLL_FOR
#undef UNROLL_FOR_A_COLS
