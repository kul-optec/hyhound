#pragma once

#include <hyhound/householder-updowndate-micro-kernels.tpp>
#include <hyhound/householder-updowndate.hpp>
#include <hyhound/loop.hpp>
#include <hyhound/lut.hpp>
#include <type_traits>

namespace hyhound::inline serial {

template <Config Conf, class UpDown>
void update_cholesky(MutableRealMatrixView L, MutableRealMatrixView A,
                     UpDown signs) {
    constexpr index_t R = Conf.block_size_r, S = Conf.block_size_s;
    constexpr index_t N       = Conf.num_blocks_r;
    constexpr bool do_packing = Conf.enable_packing;
    constexpr micro_kernels::householder::Config uConf{
        .block_size_r        = R,
        .block_size_s        = S,
        .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
    constexpr micro_kernels::householder::Config uConfR{
        .block_size_r        = R,
        .block_size_s        = R,
        .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
    assert(L.rows >= L.cols);
    assert(L.rows == A.rows);
    constinit static auto full_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            return updowndate_full<NR + 1, UpDown>;
        });
    constinit static auto diag_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            return updowndate_diag<NR + 1, UpDown>;
        });
    constinit static auto tail_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            constexpr micro_kernels::householder::Config uConf{
                .block_size_r        = NR + 1,
                .block_size_s        = S,
                .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
            return updowndate_tile_tail<uConf, UpDown>;
        });

    // Leaner accessors (without unnecessary dimensions and strides).
    micro_kernels::mut_matrix_accessor L_{L}, A_{A};
    // Workspace storage for W (upper triangular Householder representation)
    micro_kernels::householder::matrix_W_storage<> W[N];

    // Optional packing of one block row of A.
    auto A_pack_storage = [&] {
        if constexpr (do_packing) {
            index_t num_pack = R * A.cols * N;
            return std::vector<real_t>(num_pack);
        } else {
            struct Empty {};
            return Empty{};
        }
    }();
    real_t *A_pack[N];
    if constexpr (do_packing)
        for (index_t i = 0; i < N; ++i)
            A_pack[i] = &A_pack_storage[R * A.cols * i];
    auto pack_Ad = [&](index_t k) -> micro_kernels::mut_matrix_accessor {
        if constexpr (do_packing) {
            MutableRealMatrixView Ad{
                {.data = A_pack[(k / R) % N], .rows = R, .cols = A.cols}};
            Ad = A.middle_rows(k, R);
            return Ad;
        }
        return A.middle_rows(k, R);
    };

    // Process all diagonal blocks (in multiples of NR, except the last).
    index_t k;
    for (k = 0; k + R * N <= L.cols; k += R * N) {
        micro_kernels::mut_matrix_accessor Adk[N];
        // Process all rows in the diagonal block (in multiples of R)
        for (index_t kk = 0; kk < R * N; kk += R) {
            // Pack the part of A corresponding to this diagonal block
            auto &Ad = Adk[kk / R] = pack_Ad(k + kk);
            // Process blocks left of the diagonal
            for (index_t cc = 0; cc < kk; cc += R) {
                auto Ls = L_.block(k + kk, k + cc);
                updowndate_tail<uConfR, UpDown>(0, A.cols, W[cc / R], Ls,
                                                Adk[cc / R], Ad, signs);
            }
            auto Ld = L_.block(k + kk, k + kk);
            // Process the diagonal block itself
            updowndate_diag<R, UpDown>(A.cols, W[kk / R], Ld, Ad, signs);
        }
        // Process all rows below the diagonal block (in multiples of S).
        foreach_chunked(
            k + R * N, L.rows, std::integral_constant<index_t, S>(),
            [&](index_t i) {
                auto As = A_.middle_rows(i);
                // Process columns
                for (index_t cc = 0; cc < R * N; cc += R) {
                    auto Ls = L_.block(i, k + cc);
                    for (index_t c = 0; c < R; ++c)
                        __builtin_prefetch(&Ls(0, c), 0, 0); // non-temporal
                    updowndate_tail<uConf, UpDown>(0, A.cols, W[cc / R], Ls,
                                                   Adk[cc / R], As, signs);
                }
            },
            [&](index_t i, index_t rem_i) {
                auto As = A_.middle_rows(i);
                // Process columns
                for (index_t cc = 0; cc < R * N; cc += R) {
                    auto Ls = L_.block(i, k + cc);
                    for (index_t c = 0; c < R; ++c)
                        __builtin_prefetch(&Ls(0, c), 0, 0); // non-temporal
                    updowndate_tile_tail<uConf, UpDown>(rem_i, 0, A.cols,
                                                        W[cc / R], Ls,
                                                        Adk[cc / R], As, signs);
                }
            },
            LoopDir::Forward);
    }
    index_t rem_k = L.cols - k;
    assert(rem_k < R);
    if (rem_k > 0) {
        if (N != 1)
            throw std::logic_error("Not yet implemented");
        auto Ad = A_.middle_rows(k); // TODO: pack?
        auto Ld = L_.block(k, k);
        if (L.rows == L.cols) {
            full_microkernel_lut[rem_k - 1](A.cols, Ld, Ad, signs);
        } else {
            diag_microkernel_lut[rem_k - 1](A.cols, W[0], Ld, Ad, signs);
            // Process all rows below the diagonal block (in multiples of S).
            foreach_chunked_merged(
                k + rem_k, L.rows, std::integral_constant<index_t, S>(),
                [&](index_t i, index_t rem_i) {
                    auto As = A_.middle_rows(i);
                    // Process columns
                    auto Ls = L_.block(i, k);
                    for (index_t c = 0; c < R; ++c)
                        __builtin_prefetch(&Ls(0, c), 0, 0); // non-temporal
                    tail_microkernel_lut[rem_k - 1](rem_i, 0, A.cols, W[0], Ls,
                                                    Ad, As, signs);
                },
                LoopDir::Forward);
        }
    }
}

} // namespace hyhound::inline serial
