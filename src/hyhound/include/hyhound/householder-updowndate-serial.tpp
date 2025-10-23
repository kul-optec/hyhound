#pragma once

#include <hyhound/householder-updowndate-micro-kernels.tpp>
#include <hyhound/householder-updowndate.hpp>
#include <hyhound/loop.hpp>
#include <hyhound/lut.hpp>

namespace hyhound::inline serial {

template <class T, Config<T> Conf, class UpDown>
void update_cholesky(MatrixView<T> L, MatrixView<T> A, UpDown signs,
                     MatrixView<T> Ws) {
    static constexpr index_t R = Conf.block_size_r, S = Conf.block_size_s;
    static constexpr index_t N       = Conf.num_blocks_r;
    static constexpr bool do_packing = Conf.enable_packing;
    using Rt                         = index_constant<R>;
    using St                         = index_constant<S>;
    static constexpr micro_kernels::householder::Config uConf{
        .block_size_r        = R,
        .block_size_s        = S,
        .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
    static constexpr micro_kernels::householder::Config uConfR{
        .block_size_r        = R,
        .block_size_s        = R,
        .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
    assert(L.rows >= L.cols);
    assert(L.rows == A.rows);
    assert(Ws.rows == 0 || Ws.cols == L.cols);
    constinit static auto full_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            return updowndate_full<NR + 1, T, UpDown>;
        });
    constinit static auto diag_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            return updowndate_diag<NR + 1, T, UpDown>;
        });
    constinit static auto tail_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            constexpr micro_kernels::householder::Config uConf{
                .block_size_r        = NR + 1,
                .block_size_s        = S,
                .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
            return updowndate_tile_tail<uConf, T, UpDown>;
        });

    // Leaner accessors (without unnecessary dimensions and strides).
    micro_kernels::mut_matrix_accessor<T> L_{L}, A_{A}, Ws_{Ws};
    // Workspace storage for W (upper triangular Householder representation)
    micro_kernels::householder::matrix_W_storage<T> W[N];

    // Optional packing of one block row of A.
    auto A_pack_storage = [&] {
        if constexpr (do_packing) {
            index_t num_pack = R * A.cols * N;
            return std::vector<T>(num_pack);
        } else {
            struct Empty {};
            return Empty{};
        }
    }();
    T *A_pack[N];
    if constexpr (do_packing)
        for (index_t i = 0; i < N; ++i)
            A_pack[i] = &A_pack_storage[R * A.cols * i];
    auto pack_Ad = [&](index_t k,
                       index_t nk =
                           Rt{}) -> micro_kernels::mut_matrix_accessor<T> {
        if constexpr (do_packing) {
            MatrixView<T> Ad{
                {.data = A_pack[(k / R) % N], .rows = nk, .cols = A.cols}};
            Ad = A.middle_rows(k, nk);
            return Ad;
        }
        return A.middle_rows(k, nk);
    };
    auto unpack_Ad = [&](index_t k, index_t nk = Rt{}) {
        if constexpr (do_packing) {
            MatrixView<const T> Ad{
                {.data = A_pack[(k / R) % N], .rows = nk, .cols = A.cols}};
            A.middle_rows(k, nk) = Ad;
        }
    };

    // Process all diagonal blocks (in multiples of NR, except the last).
    index_t k;
    for (k = 0; k + R * N <= L.cols; k += R * N) {
        micro_kernels::mut_matrix_accessor<T> Adk[N];
        // Process all rows in the diagonal block (in multiples of R)
        for (index_t kk = 0; kk < R * N; kk += R) {
            // Pack the part of A corresponding to this diagonal block
            auto &Ad = Adk[kk / R] = pack_Ad(k + kk);
            // Process blocks left of the diagonal
            for (index_t cc = 0; cc < kk; cc += R) {
                auto Ls = L_.block(k + kk, k + cc);
                updowndate_tail<uConfR, T, UpDown>(0, A.cols, W[cc / R], Ls,
                                                   Adk[cc / R], Ad, signs);
            }
            auto Ld = L_.block(k + kk, k + kk);
            // Process the diagonal block itself
            updowndate_diag<R, T, UpDown>(A.cols, W[kk / R], Ld, Ad, signs);
            // Store W if requested
            if (Ws.rows >= R) {
                for (index_t c = 0; c < R; ++c)
                    for (index_t r = 0; r <= c; ++r)
                        Ws_(r, k + kk + c) = W[kk / R](r, c);
                unpack_Ad(k + kk);
            }
        }
        // Process all rows below the diagonal block (in multiples of S).
        foreach_chunked(
            k + R * N, L.rows, St{},
            [&](index_t i) {
                auto As = A_.middle_rows(i);
                // Process columns
                for (index_t cc = 0; cc < R * N; cc += R) {
                    auto Ls = L_.block(i, k + cc);
                    for (index_t c = 0; c < R; ++c)
                        __builtin_prefetch(&Ls(0, c), 0, 0); // non-temporal
                    updowndate_tail<uConf, T, UpDown>(0, A.cols, W[cc / R], Ls,
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
                    updowndate_tile_tail<uConf, T, UpDown>(
                        rem_i, 0, A.cols, W[cc / R], Ls, Adk[cc / R], As,
                        signs);
                }
            },
            LoopDir::Forward);
    }
    index_t rem_k = L.cols - k;
    if (rem_k > 0) {
        if (N != 1)
            throw std::logic_error("Not yet implemented");
        assert(rem_k < R);
        auto Ad = pack_Ad(k, rem_k);
        auto Ld = L_.block(k, k);
        if (L.rows == L.cols && Ws.rows < R) {
            full_microkernel_lut[rem_k - 1](A.cols, Ld, Ad, signs);
        } else {
            diag_microkernel_lut[rem_k - 1](A.cols, W[0], Ld, Ad, signs);
            // Store W if requested
            if (Ws.rows >= R) {
                for (index_t c = 0; c < rem_k; ++c)
                    for (index_t r = 0; r <= c; ++r)
                        Ws_(r, k + c) = W[0](r, c);
                unpack_Ad(k, rem_k);
            }
            // Process all rows below the diagonal block (in multiples of S).
            foreach_chunked_merged(
                k + rem_k, L.rows, St{},
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
