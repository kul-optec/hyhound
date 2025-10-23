#pragma once

#include <hyhound/householder-updowndate-micro-kernels.tpp>
#include <hyhound/householder-updowndate.hpp>
#include <hyhound/loop.hpp>
#include <hyhound/lut.hpp>
#include <type_traits>

namespace hyhound::inline serial {

template <class T, Config<T> Conf, class UpDown>
void apply_householder(MatrixView<T> L, MatrixView<T> A, UpDown signs,
                       std::type_identity_t<MatrixView<const T>> Ws,
                       std::type_identity_t<MatrixView<const T>> B) {
    static constexpr index_t R = Conf.block_size_r, S = Conf.block_size_s;
    static constexpr index_t N       = Conf.num_blocks_r;
    static constexpr bool do_packing = Conf.enable_packing;
    using Rt                         = index_constant<R>;
    using St                         = index_constant<S>;
    constexpr micro_kernels::householder::Config uConf{
        .block_size_r        = R,
        .block_size_s        = S,
        .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
    assert(L.rows == A.rows);
    assert(B.rows == L.cols);
    assert(B.cols == A.cols);
    assert(Ws.rows >= R);
    assert(Ws.cols == L.cols);
    constinit static auto tail_microkernel_lut =
        make_1d_lut<R>([]<index_t NR>(index_constant<NR>) {
            constexpr micro_kernels::householder::Config uConf{
                .block_size_r        = NR + 1,
                .block_size_s        = S,
                .prefetch_dist_col_a = Conf.prefetch_dist_col_a};
            return updowndate_tile_tail<uConf, T, UpDown>;
        });

    // Leaner accessors (without unnecessary dimensions and strides).
    micro_kernels::mut_matrix_accessor<T> L_{L}, A_{A};
    micro_kernels::matrix_accessor<T> Ws_{Ws}, B_{B};
    // Workspace storage for W (upper triangular Householder representation)
    micro_kernels::householder::matrix_W_storage<T> W[N]{};

    // Optional packing of one block row of B.
    auto B_pack_storage = [&] {
        if constexpr (do_packing) {
            index_t num_pack = R * B.cols * N;
            return std::vector<T>(num_pack);
        } else {
            struct Empty {};
            return Empty{};
        }
    }();
    T *B_pack[N];
    if constexpr (do_packing)
        for (index_t i = 0; i < N; ++i)
            B_pack[i] = &B_pack_storage[R * B.cols * i];
    auto pack_Bd = [&](index_t k,
                       index_t nk = Rt{}) -> micro_kernels::matrix_accessor<T> {
        if constexpr (do_packing) {
            MatrixView<T> Bd{
                {.data = B_pack[(k / R) % N], .rows = nk, .cols = B.cols}};
            Bd = B.middle_rows(k, nk);
            return Bd;
        }
        return B.middle_rows(k, nk);
    };

    // Process all diagonal blocks (in multiples of NR, except the last).
    index_t k;
    for (k = 0; k + R * N <= L.cols; k += R * N) {
        micro_kernels::matrix_accessor<T> Bdk[N];
        // Process all rows in the diagonal block (in multiples of R)
        for (index_t kk = 0; kk < R * N; kk += R) {
            // Pack the part of B corresponding to this diagonal block
            Bdk[kk / R] = pack_Bd(k + kk);
            // Load W
            for (index_t c = 0; c < R; ++c)
                for (index_t r = 0; r <= c; ++r)
                    W[kk / R](r, c) = Ws_(r, k + kk + c);
        }
        // Process all rows below the diagonal block (in multiples of S).
        foreach_chunked(
            0, L.rows, St{},
            [&](index_t i) {
                auto As = A_.middle_rows(i);
                // Process columns
                for (index_t cc = 0; cc < R * N; cc += R) {
                    auto Ls = L_.block(i, k + cc);
                    for (index_t c = 0; c < R; ++c)
                        __builtin_prefetch(&Ls(0, c), 0, 0); // non-temporal
                    updowndate_tail<uConf, T, UpDown>(0, A.cols, W[cc / R], Ls,
                                                      Bdk[cc / R], As, signs);
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
                        rem_i, 0, A.cols, W[cc / R], Ls, Bdk[cc / R], As,
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
        auto Bd = pack_Bd(k, rem_k);
        // Load W
        for (index_t c = 0; c < rem_k; ++c)
            for (index_t r = 0; r <= c; ++r)
                W[0](r, c) = Ws_(r, k + c);
        // Process all rows below the diagonal block (in multiples of S).
        foreach_chunked_merged(
            0, L.rows, St{},
            [&](index_t i, index_t rem_i) {
                auto As = A_.middle_rows(i);
                // Process columns
                auto Ls = L_.block(i, k);
                for (index_t c = 0; c < R; ++c)
                    __builtin_prefetch(&Ls(0, c), 0, 0); // non-temporal
                tail_microkernel_lut[rem_k - 1](rem_i, 0, A.cols, W[0], Ls, Bd,
                                                As, signs);
            },
            LoopDir::Forward);
    }
}

} // namespace hyhound::inline serial
