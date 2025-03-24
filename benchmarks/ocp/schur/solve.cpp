#include "ocp/schur.hpp"

#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

void solve(SchurFactor &factor) {
    const auto &ocp = factor.ocp;
    auto LN         = factor.LHxx(ocp.N);
    auto vN         = factor.v.col(ocp.N).topRows(ocp.nx);
    vN              = ocp.q(ocp.N);
    auto λN         = factor.λ.col(ocp.N);
    // L v = g
    guanaqo::blas::xtrsv_LNN(as_view(LN, use_index_t),
                             as_view(vN, use_index_t));
    // e + Wv
    auto WᵀN = factor.Wᵀ(ocp.N).topRows(ocp.nx);
    λN       = ocp.es.col(ocp.N);
    guanaqo::blas::xgemv_T(real_t{1}, as_view(WᵀN, use_index_t),
                           as_view(vN, use_index_t), real_t{1},
                           as_view(λN, use_index_t));
    for (index_t k = ocp.N; k-- > 0;) {
        auto vj       = factor.v.col(k);
        auto λj       = factor.λ.col(k);
        auto Lj       = factor.LH(k);
        auto Vj       = factor.V(k);
        auto λ_next   = factor.λ.col(k + 1);
        auto LΨd_next = factor.LΨd(k + 1);
        vj            = ocp.qrs.col(k);
        // L v = g
        guanaqo::blas::xtrsv_LNN(as_view(Lj, use_index_t),
                                 as_view(vj, use_index_t));
        // e + Wv
        auto Wᵀj = factor.Wᵀ(k);
        λj       = ocp.es.col(k);
        guanaqo::blas::xgemv_T(real_t{1}, as_view(Wᵀj, use_index_t),
                               as_view(vj, use_index_t), real_t{1},
                               as_view(λj, use_index_t));
        // e(j+1) + Wv(j+1) - Vv(j)
        guanaqo::blas::xgemv_N(real_t{-1}, as_view(Vj, use_index_t),
                               as_view(vj, use_index_t), real_t{1},
                               as_view(λ_next, use_index_t));
        if (k + 1 < ocp.N) {
            auto LΨs_next    = factor.LΨs(k + 1);
            auto λ_next_next = factor.λ.col(k + 2);
            guanaqo::blas::xgemv_N(real_t{-1}, as_view(LΨs_next, use_index_t),
                                   as_view(λ_next_next, use_index_t), real_t{1},
                                   as_view(λ_next, use_index_t));
        }
        // LΨd(j+1) λ(j+1) = e(j+1) + Wv(j+1) - Vv(j) - LΨs(j+1) λ(k+2)
        guanaqo::blas::xtrsv_LNN(as_view(LΨd_next, use_index_t),
                                 as_view(λ_next, use_index_t));
    }
    {
        index_t k = 0;
        auto λ0   = factor.λ.col(k);
        auto LΨd0 = factor.LΨd(k);
        if (k < ocp.N) {
            auto LΨs    = factor.LΨs(k);
            auto λ_next = factor.λ.col(k + 1);
            guanaqo::blas::xgemv<real_t, index_t>(
                CblasColMajor, CblasNoTrans, LΨs.rows(), LΨs.cols(), real_t{-1},
                LΨs.data(), LΨs.outerStride(), λ_next.data(), index_t{1},
                real_t{1}, λ0.data(), index_t{1});
        }
        // LΨd(0) λ(0) = e(0) + Wv(0) - LΨs(0) λ(1)
        guanaqo::blas::xtrsv<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, λ0.rows(),
            LΨd0.data(), LΨd0.outerStride(), λ0.data(), index_t{1});
        // LΨd(0)ᵀ λ(0) = λ(0)
        guanaqo::blas::xtrsv<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, λ0.rows(),
            LΨd0.data(), LΨd0.outerStride(), λ0.data(), index_t{1});
    }
    for (index_t k = 1; k <= ocp.N; ++k) {
        auto λj       = factor.λ.col(k);
        auto λ_prev   = factor.λ.col(k - 1);
        auto LΨdj     = factor.LΨd(k);
        auto LΨs_prev = factor.LΨs(k - 1);
        // λ(j) - LΨs(j-1)ᵀ λ(j-1)
        guanaqo::blas::xgemv_T(real_t{-1}, as_view(LΨs_prev, use_index_t),
                               as_view(λ_prev, use_index_t), real_t{1},
                               as_view(λj, use_index_t));
        // LΨd(j)ᵀ λ(j) = λ(j)
        guanaqo::blas::xtrsv_LTN(as_view(LΨdj, use_index_t),
                                 as_view(λj, use_index_t));
    }
    for (index_t k = 0; k < ocp.N; ++k) {
        // v(0) - W(0)ᵀ λ(0) + Vᵀ λ(1)
        auto Lj      = factor.LH(k);
        auto Wᵀj     = factor.Wᵀ(k);
        auto Vj      = factor.V(k);
        auto λj      = factor.λ.col(k);
        auto λj_next = factor.λ.col(k + 1);
        auto vj      = factor.v.col(k);
        guanaqo::blas::xgemv_N(real_t{1}, as_view(Wᵀj, use_index_t),
                               as_view(λj, use_index_t), real_t{-1},
                               as_view(vj, use_index_t));
        guanaqo::blas::xgemv_T(real_t{-1}, as_view(Vj, use_index_t),
                               as_view(λj_next, use_index_t), real_t{1},
                               as_view(vj, use_index_t));
        // L d = v
        guanaqo::blas::xtrsv_LTN(as_view(Lj, use_index_t),
                                 as_view(vj, use_index_t));
    }
    {
        index_t k = ocp.N;
        // v(N) - W(N)ᵀ λ(N)
        auto Lj  = factor.LHxx(k);
        auto Wᵀj = factor.Wᵀ(k).topRows(ocp.nx);
        auto λj  = factor.λ.col(k);
        auto vj  = factor.v.col(k).topRows(ocp.nx);
        guanaqo::blas::xgemv_N(real_t{1}, as_view(Wᵀj, use_index_t),
                               as_view(λj, use_index_t), real_t{-1},
                               as_view(vj, use_index_t));
        // L d = v
        guanaqo::blas::xtrsv_LTN(as_view(Lj, use_index_t),
                                 as_view(vj, use_index_t));
    }
}

} // namespace hyhound::ocp
