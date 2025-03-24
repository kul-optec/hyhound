#include <hyhound/ocp/schur.hpp>

#include <guanaqo/blas/hl-blas-interface.hpp>

namespace hyhound::ocp {

void solve(SchurFactor &factor) {
    using namespace guanaqo::blas;
    const auto &ocp = factor.ocp;
    auto LN         = factor.LHxx(ocp.N);
    auto vN         = factor.v.col(ocp.N).topRows(ocp.nx);
    vN              = ocp.q(ocp.N);
    auto λN         = factor.λ.col(ocp.N);
    // L v = g
    xtrsv_LNN(vw(LN), vw(vN));
    // e + Wv
    auto WᵀN = factor.Wᵀ(ocp.N).topRows(ocp.nx);
    λN       = ocp.es.col(ocp.N);
    xgemv_T(real_t{1}, vw(WᵀN), vw(vN), real_t{1}, vw(λN));
    for (index_t k = ocp.N; k-- > 0;) {
        auto vj       = factor.v.col(k);
        auto λj       = factor.λ.col(k);
        auto Lj       = factor.LH(k);
        auto Vj       = factor.V(k);
        auto λ_next   = factor.λ.col(k + 1);
        auto LΨd_next = factor.LΨd(k + 1);
        vj            = ocp.qrs.col(k);
        // L v = g
        xtrsv_LNN(vw(Lj), vw(vj));
        // e + Wv
        auto Wᵀj = factor.Wᵀ(k);
        λj       = ocp.es.col(k);
        xgemv_T(real_t{1}, vw(Wᵀj), vw(vj), real_t{1}, vw(λj));
        // e(j+1) + Wv(j+1) - Vv(j)
        xgemv_N(real_t{-1}, vw(Vj), vw(vj), real_t{1}, vw(λ_next));
        if (k + 1 < ocp.N) {
            auto LΨs_next    = factor.LΨs(k + 1);
            auto λ_next_next = factor.λ.col(k + 2);
            xgemv_N(real_t{-1}, vw(LΨs_next), vw(λ_next_next), real_t{1},
                    vw(λ_next));
        }
        // LΨd(j+1) λ(j+1) = e(j+1) + Wv(j+1) - Vv(j) - LΨs(j+1) λ(k+2)
        xtrsv_LNN(vw(LΨd_next), vw(λ_next));
    }
    {
        index_t k = 0;
        auto λ0   = factor.λ.col(k);
        auto LΨd0 = factor.LΨd(k);
        if (k < ocp.N) {
            auto LΨs    = factor.LΨs(k);
            auto λ_next = factor.λ.col(k + 1);
            xgemv_N(real_t{-1}, vw(LΨs), vw(λ_next), real_t{1}, vw(λ0));
        }
        // LΨd(0) λ(0) = e(0) + Wv(0) - LΨs(0) λ(1)
        xtrsv_LNN(vw(LΨd0), vw(λ0));
        // LΨd(0)ᵀ λ(0) = λ(0)
        xtrsv_LTN(vw(LΨd0), vw(λ0));
    }
    for (index_t k = 1; k <= ocp.N; ++k) {
        auto λj       = factor.λ.col(k);
        auto λ_prev   = factor.λ.col(k - 1);
        auto LΨdj     = factor.LΨd(k);
        auto LΨs_prev = factor.LΨs(k - 1);
        // λ(j) - LΨs(j-1)ᵀ λ(j-1)
        xgemv_T(real_t{-1}, vw(LΨs_prev), vw(λ_prev), real_t{1}, vw(λj));
        // LΨd(j)ᵀ λ(j) = λ(j)
        xtrsv_LTN(vw(LΨdj), vw(λj));
    }
    for (index_t k = 0; k < ocp.N; ++k) {
        // v(0) - W(0)ᵀ λ(0) + Vᵀ λ(1)
        auto Lj      = factor.LH(k);
        auto Wᵀj     = factor.Wᵀ(k);
        auto Vj      = factor.V(k);
        auto λj      = factor.λ.col(k);
        auto λj_next = factor.λ.col(k + 1);
        auto vj      = factor.v.col(k);
        xgemv_N(real_t{1}, vw(Wᵀj), vw(λj), real_t{-1}, vw(vj));
        xgemv_T(real_t{-1}, vw(Vj), vw(λj_next), real_t{1}, vw(vj));
        // L d = v
        xtrsv_LTN(vw(Lj), vw(vj));
    }
    {
        index_t k = ocp.N;
        // v(N) - W(N)ᵀ λ(N)
        auto Lj  = factor.LHxx(k);
        auto Wᵀj = factor.Wᵀ(k).topRows(ocp.nx);
        auto λj  = factor.λ.col(k);
        auto vj  = factor.v.col(k).topRows(ocp.nx);
        xgemv_N(real_t{1}, vw(Wᵀj), vw(λj), real_t{-1}, vw(vj));
        // L d = v
        xtrsv_LTN(vw(Lj), vw(vj));
    }
}

} // namespace hyhound::ocp
