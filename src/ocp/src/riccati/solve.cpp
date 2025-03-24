#include <hyhound/ocp/riccati.hpp>

#include <guanaqo/blas/hl-blas-interface.hpp>

namespace hyhound::ocp {

void solve(RiccatiFactor &factor) {
    using namespace guanaqo::blas;
    const auto &ocp     = factor.ocp;
    factor.p.col(ocp.N) = ocp.q(ocp.N);
    for (index_t k = ocp.N; k-- > 0;) {
        auto B = ocp.B(k), A = ocp.A(k);
        auto L        = factor.L(k);
        auto Luu      = L.topLeftCorner(ocp.nu, ocp.nu),
             Lxu      = L.bottomLeftCorner(ocp.nx, ocp.nu);
        auto Lxx_next = factor.Lxx(k + 1);
        auto u        = factor.u(k);
        auto Pe_next  = factor.Pe.col(k + 1);

        // lₖ = Luu⁻¹ₖ (rₖ + Bₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁))
        // -----------------------------------------
        Pe_next = ocp.es.col(k + 1);
        // Lxxₖ₊₁ᵀ eₖ₊₁
        xtrmv_LTN(vw(Lxx_next), vw(Pe_next));
        // Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁)
        xtrmv_LNN(vw(Lxx_next), vw(Pe_next));
        // Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁
        Pe_next += factor.p.col(k + 1);
        // rₖ + Bₖᵀ (Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁)
        u = ocp.r(k);
        xgemv_T(real_t{1}, vw(B), vw(Pe_next), real_t{1}, vw(u));
        // Luu⁻¹ₖ (rₖ + Bₖᵀ (Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁))
        xtrsv_LNN(vw(Luu), vw(u));

        // pₖ = qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁) - Lxuᵀ lₖ
        // ------------------------------------------
        auto p = factor.p.col(k);
        p      = ocp.q(k);
        // qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁)
        xgemv_T(real_t{1}, vw(A), vw(Pe_next), real_t{1}, vw(p));
        // qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁) - Lxuᵀ lₖ
        xgemv_N(real_t{-1}, vw(Lxu), vw(u), real_t{1}, vw(p));
    }
    {
        auto Lxx        = factor.Lxx(0);
        auto x0         = factor.x(0);
        x0              = ocp.es.col(0);
        factor.λ.col(0) = x0;
        xtrmv_LTN(vw(Lxx), vw(factor.λ.col(0)));
        xtrmv_LNN(vw(Lxx), vw(factor.λ.col(0)));
        factor.λ.col(0) += factor.p.col(0);
    }
    for (index_t k = 0; k < ocp.N; ++k) {
        auto L        = factor.L(k);
        auto Luu      = L.topLeftCorner(ocp.nu, ocp.nu),
             Lxu      = L.bottomLeftCorner(ocp.nx, ocp.nu);
        auto Lxx_next = factor.Lxx(k + 1);
        auto u = factor.u(k), x = factor.x(k), x_next = factor.x(k + 1);
        xgemv_T(real_t{-1}, vw(Lxu), vw(x), real_t{-1}, vw(u));
        xtrsv_LTN(vw(Luu), vw(u));
        x_next  = ocp.es.col(k + 1);
        auto BA = ocp.F(k);
        xgemv_N(real_t{1}, vw(BA), vw(factor.ux.col(k)), real_t{1}, vw(x_next));
        factor.λ.col(k + 1) = x_next;
        xtrmv_LTN(vw(Lxx_next), vw(factor.λ.col(k + 1)));
        xtrmv_LNN(vw(Lxx_next), vw(factor.λ.col(k + 1)));
        factor.λ.col(k + 1) += factor.p.col(k + 1);
    }
}

} // namespace hyhound::ocp
