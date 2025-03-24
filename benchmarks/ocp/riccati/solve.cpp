#include "ocp/riccati.hpp"

#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

void solve(RiccatiFactor &factor) {
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
        guanaqo::blas::xtrmv_LTN(as_view(Lxx_next, use_index_t),
                                 as_view(Pe_next, use_index_t));
        // Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁)
        guanaqo::blas::xtrmv_LNN(as_view(Lxx_next, use_index_t),
                                 as_view(Pe_next, use_index_t));
        // Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁
        Pe_next += factor.p.col(k + 1);
        // rₖ + Bₖᵀ (Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁)
        u = ocp.r(k);
        guanaqo::blas::xgemv_T(real_t{1}, as_view(B, use_index_t),
                               as_view(Pe_next, use_index_t), real_t{1},
                               as_view(u, use_index_t));
        // Luu⁻¹ₖ (rₖ + Bₖᵀ (Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁))
        guanaqo::blas::xtrsv_LNN(as_view(Luu, use_index_t),
                                 as_view(u, use_index_t));

        // pₖ = qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁) - Lxuᵀ lₖ
        // ------------------------------------------
        auto p = factor.p.col(k);
        p      = ocp.q(k);
        // qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁)
        guanaqo::blas::xgemv_T(real_t{1}, as_view(A, use_index_t),
                               as_view(Pe_next, use_index_t), real_t{1},
                               as_view(p, use_index_t));
        // qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁) - Lxuᵀ lₖ
        guanaqo::blas::xgemv_N(real_t{-1}, as_view(Lxu, use_index_t),
                               as_view(u, use_index_t), real_t{1},
                               as_view(p, use_index_t));
    }
    {
        auto Lxx        = factor.Lxx(0);
        auto x0         = factor.x(0);
        x0              = ocp.es.col(0);
        factor.λ.col(0) = x0;
        guanaqo::blas::xtrmv_LTN(as_view(Lxx, use_index_t),
                                 as_view(factor.λ.col(0), use_index_t));
        guanaqo::blas::xtrmv_LNN(as_view(Lxx, use_index_t),
                                 as_view(factor.λ.col(0), use_index_t));
        factor.λ.col(0) += factor.p.col(0);
    }
    for (index_t k = 0; k < ocp.N; ++k) {
        auto L        = factor.L(k);
        auto Luu      = L.topLeftCorner(ocp.nu, ocp.nu),
             Lxu      = L.bottomLeftCorner(ocp.nx, ocp.nu);
        auto Lxx_next = factor.Lxx(k + 1);
        auto u = factor.u(k), x = factor.x(k), x_next = factor.x(k + 1);
        guanaqo::blas::xgemv_T(real_t{-1}, as_view(Lxu, use_index_t),
                               as_view(x, use_index_t), real_t{-1},
                               as_view(u, use_index_t));
        guanaqo::blas::xtrsv_LTN(as_view(Luu, use_index_t),
                                 as_view(u, use_index_t));
        x_next  = ocp.es.col(k + 1);
        auto BA = ocp.F(k);
        guanaqo::blas::xgemv_N(real_t{1}, as_view(BA, use_index_t),
                               as_view(factor.ux.col(k), use_index_t),
                               real_t{1}, as_view(x_next, use_index_t));
        factor.λ.col(k + 1) = x_next;
        guanaqo::blas::xtrmv_LTN(as_view(Lxx_next, use_index_t),
                                 as_view(factor.λ.col(k + 1), use_index_t));
        guanaqo::blas::xtrmv_LNN(as_view(Lxx_next, use_index_t),
                                 as_view(factor.λ.col(k + 1), use_index_t));
        factor.λ.col(k + 1) += factor.p.col(k + 1);
    }
}

} // namespace hyhound::ocp
