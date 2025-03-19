#include "ocp/riccati.hpp"

#include <hyhound/linalg/blas-interface.hpp>

namespace hyhound::ocp {

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
        linalg::xtrmv<real_t, index_t>(CblasColMajor, CblasLower, CblasTrans,
                                       CblasNonUnit, Lxx_next.rows(),
                                       Lxx_next.data(), Lxx_next.outerStride(),
                                       Pe_next.data(), index_t{1});
        // Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁)
        linalg::xtrmv<real_t, index_t>(CblasColMajor, CblasLower, CblasNoTrans,
                                       CblasNonUnit, Lxx_next.rows(),
                                       Lxx_next.data(), Lxx_next.outerStride(),
                                       Pe_next.data(), index_t{1});
        // Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁
        Pe_next += factor.p.col(k + 1);
        // rₖ + Bₖᵀ (Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁)
        u = ocp.r(k);
        linalg::xgemv<real_t, index_t>(
            CblasColMajor, CblasTrans, B.rows(), B.cols(), real_t{1}, B.data(),
            B.outerStride(), Pe_next.data(), index_t{1}, real_t{1}, u.data(),
            index_t{1});
        // Luu⁻¹ₖ (rₖ + Bₖᵀ (Lxxₖ₊₁ (Lxxₖ₊₁ᵀ eₖ₊₁) + pₖ₊₁))
        linalg::xtrsv<real_t, index_t>(CblasColMajor, CblasLower, CblasNoTrans,
                                       CblasNonUnit, Luu.rows(), Luu.data(),
                                       Luu.outerStride(), u.data(), index_t{1});

        // pₖ = qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁) - Lxuᵀ lₖ
        // ------------------------------------------
        auto p = factor.p.col(k);
        p      = ocp.q(k);
        // qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁)
        linalg::xgemv<real_t, index_t>(
            CblasColMajor, CblasTrans, A.rows(), A.cols(), real_t{1}, A.data(),
            A.outerStride(), Pe_next.data(), index_t{1}, real_t{1}, p.data(),
            index_t{1});
        // qₖ + Aₖᵀ (Pₖ₊₁ eₖ₊₁ + pₖ₊₁) - Lxuᵀ lₖ
        linalg::xgemv<real_t, index_t>(CblasColMajor, CblasNoTrans, Lxu.rows(),
                                       Lxu.cols(), real_t{-1}, Lxu.data(),
                                       Lxu.outerStride(), u.data(), index_t{1},
                                       real_t{1}, p.data(), index_t{1});
    }
    {
        auto Lxx        = factor.Lxx(0);
        auto x0         = factor.x(0);
        x0              = ocp.es.col(0);
        factor.λ.col(0) = x0;
        linalg::xtrmv<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, Lxx.rows(),
            Lxx.data(), Lxx.outerStride(), factor.λ.col(0).data(), index_t{1});
        linalg::xtrmv<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, Lxx.rows(),
            Lxx.data(), Lxx.outerStride(), factor.λ.col(0).data(), index_t{1});
        factor.λ.col(0) += factor.p.col(0);
    }
    for (index_t k = 0; k < ocp.N; ++k) {
        auto L        = factor.L(k);
        auto Luu      = L.topLeftCorner(ocp.nu, ocp.nu),
             Lxu      = L.bottomLeftCorner(ocp.nx, ocp.nu);
        auto Lxx_next = factor.Lxx(k + 1);
        auto u = factor.u(k), x = factor.x(k), x_next = factor.x(k + 1);
        linalg::xgemv<real_t, index_t>(CblasColMajor, CblasTrans, Lxu.rows(),
                                       Lxu.cols(), real_t{-1}, Lxu.data(),
                                       Lxu.outerStride(), x.data(), index_t{1},
                                       real_t{-1}, u.data(), index_t{1});
        linalg::xtrsv<real_t, index_t>(CblasColMajor, CblasLower, CblasTrans,
                                       CblasNonUnit, Luu.rows(), Luu.data(),
                                       Luu.outerStride(), u.data(), index_t{1});
        x_next  = ocp.es.col(k + 1);
        auto BA = ocp.F(k);
        linalg::xgemv<real_t, index_t>(
            CblasColMajor, CblasNoTrans, BA.rows(), BA.cols(), real_t{1},
            BA.data(), BA.outerStride(), factor.ux.col(k).data(), index_t{1},
            real_t{1}, x_next.data(), index_t{1});
        factor.λ.col(k + 1) = x_next;
        linalg::xtrmv<real_t, index_t>(CblasColMajor, CblasLower, CblasTrans,
                                       CblasNonUnit, Lxx_next.rows(),
                                       Lxx_next.data(), Lxx_next.outerStride(),
                                       factor.λ.col(k + 1).data(), index_t{1});
        linalg::xtrmv<real_t, index_t>(CblasColMajor, CblasLower, CblasNoTrans,
                                       CblasNonUnit, Lxx_next.rows(),
                                       Lxx_next.data(), Lxx_next.outerStride(),
                                       factor.λ.col(k + 1).data(), index_t{1});
        factor.λ.col(k + 1) += factor.p.col(k + 1);
    }
}

} // namespace hyhound::ocp
