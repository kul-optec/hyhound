#pragma once

#include <hyhound/config.hpp>
#include <hyhound/ocp/export.h>
#include <Eigen/Core>

#include <algorithm>
#include <random>

#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

using mat = Eigen::MatrixX<real_t>;
using vec = Eigen::VectorX<real_t>;

template <class T>
auto vw(T &&t) {
    return as_view(std::forward<T>(t), guanaqo::with_index_type<index_t>);
}

struct HYHOUND_OCP_EXPORT OCPDataRiccati {
    index_t N  = 31;
    index_t nx = 40, nu = 10, ny = 10;
    mat BAs  = mat::Zero(nx, (nu + nx) * N);
    mat DCs  = mat::Zero(ny, (nu + nx) * (N + 1));
    mat RSQs = mat::Zero(nu + nx, (nu + nx) * (N + 1));
    mat rqs  = mat::Zero(nu + nx, N + 1);
    mat es   = mat::Zero(nx, N + 1);

    auto F(index_t j) { return BAs.middleCols((nu + nx) * j, nu + nx); }
    auto G(index_t j) { return DCs.middleCols((nu + nx) * j, nu + nx); }
    auto B(index_t j) { return BAs.middleCols((nu + nx) * j, nu); }
    auto A(index_t j) { return BAs.middleCols((nu + nx) * j + nu, nx); }
    auto D(index_t j) { return DCs.middleCols((nu + nx) * j, nu); }
    auto C(index_t j) { return DCs.middleCols((nu + nx) * j + nu, nx); }
    auto H(index_t j) { return RSQs.middleCols((nu + nx) * j, nu + nx); }
    auto Q(index_t j) { return H(j).bottomRightCorner(nx, nx); }
    auto R(index_t j) { return H(j).topLeftCorner(nu, nu); }
    auto S(index_t j) { return H(j).bottomLeftCorner(nx, nu); }
    auto q(index_t j) { return rqs.col(j).bottomRows(nx); }
    auto r(index_t j) { return rqs.col(j).topRows(nu); }
    auto F(index_t j) const { return BAs.middleCols((nu + nx) * j, nu + nx); }
    auto G(index_t j) const { return DCs.middleCols((nu + nx) * j, nu + nx); }
    auto B(index_t j) const { return BAs.middleCols((nu + nx) * j, nu); }
    auto A(index_t j) const { return BAs.middleCols((nu + nx) * j + nu, nx); }
    auto D(index_t j) const { return DCs.middleCols((nu + nx) * j, nu); }
    auto C(index_t j) const { return DCs.middleCols((nu + nx) * j + nu, nx); }
    auto H(index_t j) const { return RSQs.middleCols((nu + nx) * j, nu + nx); }
    auto Q(index_t j) const { return H(j).bottomRightCorner(nx, nx); }
    auto R(index_t j) const { return H(j).topLeftCorner(nu, nu); }
    auto S(index_t j) const { return H(j).bottomLeftCorner(nx, nu); }
    auto q(index_t j) const { return rqs.col(j).bottomRows(nx); }
    auto r(index_t j) const { return rqs.col(j).topRows(nu); }

    void init_random(uint_fast32_t seed = 0) {
        std::mt19937 rng{seed};
        std::uniform_real_distribution<real_t> uni{-1, 1};
        for (index_t i = 0; i < N; ++i) {
            auto Ai = A(i), Bi = B(i), Ci = C(i), Di = D(i);
            auto Hi = H(i);
            std::ranges::generate(Ai.reshaped(), [&] { return uni(rng) / 2; });
            std::ranges::generate(Bi.reshaped(), [&] { return uni(rng); });
            std::ranges::generate(Ci.reshaped(), [&] { return uni(rng); });
            std::ranges::generate(Di.reshaped(), [&] { return uni(rng); });
            std::ranges::generate(Hi.reshaped(), [&] { return uni(rng); });
            auto nux = nu + nx;
            Hi += static_cast<real_t>(nux) * mat::Identity(nux, nux);
            Hi.triangularView<Eigen::StrictlyUpper>() =
                Hi.triangularView<Eigen::StrictlyLower>().transpose();
        }
        auto Ci = C(N);
        auto Qi = Q(N);
        std::ranges::generate(Ci.reshaped(), [&] { return uni(rng); });
        std::ranges::generate(Qi.reshaped(), [&] { return uni(rng); });
        Qi += static_cast<real_t>(nx) * mat::Identity(nx, nx);
        Qi.triangularView<Eigen::StrictlyUpper>() =
            Qi.triangularView<Eigen::StrictlyLower>().transpose();
        // Random right-hand sides
        std::ranges::generate(rqs.reshaped(), [&] { return uni(rng); });
        std::ranges::generate(es.reshaped(), [&] { return uni(rng); });
    }
};

struct HYHOUND_OCP_EXPORT RiccatiFactor {
    const OCPDataRiccati &ocp;
    mat Ls = mat::Zero(ocp.nu + ocp.nx, (ocp.nu + ocp.nx) * (ocp.N + 1));
    mat Vᵀ = mat::Zero(ocp.nx, ocp.nu + ocp.nx);
    mat ΣG = mat::Zero(ocp.ny, ocp.nu + ocp.nx);
    mat p  = mat::Zero(ocp.nx, ocp.N + 1);
    mat Pe = mat::Zero(ocp.nx, ocp.N + 1);
    mat ux = mat::Zero(ocp.nu + ocp.nx, ocp.N + 1);
    mat λ  = mat::Zero(ocp.nx, ocp.N + 1);
    mat ΦY = mat::Zero(ocp.nu + ocp.nx, ocp.ny *(ocp.N + 1));
    mat YΦ = mat::Zero(ocp.nu + ocp.nx, ocp.ny *(ocp.N + 1));
    mat S  = mat::Zero(1, ocp.ny *(ocp.N + 1));

    auto L(index_t j) {
        return Ls.middleCols((ocp.nu + ocp.nx) * j, ocp.nu + ocp.nx);
    }
    auto Lxx(index_t j) { return L(j).bottomRightCorner(ocp.nx, ocp.nx); }
    auto u(index_t j) { return ux.col(j).topRows(ocp.nu); }
    auto x(index_t j) { return ux.col(j).bottomRows(ocp.nx); }
};

HYHOUND_OCP_EXPORT void factor(RiccatiFactor &factor, Eigen::Ref<const mat> Σ);
HYHOUND_OCP_EXPORT void update(RiccatiFactor &factor, Eigen::Ref<const mat> ΔΣ);
HYHOUND_OCP_EXPORT void solve(RiccatiFactor &factor);

} // namespace hyhound::ocp
