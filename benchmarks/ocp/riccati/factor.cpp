#include "ocp/riccati.hpp"

#include <guanaqo/blas/hl-blas-interface.hpp>

namespace hyhound::ocp {

void factor(RiccatiFactor &factor, Eigen::Ref<const mat> Σ) {
    const auto &ocp = factor.ocp;
    /* j = N */ {
        auto CN   = ocp.C(ocp.N);
        auto ΣCN  = factor.ΣG.rightCols(ocp.nx);
        auto LxxN = factor.Lxx(ocp.N);
        // H(N) = Q(N) + C(N)ᵀΣ(N) G(N)
        ΣCN = Σ.col(ocp.N).asDiagonal() * CN;
        LxxN.triangularView<Eigen::Lower>() =
            ocp.Q(ocp.N).triangularView<Eigen::Lower>();
        guanaqo::blas::xgemmt<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, LxxN.rows(),
            CN.rows(), real_t{1}, CN.data(), CN.outerStride(), ΣCN.data(),
            ΣCN.outerStride(), real_t{1}, LxxN.data(), LxxN.outerStride());
        // Lxx(N) = chol(H(N))
        index_t info = 0;
        guanaqo::blas::xpotrf<real_t, index_t>("L", LxxN.rows(), LxxN.data(),
                                               LxxN.outerStride(), &info);
        if (info != 0)
            throw std::runtime_error(std::format(
                "Cholesky failed at stage {} with status {}", ocp.N, info));
    }
    for (index_t j = ocp.N; j-- > 0;) {
        auto Gj       = ocp.G(j);
        auto &ΣGj     = factor.ΣG;
        auto Lj       = factor.L(j);
        auto Lxx_next = factor.Lxx(j + 1);
        auto &Vᵀ      = factor.Vᵀ;
        // V(j) = F(j)ᵀ Lxx(j+1), Vᵀ = Lxx(j+1)ᵀ F(j)
        Vᵀ = ocp.F(j);
        guanaqo::blas::xtrmm<real_t, index_t>(
            CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit,
            Vᵀ.rows(), Vᵀ.cols(), real_t{1}, Lxx_next.data(),
            Lxx_next.outerStride(), Vᵀ.data(), Vᵀ.outerStride());
        // H(j) = Hl(j) + G(j)ᵀΣ(j) G(j) + V(j)Vᵀ(j)
        ΣGj = Σ.col(j).asDiagonal() * Gj;
        Lj.triangularView<Eigen::Lower>() =
            ocp.H(j).triangularView<Eigen::Lower>();
        guanaqo::blas::xgemmt<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, Lj.rows(),
            ocp.ny, real_t{1}, Gj.data(), Gj.outerStride(), ΣGj.data(),
            ΣGj.outerStride(), real_t{1}, Lj.data(), Lj.outerStride());
        guanaqo::blas::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, Lj.rows(), Vᵀ.rows(),
            real_t{1}, Vᵀ.data(), Vᵀ.outerStride(), real_t{1}, Lj.data(),
            Lj.outerStride());
        // L(j) = chol(H(j))
        index_t info;
        guanaqo::blas::xpotrf<real_t, index_t>("L", Lj.rows(), Lj.data(),
                                               Lj.outerStride(), &info);
        if (info != 0)
            throw std::runtime_error(std::format(
                "Cholesky failed at stage {} with status {}", j, info));
    }
}

} // namespace hyhound::ocp
