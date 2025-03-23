#include "ocp/schur.hpp"

#include <guanaqo/blas/hl-blas-interface.hpp>

namespace hyhound::ocp {

void factor(SchurFactor &factor, Eigen::Ref<const mat> Σ) {
    const auto &ocp = factor.ocp;
    auto CN         = ocp.C(ocp.N);
    auto ΣCN        = factor.ΣG.rightCols(ocp.nx);
    auto LxxN       = factor.LHxx(ocp.N);
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
    auto WxN                           = factor.Wᵀ(ocp.N).topRows(ocp.nx);
    WxN.triangularView<Eigen::Lower>() = LxxN.triangularView<Eigen::Lower>();
    guanaqo::blas::xtrtri<real_t, index_t>("L", "N", WxN.rows(), WxN.data(),
                                           WxN.outerStride(), &info);
    if (info != 0)
        throw std::runtime_error(std::format(
            "Inverse failed at stage {} with status {}", ocp.N, info));
    auto LΨdN                           = factor.LΨd(ocp.N);
    LΨdN.triangularView<Eigen::Lower>() = WxN.triangularView<Eigen::Lower>();
    guanaqo::blas::xlauum<real_t, index_t>("L", LΨdN.rows(), LΨdN.data(),
                                           LΨdN.outerStride(), &info);

    for (index_t j = ocp.N; j-- > 0;) {
        auto Gj       = ocp.G(j);
        auto &ΣGj     = factor.ΣG;
        auto Lj       = factor.LH(j);
        auto LΨd_next = factor.LΨd(j + 1);
        auto LΨdj     = factor.LΨd(j);
        auto LΨsj     = factor.LΨs(j);
        // H(j) = Hl(j) + G(j)ᵀΣ(j) G(j)
        ΣGj = Σ.col(j).asDiagonal() * Gj;
        Lj.triangularView<Eigen::Lower>() =
            ocp.H(j).triangularView<Eigen::Lower>();
        guanaqo::blas::xgemmt<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, Lj.rows(),
            ocp.ny, real_t{1}, Gj.data(), Gj.outerStride(), ΣGj.data(),
            ΣGj.outerStride(), real_t{1}, Lj.data(), Lj.outerStride());
        guanaqo::blas::xpotrf<real_t, index_t>("L", Lj.rows(), Lj.data(),
                                               Lj.outerStride(), &info);
        if (info != 0)
            throw std::runtime_error(std::format(
                "Cholesky failed at stage {} with status {}", j, info));
        // V = F L⁻ᵀ
        auto Vj = factor.V(j);
        Vj      = ocp.F(j);
        guanaqo::blas::xtrsm<real_t, index_t>(
            CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
            Vj.rows(), Vj.cols(), real_t{1}, Lj.data(), Lj.outerStride(),
            Vj.data(), Vj.outerStride());
        // VVᵀ
        guanaqo::blas::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, LΨd_next.rows(), Vj.cols(),
            real_t{1}, Vj.data(), Vj.outerStride(), real_t{1}, LΨd_next.data(),
            LΨd_next.outerStride());
        // chol(Θ + VVᵀ)
        guanaqo::blas::xpotrf<real_t, index_t>("L", LΨd_next.rows(),
                                               LΨd_next.data(),
                                               LΨd_next.outerStride(), &info);
        if (info != 0)
            throw std::runtime_error(std::format(
                "Cholesky Ψ failed at stage {} with status {}", j + 1, info));
        // W = (I 0) L⁻ᵀ
        auto Wᵀj = factor.Wᵀ(j);
        Wᵀj      = Lj.topLeftCorner(ocp.nx + ocp.nu, ocp.nx);
        guanaqo::blas::xtrtri<real_t, index_t>("L", "N", ocp.nx, Wᵀj.data(),
                                               Wᵀj.outerStride(), &info);
        if (info != 0)
            throw std::runtime_error(std::format(
                "Inverse failed at stage {} with status {}", j, info));
        auto Wᵀxj = Wᵀj.topRows(ocp.nx), Wᵀuj = Wᵀj.bottomRows(ocp.nu);
        guanaqo::blas::xtrmm<real_t, index_t>(
            CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit,
            Wᵀuj.rows(), Wᵀuj.cols(), real_t(-1), Wᵀxj.data(),
            Wᵀxj.outerStride(), Wᵀuj.data(), Wᵀuj.outerStride());
        auto Luuj = Lj.bottomRightCorner(ocp.nu, ocp.nu);
        guanaqo::blas::xtrsm<real_t, index_t>(
            CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
            Wᵀuj.rows(), Wᵀuj.cols(), real_t{1}, Luuj.data(),
            Luuj.outerStride(), Wᵀuj.data(), Wᵀuj.outerStride());
        // -WVᵀ
        guanaqo::blas::xgemm<real_t, index_t>(
            CblasColMajor, CblasTrans, CblasTrans, LΨsj.rows(), LΨsj.cols(),
            Wᵀj.rows(), real_t{-1}, Wᵀj.data(), Wᵀj.outerStride(), Vj.data(),
            Vj.outerStride(), real_t{0}, LΨsj.data(), LΨsj.outerStride());
        guanaqo::blas::xtrsm<real_t, index_t>(
            CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
            LΨsj.rows(), LΨsj.cols(), real_t{1}, LΨd_next.data(),
            LΨd_next.outerStride(), LΨsj.data(), LΨsj.outerStride());
        // WWᵀ
        guanaqo::blas::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, LΨdj.rows(), Wᵀj.rows(),
            real_t{1}, Wᵀj.data(), Wᵀj.outerStride(), real_t{0}, LΨdj.data(),
            LΨdj.outerStride());
        // -LΨs LΨsᵀ
        guanaqo::blas::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, LΨdj.rows(), LΨsj.cols(),
            real_t{-1}, LΨsj.data(), LΨsj.outerStride(), real_t{1}, LΨdj.data(),
            LΨdj.outerStride());
    }
    // chol(Θ)
    auto LΨd0 = factor.LΨd(0);
    guanaqo::blas::xpotrf<real_t, index_t>("L", LΨd0.rows(), LΨd0.data(),
                                           LΨd0.outerStride(), &info);
    if (info != 0)
        throw std::runtime_error(std::format(
            "Cholesky Ψ failed at stage {} with status {}", 0, info));
}

} // namespace hyhound::ocp
