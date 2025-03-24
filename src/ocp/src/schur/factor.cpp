#include <hyhound/ocp/schur.hpp>

#include <guanaqo/blas/hl-blas-interface.hpp>

namespace hyhound::ocp {

void factor(SchurFactor &factor, Eigen::Ref<const mat> Σ) {
    using namespace guanaqo::blas;
    const auto &ocp = factor.ocp;
    auto CN         = ocp.C(ocp.N);
    auto ΣCN        = factor.ΣG.rightCols(ocp.nx);
    auto LxxN       = factor.LHxx(ocp.N);
    // H(N) = Q(N) + C(N)ᵀΣ(N) G(N)
    ΣCN = Σ.col(ocp.N).asDiagonal() * CN;
    LxxN.triangularView<Eigen::Lower>() =
        ocp.Q(ocp.N).triangularView<Eigen::Lower>();
    xgemmt_LTN(real_t{1}, vw(CN), vw(ΣCN), real_t{1}, vw(LxxN));
    // Lxx(N) = chol(H(N))
    xpotrf_L(vw(LxxN));
    auto WxN                           = factor.Wᵀ(ocp.N).topRows(ocp.nx);
    WxN.triangularView<Eigen::Lower>() = LxxN.triangularView<Eigen::Lower>();
    xtrtri_LN(vw(WxN));
    auto LΨdN                           = factor.LΨd(ocp.N);
    LΨdN.triangularView<Eigen::Lower>() = WxN.triangularView<Eigen::Lower>();
    xlauum_L(vw(LΨdN));

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
        xgemmt_LTN(real_t{1}, vw(Gj), vw(ΣGj), real_t{1}, vw(Lj));
        xpotrf_L(vw(Lj));
        // V = F L⁻ᵀ
        auto Vj = factor.V(j);
        Vj      = ocp.F(j);
        xtrsm_RLTN(real_t{1}, vw(Lj), vw(Vj));
        // VVᵀ
        xsyrk_LN(real_t{1}, vw(Vj), //
                 real_t{1}, vw(LΨd_next));
        // chol(Θ + VVᵀ)
        xpotrf_L(vw(LΨd_next));
        // W = (I 0) L⁻ᵀ
        auto Wᵀj  = factor.Wᵀ(j);
        Wᵀj       = Lj.topLeftCorner(ocp.nx + ocp.nu, ocp.nx);
        auto Wᵀxj = Wᵀj.topRows(ocp.nx), Wᵀuj = Wᵀj.bottomRows(ocp.nu);
        xtrtri_LN(vw(Wᵀxj));
        xtrmm_RLNN(real_t{-1}, vw(Wᵀxj), vw(Wᵀuj));
        auto Luuj = Lj.bottomRightCorner(ocp.nu, ocp.nu);
        xtrsm_LLNN(real_t{1}, vw(Luuj), vw(Wᵀuj));
        // -WVᵀ
        xgemm_TT(real_t{-1}, vw(Wᵀj), vw(Vj), real_t{0}, vw(LΨsj));
        xtrsm_RLTN(real_t{1}, vw(LΨd_next), vw(LΨsj));
        // WWᵀ
        xsyrk_LT(real_t{1}, vw(Wᵀj), //
                 real_t{0}, vw(LΨdj));
        // -LΨs LΨsᵀ
        xsyrk_LN(real_t{-1}, vw(LΨsj), //
                 real_t{1}, vw(LΨdj));
    }
    // chol(Θ)
    xpotrf_L(vw(factor.LΨd(0)));
}

} // namespace hyhound::ocp
