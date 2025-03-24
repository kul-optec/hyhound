#include "ocp/schur.hpp"

#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

void factor(SchurFactor &factor, Eigen::Ref<const mat> Σ) {
    const auto &ocp = factor.ocp;
    auto CN         = ocp.C(ocp.N);
    auto ΣCN        = factor.ΣG.rightCols(ocp.nx);
    auto LxxN       = factor.LHxx(ocp.N);
    // H(N) = Q(N) + C(N)ᵀΣ(N) G(N)
    ΣCN = Σ.col(ocp.N).asDiagonal() * CN;
    LxxN.triangularView<Eigen::Lower>() =
        ocp.Q(ocp.N).triangularView<Eigen::Lower>();
    guanaqo::blas::xgemmt_LTN(real_t{1}, as_view(CN, use_index_t),
                              as_view(ΣCN, use_index_t), real_t{1},
                              as_view(LxxN, use_index_t));
    // Lxx(N) = chol(H(N))
    guanaqo::blas::xpotrf_L(as_view(LxxN, use_index_t));
    auto WxN                           = factor.Wᵀ(ocp.N).topRows(ocp.nx);
    WxN.triangularView<Eigen::Lower>() = LxxN.triangularView<Eigen::Lower>();
    guanaqo::blas::xtrtri_LN(as_view(WxN, use_index_t));
    auto LΨdN                           = factor.LΨd(ocp.N);
    LΨdN.triangularView<Eigen::Lower>() = WxN.triangularView<Eigen::Lower>();
    guanaqo::blas::xlauum_L(as_view(LΨdN, use_index_t));

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
        guanaqo::blas::xgemmt_LTN(real_t{1}, as_view(Gj, use_index_t),
                                  as_view(ΣGj, use_index_t), real_t{1},
                                  as_view(Lj, use_index_t));
        guanaqo::blas::xpotrf_L(as_view(Lj, use_index_t));
        // V = F L⁻ᵀ
        auto Vj = factor.V(j);
        Vj      = ocp.F(j);
        guanaqo::blas::xtrsm_RLTN(real_t{1}, as_view(Lj, use_index_t),
                                  as_view(Vj, use_index_t));
        // VVᵀ
        guanaqo::blas::xsyrk_LN(real_t{1}, as_view(Vj, use_index_t), //
                                real_t{1}, as_view(LΨd_next, use_index_t));
        // chol(Θ + VVᵀ)
        guanaqo::blas::xpotrf_L(as_view(LΨd_next, use_index_t));
        // W = (I 0) L⁻ᵀ
        auto Wᵀj  = factor.Wᵀ(j);
        Wᵀj       = Lj.topLeftCorner(ocp.nx + ocp.nu, ocp.nx);
        auto Wᵀxj = Wᵀj.topRows(ocp.nx), Wᵀuj = Wᵀj.bottomRows(ocp.nu);
        guanaqo::blas::xtrtri_LN(as_view(Wᵀxj, use_index_t));
        guanaqo::blas::xtrmm_RLNN(real_t{-1}, as_view(Wᵀxj, use_index_t),
                                  as_view(Wᵀuj, use_index_t));
        auto Luuj = Lj.bottomRightCorner(ocp.nu, ocp.nu);
        guanaqo::blas::xtrsm_LLNN(real_t{1}, as_view(Luuj, use_index_t),
                                  as_view(Wᵀuj, use_index_t));
        // -WVᵀ
        guanaqo::blas::xgemm_TT(real_t{-1}, as_view(Wᵀj, use_index_t),
                                as_view(Vj, use_index_t), real_t{0},
                                as_view(LΨsj, use_index_t));
        guanaqo::blas::xtrsm_RLTN(real_t{1}, as_view(LΨd_next, use_index_t),
                                  as_view(LΨsj, use_index_t));
        // WWᵀ
        guanaqo::blas::xsyrk_LT(real_t{1}, as_view(Wᵀj, use_index_t), //
                                real_t{0}, as_view(LΨdj, use_index_t));
        // -LΨs LΨsᵀ
        guanaqo::blas::xsyrk_LN(real_t{-1}, as_view(LΨsj, use_index_t), //
                                real_t{1}, as_view(LΨdj, use_index_t));
    }
    // chol(Θ)
    guanaqo::blas::xpotrf_L(as_view(factor.LΨd(0), use_index_t));
}

} // namespace hyhound::ocp
