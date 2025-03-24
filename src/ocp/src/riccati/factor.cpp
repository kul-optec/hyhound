#include <hyhound/ocp/riccati.hpp>

#include <guanaqo/blas/hl-blas-interface.hpp>

namespace hyhound::ocp {

void factor(RiccatiFactor &factor, Eigen::Ref<const mat> Σ) {
    using namespace guanaqo::blas;
    const auto &ocp = factor.ocp;
    /* j = N */ {
        auto CN   = ocp.C(ocp.N);
        auto ΣCN  = factor.ΣG.rightCols(ocp.nx);
        auto LxxN = factor.Lxx(ocp.N);
        // H(N) = Q(N) + C(N)ᵀΣ(N) G(N)
        ΣCN = Σ.col(ocp.N).asDiagonal() * CN;
        LxxN.triangularView<Eigen::Lower>() =
            ocp.Q(ocp.N).triangularView<Eigen::Lower>();
        xgemmt_LTN(real_t{1}, vw(CN), vw(ΣCN), real_t{1}, vw(LxxN));
        // Lxx(N) = chol(H(N))
        xpotrf_L(vw(LxxN));
    }
    for (index_t j = ocp.N; j-- > 0;) {
        auto Gj       = ocp.G(j);
        auto &ΣGj     = factor.ΣG;
        auto Lj       = factor.L(j);
        auto Lxx_next = factor.Lxx(j + 1);
        auto &Vᵀ      = factor.Vᵀ;
        // V(j) = F(j)ᵀ Lxx(j+1), Vᵀ = Lxx(j+1)ᵀ F(j)
        Vᵀ = ocp.F(j);
        xtrmm_LLTN(real_t{1}, vw(Lxx_next), vw(Vᵀ));
        // H(j) = Hl(j) + G(j)ᵀΣ(j) G(j) + V(j)Vᵀ(j)
        ΣGj = Σ.col(j).asDiagonal() * Gj;
        Lj.triangularView<Eigen::Lower>() =
            ocp.H(j).triangularView<Eigen::Lower>();
        xgemmt_LTN(real_t{1}, vw(Gj), vw(ΣGj), real_t{1}, vw(Lj));
        xsyrk_LT(real_t{1}, vw(Vᵀ), real_t{1}, vw(Lj));
        // L(j) = chol(H(j))
        xpotrf_L(vw(Lj));
    }
}

} // namespace hyhound::ocp
