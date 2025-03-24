#include "ocp/riccati.hpp"

#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

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
        guanaqo::blas::xgemmt_LTN(real_t{1}, as_view(CN, use_index_t),
                                  as_view(ΣCN, use_index_t), real_t{1},
                                  as_view(LxxN, use_index_t));
        // Lxx(N) = chol(H(N))
        guanaqo::blas::xpotrf_L(as_view(LxxN, use_index_t));
    }
    for (index_t j = ocp.N; j-- > 0;) {
        auto Gj       = ocp.G(j);
        auto &ΣGj     = factor.ΣG;
        auto Lj       = factor.L(j);
        auto Lxx_next = factor.Lxx(j + 1);
        auto &Vᵀ      = factor.Vᵀ;
        // V(j) = F(j)ᵀ Lxx(j+1), Vᵀ = Lxx(j+1)ᵀ F(j)
        Vᵀ = ocp.F(j);
        guanaqo::blas::xtrmm<real_t, index_t>( // TODO
            CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit,
            Vᵀ.rows(), Vᵀ.cols(), real_t{1}, Lxx_next.data(),
            Lxx_next.outerStride(), Vᵀ.data(), Vᵀ.outerStride());
        // H(j) = Hl(j) + G(j)ᵀΣ(j) G(j) + V(j)Vᵀ(j)
        ΣGj = Σ.col(j).asDiagonal() * Gj;
        Lj.triangularView<Eigen::Lower>() =
            ocp.H(j).triangularView<Eigen::Lower>();
        guanaqo::blas::xgemmt_LTN(real_t{1}, as_view(Gj, use_index_t),
                                  as_view(ΣGj, use_index_t), real_t{1},
                                  as_view(Lj, use_index_t));
        guanaqo::blas::xsyrk_LT(real_t{1}, as_view(Vᵀ, use_index_t), real_t{1},
                                as_view(Lj, use_index_t));
        // L(j) = chol(H(j))
        guanaqo::blas::xpotrf_L(as_view(Lj, use_index_t));
    }
}

} // namespace hyhound::ocp
