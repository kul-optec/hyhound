#include <hyhound/ocp/riccati.hpp>

#include <hyhound/householder-updowndate.hpp>
#include <hyhound/updown.hpp>
#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/span.hpp>

namespace hyhound::ocp {

static auto with_signature_matrix(auto &&S) {
    return hyhound::UpDowndate{guanaqo::as_span(S.reshaped())};
}

void update(RiccatiFactor &factor, Eigen::Ref<const mat> ΔΣ) {
    using namespace guanaqo::blas;
    using std::abs;
    using std::copysign;
    using std::sqrt;
    const auto &ocp = factor.ocp;
    index_t nJ      = 0;
    {
        auto CN   = ocp.C(ocp.N);
        auto LxxN = factor.Lxx(ocp.N);
        auto *ΦN  = (ocp.N & 1) == 0 ? &factor.ΦY : &factor.YΦ;
        auto *YN  = (ocp.N & 1) == 1 ? &factor.ΦY : &factor.YΦ;
        for (index_t i = 0; i < ocp.ny; ++i) {
            auto Σji = ΔΣ.col(ocp.N)(i);
            if (Σji == 0)
                continue;
            YN->col(nJ).bottomRows(ocp.nx) =
                CN.row(i).transpose() * sqrt(abs(Σji));
            factor.S(0, nJ) = copysign(real_t{0}, Σji);
            ++nJ;
        }
        if (nJ > 0) {
            auto YNJ = YN->leftCols(nJ).bottomRows(ocp.nx);
            auto ΦNJ = ΦN->leftCols(nJ);
            auto SNJ = factor.S.leftCols(nJ);
            xgemm_TN(real_t{1}, vw(ocp.F(ocp.N - 1)), vw(YNJ), real_t{0},
                     vw(ΦNJ));
            update_cholesky(vw(LxxN), vw(YNJ), with_signature_matrix(SNJ));
        }
    }
    for (index_t j = ocp.N; j-- > 0;) {
        auto *ΦN = (j & 1) == 0 ? &factor.ΦY : &factor.YΦ;
        auto *YN = (j & 1) == 1 ? &factor.ΦY : &factor.YΦ;
        for (index_t i = 0; i < ocp.ny; ++i) {
            auto Σji = ΔΣ.col(j)(i);
            if (Σji == 0)
                continue;
            YN->col(nJ)     = ocp.G(j).row(i).transpose() * sqrt(abs(Σji));
            factor.S(0, nJ) = copysign(real_t{0}, Σji);
            ++nJ;
        }
        if (nJ == 0)
            continue;
        auto YjJ    = YN->leftCols(nJ);
        auto YjJx   = YjJ.bottomRows(ocp.nx);
        auto ΦjJ    = ΦN->leftCols(nJ);
        auto SjJ    = factor.S.leftCols(nJ);
        auto Lj     = factor.L(j);
        auto Luuxuj = Lj.leftCols(ocp.nu);
        auto Lxxj   = Lj.bottomRightCorner(ocp.nx, ocp.nx);
        update_cholesky(vw(Luuxuj), vw(YjJ), with_signature_matrix(SjJ));
        if (j > 0)
            xgemm_TN(real_t{1}, vw(ocp.F(j - 1)), vw(YjJx), real_t{0}, vw(ΦjJ));
        update_cholesky(vw(Lxxj), vw(YjJx), with_signature_matrix(SjJ));
    }
}

} // namespace hyhound::ocp
