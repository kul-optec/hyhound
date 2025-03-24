#include "ocp/riccati.hpp"

#include <hyhound/householder-updowndate.hpp>
#include <hyhound/updown.hpp>
#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/span.hpp>
#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

void update(RiccatiFactor &factor, Eigen::Ref<const mat> ΔΣ) {
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
            guanaqo::blas::xgemm_TN(real_t{1},
                                    as_view(ocp.F(ocp.N - 1), use_index_t),
                                    as_view(YNJ, use_index_t), real_t{0},
                                    as_view(ΦNJ, use_index_t));
            update_cholesky(
                guanaqo::as_view(LxxN, use_index_t),
                guanaqo::as_view(YNJ, use_index_t),
                hyhound::UpDowndate{guanaqo::as_span(SNJ.reshaped())});
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
        update_cholesky(guanaqo::as_view(Luuxuj, use_index_t),
                        guanaqo::as_view(YjJ, use_index_t),
                        hyhound::UpDowndate{guanaqo::as_span(SjJ.reshaped())});
        if (j > 0)
            guanaqo::blas::xgemm_TN(real_t{1},
                                    as_view(ocp.F(j - 1), use_index_t),
                                    as_view(YjJx, use_index_t), real_t{0},
                                    as_view(ΦjJ, use_index_t));
        update_cholesky(guanaqo::as_view(Lxxj, use_index_t),
                        guanaqo::as_view(YjJx, use_index_t),
                        hyhound::UpDowndate{guanaqo::as_span(SjJ.reshaped())});
    }
}

} // namespace hyhound::ocp
