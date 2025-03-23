#include "ocp/riccati.hpp"

#include <hyhound/householder-updowndate.hpp>
#include <hyhound/updown.hpp>
#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/span.hpp>
#include <guanaqo/eigen/view.hpp>

namespace hyhound::ocp {

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
            auto YNJ  = YN->leftCols(nJ).bottomRows(ocp.nx);
            auto ΦNJ  = ΦN->leftCols(nJ);
            auto SNJ  = factor.S.leftCols(nJ);
            auto FNm1 = ocp.F(ocp.N - 1);
            guanaqo::blas::xgemm<real_t, index_t>(
                CblasColMajor, CblasTrans, CblasNoTrans, ΦNJ.rows(), ΦNJ.cols(),
                FNm1.rows(), real_t{1}, FNm1.data(), FNm1.outerStride(),
                YNJ.data(), YNJ.outerStride(), real_t{0}, ΦNJ.data(),
                ΦNJ.outerStride());
            update_cholesky(
                guanaqo::as_view(LxxN, guanaqo::with_index_type<index_t>),
                guanaqo::as_view(YNJ, guanaqo::with_index_type<index_t>),
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
        update_cholesky(
            guanaqo::as_view(Luuxuj, guanaqo::with_index_type<index_t>),
            guanaqo::as_view(YjJ, guanaqo::with_index_type<index_t>),
            hyhound::UpDowndate{guanaqo::as_span(SjJ.reshaped())});
        if (j > 0) {
            auto FNm1 = ocp.F(j - 1);
            guanaqo::blas::xgemm<real_t, index_t>(
                CblasColMajor, CblasTrans, CblasNoTrans, ΦjJ.rows(), ΦjJ.cols(),
                FNm1.rows(), real_t{1}, FNm1.data(), FNm1.outerStride(),
                YjJx.data(), YjJx.outerStride(), real_t{0}, ΦjJ.data(),
                ΦjJ.outerStride());
        }
        update_cholesky(
            guanaqo::as_view(Lxxj, guanaqo::with_index_type<index_t>),
            guanaqo::as_view(YjJx, guanaqo::with_index_type<index_t>),
            hyhound::UpDowndate{guanaqo::as_span(SjJ.reshaped())});
    }
}

} // namespace hyhound::ocp
