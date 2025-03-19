#include "ocp/schur.hpp"

#include <hyhound/linalg/blas-interface.hpp>

namespace hyhound::ocp {

static void update_schur_rank_one(SchurFactor &factor,
                                  index_t j /* stage index */,
                                  index_t i /* constraint index */,
                                  real_t Σji /* penalty factor change */) {
    const auto &ocp = factor.ocp;
    using std::abs;
    using std::copysign;
    using std::sqrt;
    auto Hj = j == ocp.N ? Eigen::Ref<mat>(factor.LHxx(j))
                         : Eigen::Ref<mat>(factor.LH(j));
    if (j == ocp.N)
        factor.e̅.topRows(ocp.nx) = ocp.C(j).row(i).transpose() * sqrt(abs(Σji));
    else
        factor.e̅ = ocp.G(j).row(i).transpose() * sqrt(abs(Σji));
    factor.ẽ = factor.e̅;
    // e̅ = H ẽ
    linalg::xtrsv<real_t, index_t>(
        CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, Hj.rows(),
        Hj.data(), Hj.outerStride(), factor.ẽ.data(), index_t{1});
    linalg::xtrsv<real_t, index_t>(
        CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, Hj.rows(),
        Hj.data(), Hj.outerStride(), factor.ẽ.data(), index_t{1});
    factor.ẽ *= (1 / sqrt(1 + copysign(factor.ẽ.dot(factor.e̅), Σji)));
    if (j == ocp.N) {
        factor.ψ.topRows(ocp.nx) = -factor.ẽ.topRows(ocp.nx);
        factor.ψ.bottomRows(ocp.nx).setZero();
    } else {
        auto Fj = ocp.F(j);
        linalg::xgemv<real_t, index_t>(
            CblasColMajor, CblasNoTrans, Fj.rows(), Fj.cols(), real_t{1},
            Fj.data(), Fj.outerStride(), factor.ẽ.data(), index_t{1}, real_t{0},
            factor.ψ.data(), index_t{1});
        factor.ψ.bottomRows(ocp.nx) = -factor.ẽ.topRows(ocp.nx);
    }
    real_t α                    = -copysign(real_t{1}, Σji);
    index_t last_affected_stage = j == ocp.N ? j : j + 1;
    for (index_t jj = last_affected_stage + 1; jj-- > 0;) {
        // Gill et al. Algorithm C1 for updating LΨd and LΨs
        for (index_t r = 0; r < ocp.nx; ++r) {
            real_t p     = factor.ψ(r);
            real_t λ     = factor.LΨd(jj)(r, r);
            real_t d     = λ * λ;
            real_t d̃     = d + α * (p * p);
            real_t λ̃     = copysign(sqrt(d̃), λ);
            real_t inv_d̃ = 1 / d̃;
            real_t b     = p * α * inv_d̃;

            auto ad = factor.ψ.middleRows(r, ocp.nx - r);
            auto ld = factor.LΨd(jj).col(r).bottomRows(ocp.nx - r);
            ad -= (p / λ) * ld;
            ld    = (λ̃ / λ) * ld + (b * λ̃) * ad;
            ld(0) = λ̃;
            if (jj > 0) {
                auto as = factor.ψ.bottomRows(ocp.nx);
                auto ls = factor.LΨs(jj - 1).col(r);
                as -= (p / λ) * ls;
                ls = (λ̃ / λ) * ls + (b * λ̃) * as;
            }
            α *= d * inv_d̃;
        }
        if (jj > 0) {
            factor.ψ.topRows(ocp.nx) = factor.ψ.bottomRows(ocp.nx);
            factor.ψ.bottomRows(ocp.nx).setZero();
        }
    }
    // Gill et al. Algorithm C1 for updating LH
    α = copysign(real_t{1}, Σji);
    for (index_t r = 0; r < Hj.cols(); ++r) {
        real_t p     = factor.e̅(r);
        real_t λ     = Hj(r, r);
        real_t d     = λ * λ;
        real_t d̃     = d + α * (p * p);
        real_t λ̃     = copysign(sqrt(d̃), λ);
        real_t inv_d̃ = 1 / d̃;
        real_t b     = p * α * inv_d̃;

        auto ad = factor.e̅.middleRows(r, Hj.rows() - r);
        auto ld = Hj.col(r).middleRows(r, Hj.rows() - r);
        ad -= (p / λ) * ld;
        ld    = (λ̃ / λ) * ld + (b * λ̃) * ad;
        ld(0) = λ̃;
        α *= d * inv_d̃;
    }
}

void update(SchurFactor &factor, Eigen::Ref<const mat> ΔΣ) {
    const auto &ocp = factor.ocp;
    for (index_t j = ocp.N + 1; j-- > 0;) {
        for (index_t i = 0; i < ocp.ny; ++i) {
            real_t Σji = ΔΣ.col(j)(i);
            if (Σji != 0)
                update_schur_rank_one(factor, j, i, Σji);
        }
    }
}

} // namespace hyhound::ocp
