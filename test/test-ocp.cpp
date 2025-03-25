#include <hyhound/ocp/riccati.hpp>
#include <hyhound/ocp/schur.hpp>

#include <gtest/gtest.h>
#include <limits>

TEST(OCP, factorization) {
    using std::exp2;
    using std::pow;
    using namespace hyhound;
    using namespace hyhound::ocp;
    std::mt19937 rng{321};
    std::normal_distribution<real_t> nrml{0, 5};
    std::bernoulli_distribution bern{0.8};
    const real_t ε = pow(std::numeric_limits<real_t>::epsilon(), real_t(0.5));

    // Generate random OCP
    OCPDataRiccati ocp{.N = 5, .nx = 13, .nu = 11, .ny = 8};
    ocp.init_random(123);

    auto rand_penalty = [&] { return bern(rng) ? exp2(nrml(rng)) : real_t{0}; };
    mat Σ = mat::Zero(ocp.ny, ocp.N + 1), Σ̃ = mat::Zero(ocp.ny, ocp.N + 1);
    std::ranges::generate(Σ.reshaped(), rand_penalty);
    std::ranges::generate(Σ̃.reshaped(), rand_penalty);

    // Solution using Riccati recursion
    RiccatiFactor factor_ric{.ocp = ocp};
    factor(factor_ric, Σ);
    solve(factor_ric);

    // Solution using Schur complement method
    auto ocp_sch = OCPDataSchur::from_riccati(ocp);
    SchurFactor factor_sch{.ocp = ocp_sch};
    factor(factor_sch, Σ);
    solve(factor_sch);

    // λ
    for (index_t j = 0; j < ocp.N + 1; ++j) {
        auto err = factor_sch.λ.col(j) - factor_ric.λ.col(j);
        EXPECT_TRUE(err.allFinite()) << " λ(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), 10 * ε) << " λ(" << j << ")";
    }
    // u
    for (index_t j = 0; j < ocp.N; ++j) {
        auto err = factor_sch.u(j) - factor_ric.u(j);
        EXPECT_TRUE(err.allFinite()) << " u(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), ε) << " u(" << j << ")";
    }
    // x
    for (index_t j = 0; j < ocp.N + 1; ++j) {
        auto err = factor_sch.x(j) - factor_ric.x(j);
        EXPECT_TRUE(err.allFinite()) << " x(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), ε) << " x(" << j << ")";
    }

    // Factorization update using Riccati recursion
    mat ΔΣ = Σ̃ - Σ;
    update(factor_ric, ΔΣ);
    RiccatiFactor new_factor_ric{.ocp = ocp};
    factor(new_factor_ric, Σ̃);

    // L
    for (index_t j = 0; j < ocp.N; ++j) {
        auto err = factor_ric.L(j) - new_factor_ric.L(j);
        EXPECT_TRUE(err.allFinite()) << " L(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), 10 * ε) << " L(" << j << ")";
    }
    {
        const auto j = ocp.N;
        auto err     = factor_ric.Lxx(j) - new_factor_ric.Lxx(j);
        EXPECT_TRUE(err.allFinite()) << " L(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), ε) << " L(" << j << ")";
    }

    // Factorization update using Schur complement method
    update(factor_sch, ΔΣ);
    SchurFactor new_factor_sch{.ocp = ocp_sch};
    factor(new_factor_sch, Σ̃);

    // LH
    for (index_t j = 0; j < ocp.N; ++j) {
        auto err = factor_sch.LH(j) - new_factor_sch.LH(j);
        EXPECT_TRUE(err.allFinite()) << " LH(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), 10 * ε) << " LH(" << j << ")";
    }
    {
        const auto j = ocp.N;
        auto err     = factor_sch.LHxx(j) - new_factor_sch.LHxx(j);
        EXPECT_TRUE(err.allFinite()) << " LH(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), 10 * ε) << " LH(" << j << ")";
    }

    // LΨd
    for (index_t j = 0; j < ocp.N + 1; ++j) {
        auto err = factor_sch.LΨd(j) - new_factor_sch.LΨd(j);
        EXPECT_TRUE(err.allFinite()) << " LΨd(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), ε) << " LΨd(" << j << ")";
    }
    // LΨs
    for (index_t j = 0; j < ocp.N; ++j) {
        auto err = factor_sch.LΨs(j) - new_factor_sch.LΨs(j);
        EXPECT_TRUE(err.allFinite()) << " LΨs(" << j << ")";
        EXPECT_LE(err.lpNorm<Eigen::Infinity>(), ε) << " LΨs(" << j << ")";
    }
}
