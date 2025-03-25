#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <random>

#include <hyhound/householder-updowndate.hpp>
#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/span.hpp>
#include <guanaqo/eigen/view.hpp>
using guanaqo::as_span;

#include <Eigen/Cholesky>
#include <Eigen/Core>

#if GUANAQO_WITH_OPENMP
#include <omp.h>
#endif

namespace hyhound {
namespace {

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

struct ProblemMatrices {
    Eigen::MatrixX<real_t> K̃, K, L, A;
    Eigen::VectorX<real_t> S;
};

ProblemMatrices generate_problem(index_t m, index_t n) {
#if GUANAQO_WITH_OPENMP
    int old_num_threads = omp_get_max_threads();
    omp_set_num_threads(std::thread::hardware_concurrency());
#endif

    std::mt19937 rng{12345};
    std::uniform_real_distribution<real_t> dist(-1, 1);
    std::bernoulli_distribution brnl(0.5);
    ProblemMatrices mat;
    mat.K̃.resize(n, n), mat.K.resize(n, n), mat.L.resize(n, n);
    mat.A.resize(n, m), mat.S.resize(m);
    std::ranges::generate(mat.K.reshaped(), [&] { return dist(rng); });
    std::ranges::generate(mat.A.reshaped(), [&] { return dist(rng); });
    Eigen::MatrixX<real_t> Ad(n, m), Au(n, m);
    Eigen::Index ju = 0, jd = 0;
    for (index_t i = 0; i < m; ++i) {
        if (brnl(rng)) {
            Au.col(ju++) = mat.A.col(i);
            mat.S(i)     = +real_t(0);
        } else {
            Ad.col(jd++) = mat.A.col(i);
            mat.S(i)     = -real_t(0);
        }
    }
    Au = Au.leftCols(ju).eval();
    Ad = Ad.leftCols(jd).eval();
    if constexpr (std::is_same_v<real_t, float>)
        mat.K += 10 * Eigen::MatrixX<real_t>::Identity(n, n);
    // K̃ ← KᵀK
    guanaqo::blas::xsyrk_LT(real_t{1}, as_view(mat.K, use_index_t), //
                            real_t{0}, as_view(mat.K̃, use_index_t));
    // K ← K̃
    mat.K = mat.K̃;
    // K += A₋A₋ᵀ
    guanaqo::blas::xsyrk_LN(real_t{1}, as_view(Ad, use_index_t), //
                            real_t{1}, as_view(mat.K, use_index_t));
    // K̃ += A₊A₊ᵀ
    guanaqo::blas::xsyrk_LN(real_t{1}, as_view(Au, use_index_t), //
                            real_t{1}, as_view(mat.K̃, use_index_t));
    // L = chol(K)
    mat.L = mat.K;
    guanaqo::blas::xpotrf_L(as_view(mat.L, use_index_t));
    mat.L.triangularView<Eigen::StrictlyUpper>().setZero();
    mat.K̃.triangularView<Eigen::StrictlyUpper>() =
        mat.K̃.triangularView<Eigen::StrictlyLower>().transpose();
    mat.K.triangularView<Eigen::StrictlyUpper>() =
        mat.K.triangularView<Eigen::StrictlyLower>().transpose();

#if GUANAQO_WITH_OPENMP
    omp_set_num_threads(old_num_threads);
#endif

    return mat;
}

real_t calculate_error(const ProblemMatrices &matrices,
                       const Eigen::Ref<const Eigen::MatrixX<real_t>> &L̃) {
    Eigen::MatrixX<real_t> E = matrices.K̃;
#if GUANAQO_WITH_OPENMP
    int old_num_threads = omp_get_max_threads();
    omp_set_num_threads(std::thread::hardware_concurrency());
#endif
    guanaqo::blas::xsyrk_LN(real_t{-1}, as_view(L̃, use_index_t), //
                            real_t{1}, as_view(E, use_index_t));
#if GUANAQO_WITH_OPENMP
    omp_set_num_threads(old_num_threads);
#endif
    E.triangularView<Eigen::StrictlyUpper>().setZero();
    return E.lpNorm<Eigen::Infinity>();
}

} // namespace
} // namespace hyhound

using hyhound::index_t;
using hyhound::real_t;
using hyhound::use_index_t;

const auto ε = std::pow(std::numeric_limits<real_t>::epsilon(), 0.5);

struct HyHUpDown : testing::TestWithParam<index_t> {};

TEST_P(HyHUpDown, VariousSizes) {
    index_t n = GetParam();
    for (index_t m : {1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 17, 31, 32}) {
        auto matrices            = hyhound::generate_problem(m, n);
        Eigen::MatrixX<real_t> L̃ = matrices.L;
        Eigen::MatrixX<real_t> Ã = matrices.A;
        Eigen::VectorX<real_t> S̃ = matrices.S;
        update_cholesky(as_view(L̃, use_index_t), as_view(Ã, use_index_t),
                        hyhound::UpDowndate{as_span(S̃)});
        real_t residual = hyhound::calculate_error(matrices, L̃);
        EXPECT_LE(residual, ε) << "m=" << m;
    }
}

INSTANTIATE_TEST_SUITE_P(Cholundate, HyHUpDown,
                         testing::Range<index_t>(1, 256));
