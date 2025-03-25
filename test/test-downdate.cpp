#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <random>

#include <hyhound/householder-updowndate.hpp>
#include <guanaqo/blas/hl-blas-interface.hpp>
#include <guanaqo/eigen/view.hpp>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#if GUANAQO_WITH_OPENMP
#include <omp.h>
#endif

namespace hyhound {
namespace {

struct ProblemMatrices {
    Eigen::MatrixX<real_t> K̃, K, L, A;
};

constexpr auto use_index_t = guanaqo::with_index_type<index_t>;

ProblemMatrices generate_problem(index_t m, index_t n, index_t l = 0) {
#if GUANAQO_WITH_OPENMP
    int old_num_threads = omp_get_max_threads();
    omp_set_num_threads(std::thread::hardware_concurrency());
#endif
    if (l == 0)
        l = n;
    assert(l >= n);

    std::mt19937 rng{12345};
    std::uniform_real_distribution<real_t> dist(-1, 1);
    ProblemMatrices mat;
    mat.K̃.resize(l, l), mat.K.resize(l, l), mat.L.resize(l, n);
    mat.A.resize(l, m);
    std::ranges::generate(mat.K.reshaped(), [&] { return dist(rng); });
    std::ranges::generate(mat.A.reshaped(), [&] { return dist(rng); });
    guanaqo::blas::xsyrk_LT(real_t{1}, as_view(mat.K, use_index_t), //
                            real_t{0}, as_view(mat.K̃, use_index_t));
    mat.K = mat.K̃;
    guanaqo::blas::xsyrk_LN(real_t{1}, as_view(mat.A, use_index_t), //
                            real_t{1}, as_view(mat.K, use_index_t));
    mat.L = mat.K.leftCols(n);
    guanaqo::blas::xpotrf_L(as_view(mat.L.topRows(n), use_index_t));
    if (l > n) {
        guanaqo::blas::xtrsm_RLTN(
            real_t{1}, as_view(mat.L.topRows(n), use_index_t),
            as_view(mat.L.bottomRows(l - n), use_index_t));
    }
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
    const auto n             = static_cast<index_t>(L̃.cols()),
               l             = static_cast<index_t>(L̃.rows());
    Eigen::MatrixX<real_t> E = matrices.K̃;
#if GUANAQO_WITH_OPENMP
    int old_num_threads = omp_get_max_threads();
    omp_set_num_threads(std::thread::hardware_concurrency());
#endif
    guanaqo::blas::xsyrk_LN(                            //
        real_t{-1}, as_view(L̃.topRows(n), use_index_t), //
        real_t{1}, as_view(E.topLeftCorner(n, n), use_index_t));
    if (l > n) {
        guanaqo::blas::xtrsm_RLTN(
            real_t{1}, as_view(L̃.topRows(n), use_index_t),
            as_view(E.bottomLeftCorner(l - n, n), use_index_t));
        E.bottomLeftCorner(l - n, n) -= L̃.bottomRows(l - n);
    }
#if GUANAQO_WITH_OPENMP
    omp_set_num_threads(old_num_threads);
#endif
    E.triangularView<Eigen::StrictlyUpper>().setZero();
    return E.leftCols(n).lpNorm<Eigen::Infinity>();
}

} // namespace
} // namespace hyhound

using hyhound::index_t;
using hyhound::real_t;
using hyhound::use_index_t;

const auto ε = 10 * std::pow(std::numeric_limits<real_t>::epsilon(), 0.5);

struct HyHDown : testing::TestWithParam<index_t> {};

TEST_P(HyHDown, VariousSizes) {
    index_t n = GetParam();
    for (index_t m : {1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 17, 31, 32}) {
        auto matrices            = hyhound::generate_problem(m, n);
        Eigen::MatrixX<real_t> L̃ = matrices.L;
        Eigen::MatrixX<real_t> Ã = matrices.A;
        hyhound::update_cholesky(as_view(L̃, use_index_t),
                                 as_view(Ã, use_index_t), hyhound::Downdate());
        real_t residual = hyhound::calculate_error(matrices, L̃);
        EXPECT_LE(residual, ε) << "m=" << m;
    }
}

struct HyHDownRect : testing::TestWithParam<index_t> {};

TEST_P(HyHDownRect, VariousSizes) {
    index_t n       = GetParam();
    const index_t m = 13;
    for (index_t l = n; l < n + 7; ++l) {
        auto matrices            = hyhound::generate_problem(m, n, l);
        Eigen::MatrixX<real_t> L̃ = matrices.L;
        Eigen::MatrixX<real_t> Ã = matrices.A;
        hyhound::update_cholesky(as_view(L̃, use_index_t),
                                 as_view(Ã, use_index_t), hyhound::Downdate());
        real_t residual = hyhound::calculate_error(matrices, L̃);
        EXPECT_LE(residual, ε) << "l=" << l;
    }
}

INSTANTIATE_TEST_SUITE_P(Cholundate, HyHDown, testing::Range<index_t>(1, 256));
INSTANTIATE_TEST_SUITE_P(Cholundate, HyHDownRect,
                         testing::Range<index_t>(1, 256));
