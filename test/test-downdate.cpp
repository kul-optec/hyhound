#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <random>

#include <hyhound/householder-updowndate.hpp>
#include <hyhound/linalg/blas-interface.hpp>
#include <guanaqo/eigen/view.hpp>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#if HYHOUND_WITH_OPENMP
#include <omp.h>
#endif

namespace hyhound {
namespace {

struct ProblemMatrices {
    Eigen::MatrixXd K̃, K, L, A;
};

ProblemMatrices generate_problem(index_t m, index_t n, index_t l = 0) {
#if HYHOUND_WITH_OPENMP
    int old_num_threads = omp_get_max_threads();
    omp_set_num_threads(std::thread::hardware_concurrency());
#endif
    if (l == 0)
        l = n;
    assert(l >= n);

    std::mt19937 rng{12345};
    std::uniform_real_distribution<> dist(0.0, 1.0);
    ProblemMatrices mat;
    mat.K̃.resize(l, l), mat.K.resize(l, l), mat.L.resize(l, n);
    mat.A.resize(l, m);
    std::ranges::generate(mat.K.reshaped(), [&] { return dist(rng); });
    std::ranges::generate(mat.A.reshaped(), [&] { return dist(rng); });
    const auto ldK = static_cast<index_t>(mat.K.outerStride()),
               ldA = static_cast<index_t>(mat.A.outerStride());
    linalg::xsyrk<real_t, index_t>(CblasColMajor, CblasLower, CblasTrans, n, n,
                                   1, mat.K.data(), ldK, 0, mat.K̃.data(), ldK);
    mat.K = mat.K̃;
    linalg::xsyrk<real_t, index_t>(CblasColMajor, CblasLower, CblasNoTrans, l,
                                   m, 1, mat.A.data(), ldA, 1, mat.K.data(),
                                   ldK);
    mat.L          = mat.K.leftCols(n);
    const auto ldL = static_cast<index_t>(mat.L.outerStride());
    index_t info   = 0;
    linalg::xpotrf<real_t, index_t>("L", n, mat.L.data(), ldL, &info);
    if (l > n) {
        linalg::xtrsm<real_t, index_t>(CblasColMajor, CblasRight, CblasLower,
                                       CblasTrans, CblasNonUnit, l - n, n, 1,
                                       mat.L.data(), ldL,
                                       mat.L.bottomRows(l - n).data(), ldL);
    }
    mat.L.triangularView<Eigen::StrictlyUpper>().setZero();
    mat.K̃.triangularView<Eigen::StrictlyUpper>() =
        mat.K̃.triangularView<Eigen::StrictlyLower>().transpose();
    mat.K.triangularView<Eigen::StrictlyUpper>() =
        mat.K.triangularView<Eigen::StrictlyLower>().transpose();

#if HYHOUND_WITH_OPENMP
    omp_set_num_threads(old_num_threads);
#endif

    return mat;
}

real_t calculate_error(const ProblemMatrices &matrices,
                       const Eigen::Ref<const Eigen::MatrixX<real_t>> &L̃) {
    Eigen::MatrixXd E = matrices.K̃;
    const auto n      = static_cast<index_t>(L̃.cols()),
               l      = static_cast<index_t>(L̃.rows()),
               ldL̃    = static_cast<index_t>(L̃.outerStride()),
               ldE    = static_cast<index_t>(E.outerStride());
#if HYHOUND_WITH_OPENMP
    int old_num_threads = omp_get_max_threads();
    omp_set_num_threads(std::thread::hardware_concurrency());
#endif
    linalg::xsyrk<real_t, index_t>(CblasColMajor, CblasLower, CblasNoTrans, n,
                                   n, -1, L̃.data(), ldL̃, 1, E.data(), ldE);
    if (l > n) {
        linalg::xtrsm<real_t, index_t>(
            CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
            l - n, n, 1, L̃.data(), ldL̃, E.bottomRows(l - n).data(), ldE);
        E.bottomLeftCorner(l - n, n) -= L̃.bottomRows(l - n);
    }
#if HYHOUND_WITH_OPENMP
    omp_set_num_threads(old_num_threads);
#endif
    E.triangularView<Eigen::StrictlyUpper>().setZero();
    return E.leftCols(n).lpNorm<Eigen::Infinity>();
}

} // namespace
} // namespace hyhound

using hyhound::index_t;
using hyhound::real_t;

const auto ε = 10 * std::pow(std::numeric_limits<real_t>::epsilon(), 0.5);

struct HyHDown : testing::TestWithParam<index_t> {};

TEST_P(HyHDown, VariousSizes) {
    index_t n = GetParam();
    for (index_t m : {1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 17, 31, 32}) {
        auto matrices     = hyhound::generate_problem(m, n);
        Eigen::MatrixXd L̃ = matrices.L;
        Eigen::MatrixXd Ã = matrices.A;
        hyhound::update_cholesky(as_view(L̃, guanaqo::with_index_type<index_t>),
                                 as_view(Ã, guanaqo::with_index_type<index_t>),
                                 hyhound::Downdate());
        real_t residual = hyhound::calculate_error(matrices, L̃);
        EXPECT_LE(residual, ε) << "m=" << m;
    }
}

struct HyHDownRect : testing::TestWithParam<index_t> {};

TEST_P(HyHDownRect, VariousSizes) {
    index_t n       = GetParam();
    const index_t m = 13;
    for (index_t l = n; l < n + 7; ++l) {
        auto matrices     = hyhound::generate_problem(m, n, l);
        Eigen::MatrixXd L̃ = matrices.L;
        Eigen::MatrixXd Ã = matrices.A;
        hyhound::update_cholesky(as_view(L̃, guanaqo::with_index_type<index_t>),
                                 as_view(Ã, guanaqo::with_index_type<index_t>),
                                 hyhound::Downdate());
        real_t residual = hyhound::calculate_error(matrices, L̃);
        EXPECT_LE(residual, ε) << "l=" << l;
    }
}

INSTANTIATE_TEST_SUITE_P(Cholundate, HyHDown, testing::Range<index_t>(1, 256));
INSTANTIATE_TEST_SUITE_P(Cholundate, HyHDownRect,
                         testing::Range<index_t>(1, 256));
