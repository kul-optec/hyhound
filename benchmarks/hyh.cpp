#include <hyhound/householder-updowndate.hpp>
#include <hyhound/linalg/blas-interface.hpp>
#include <hyhound-version.h>

#include <guanaqo/eigen/view.hpp>
#include <guanaqo/preprocessor.h>

#include <Eigen/Dense>
#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdlib>
#include <format>
#include <map>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

using hyhound::index_t;
using hyhound::real_t;
static constexpr auto use_index_t = guanaqo::with_index_type<index_t>;
using std::pow;

#if HYHOUND_WITH_OPENMP
#include <omp.h>
#include <thread>
#endif

struct ProblemMatrices {
    Eigen::MatrixXd K̃, K, L, A;
};
using cache_t = std::map<std::pair<index_t, index_t>, ProblemMatrices>;
std::mutex cache_mtx;
cache_t cache;

struct CholeskyFixture : benchmark::Fixture {
    index_t m, n;
    Eigen::MatrixXd L̃;
    cache_t::const_iterator matrices;

    static cache_t::const_iterator generate_problem(index_t m, index_t n) {
        std::lock_guard lck{cache_mtx};
        auto [it, inserted] =
            cache.try_emplace(std::pair{m, n}, ProblemMatrices{});
        auto &mat = it->second;
        if (!inserted)
            return it;

#if HYHOUND_WITH_OPENMP
        int old_num_threads = omp_get_max_threads();
        omp_set_num_threads(std::thread::hardware_concurrency() / 2);
#endif

        std::mt19937 rng{12345};
        std::uniform_real_distribution<> dist(0.0, 1.0);
        mat.K̃.resize(n, n), mat.K.resize(n, n), mat.L.resize(n, n);
        mat.A.resize(n, m);
        std::ranges::generate(mat.K.reshaped(), [&] { return dist(rng); });
        std::ranges::generate(mat.A.reshaped(), [&] { return dist(rng); });
        const auto ldK = static_cast<index_t>(mat.K.outerStride()),
                   ldA = static_cast<index_t>(mat.A.outerStride());
        hyhound::linalg::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasTrans, n, n, 1, mat.K.data(), ldK,
            0, mat.K̃.data(), ldK);
        mat.K = mat.K̃;
        hyhound::linalg::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, n, m, 1, mat.A.data(), ldA,
            1, mat.K.data(), ldK);
        mat.L          = mat.K;
        const auto ldL = static_cast<index_t>(mat.L.outerStride());
        index_t info   = 0;
        hyhound::linalg::xpotrf<real_t, index_t>("L", n, mat.L.data(), ldL,
                                                 &info);
        mat.L.triangularView<Eigen::StrictlyUpper>().setZero();
        mat.K̃.triangularView<Eigen::StrictlyUpper>() =
            mat.K̃.triangularView<Eigen::StrictlyLower>().transpose();
        mat.K.triangularView<Eigen::StrictlyUpper>() =
            mat.K.triangularView<Eigen::StrictlyLower>().transpose();

#if HYHOUND_WITH_OPENMP
        omp_set_num_threads(old_num_threads);
#endif
        using namespace std::chrono_literals;
        // Prevent down-clocking after heavy multi-threaded code
        std::this_thread::sleep_for(1000ms);

        return it;
    }

    void SetUp(benchmark::State &state) final {
        m        = static_cast<index_t>(state.range(0));
        n        = static_cast<index_t>(state.range(1));
        matrices = generate_problem(m, n);
        L̃.resize(n, n);
    }

    void TearDown(benchmark::State &state) final {
        Eigen::MatrixXd E = matrices->second.K̃;
        const auto n      = static_cast<index_t>(L̃.rows()),
                   ldL̃    = static_cast<index_t>(L̃.outerStride()),
                   ldE    = static_cast<index_t>(E.outerStride());
#if HYHOUND_WITH_OPENMP
        int old_num_threads = omp_get_max_threads();
        omp_set_num_threads(std::thread::hardware_concurrency() / 2);
#endif
        hyhound::linalg::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, n, n, -1, L̃.data(), ldL̃, 1,
            E.data(), ldE);
#if HYHOUND_WITH_OPENMP
        omp_set_num_threads(old_num_threads);
#endif
        E.triangularView<Eigen::StrictlyUpper>().setZero();
        real_t r          = E.lpNorm<Eigen::Infinity>();
        const char *color = r < 1e-9 ? "" : "\x1b[0;31m";
        const char *reset = r < 1e-9 ? "" : "\x1b[0m";
        state.SetLabel(std::format("{}resid={:7e}{}", color, r, reset));
        compute_flops(state);
    }

    virtual void compute_flops(benchmark::State &state) const {
        auto nd = static_cast<double>(n), md = static_cast<double>(m);
        auto flop_cnt                 = pow(nd, 3) / 6 + md * pow(nd, 2) / 2;
        state.counters["GFLOP count"] = {1e-9 * flop_cnt};
        state.counters["GFLOPS"]      = {
            1e-9 * flop_cnt, benchmark::Counter::kIsIterationInvariantRate};
    }

    template <auto Func>
    void runUpdateBenchmark(benchmark::State &state) {
        Eigen::MatrixXd Ã(m, n);
        for (auto _ : state) {
            state.PauseTiming();
            Ã = matrices->second.A;
            L̃ = matrices->second.L;
            state.ResumeTiming();
            benchmark::DoNotOptimize(Ã.data());
            benchmark::DoNotOptimize(L̃.data());
            Func(as_view(L̃, use_index_t), as_view(Ã, use_index_t), {});
            benchmark::ClobberMemory();
        }
    }
};

template <index_t BsL, index_t BsA, auto...>
struct BlockedFixture : CholeskyFixture {
    void compute_flops(benchmark::State &state) const override {
        auto nd = static_cast<double>(n), md = static_cast<double>(m);
        auto flop_cnt                 = md * pow(nd, 2) + pow(nd, 2) + md * nd;
        state.counters["GFLOP count"] = {1e-9 * flop_cnt};
        state.counters["GFLOPS"]      = {
            1e-9 * flop_cnt, benchmark::Counter::kIsIterationInvariantRate};
    }
};

BENCHMARK_DEFINE_F(CholeskyFixture, full_factorize)(benchmark::State &state) {
    for (auto _ : state) {
        state.PauseTiming();
        L̃ = matrices->second.K.triangularView<Eigen::Lower>();
        state.ResumeTiming();
        benchmark::DoNotOptimize(matrices->second.A.data());
        benchmark::DoNotOptimize(L̃.data());
        const auto ldA = static_cast<index_t>(matrices->second.A.outerStride());
        const auto ldL = static_cast<index_t>(L̃.outerStride());
        hyhound::linalg::xsyrk<real_t, index_t>(
            CblasColMajor, CblasLower, CblasNoTrans, n, m, -1,
            matrices->second.A.data(), ldA, 1, L̃.data(), ldL);
        index_t info = 0;
        hyhound::linalg::xpotrf<real_t, index_t>("L", n, L̃.data(), ldL, &info);
        benchmark::ClobberMemory();
    }
}

std::vector<::benchmark::internal::Benchmark *> benchmarks;

#define BM_BLK_REGISTER_F(BaseClass, Method)                                   \
    BM_BLK_PRIVATE_REGISTER_F(BENCHMARK_PRIVATE_CONCAT_NAME(BaseClass, Method))
#define BM_BLK_PRIVATE_REGISTER_F(TestName)                                    \
    BENCHMARK_PRIVATE_DECLARE(TestName) = [] {                                 \
        using ::benchmark::internal::RegisterBenchmarkInternal;                \
        auto *bm = RegisterBenchmarkInternal(new TestName());                  \
        benchmarks.push_back(bm);                                              \
        return bm;                                                             \
    }()
#define BM_BLK_IMPL_NAME(name, ...)                                            \
    GUANAQO_CONCATENATE_TOKENS(name, GUANAQO_JOIN_TOKENS(__VA_ARGS__))
#define BM_BLK_NAME(name, ...)                                                 \
    #name "<" GUANAQO_JOIN_STRINGS(", ", __VA_ARGS__) ">"
#define BENCHMARK_BLOCKED(name, func, updown, ...)                             \
    BENCHMARK_TEMPLATE_DEFINE_F(                                               \
        BlockedFixture, BM_BLK_IMPL_NAME(name, __VA_ARGS__), __VA_ARGS__)      \
    (benchmark::State & state) {                                               \
        this->runUpdateBenchmark<func<{__VA_ARGS__}, updown>>(state);          \
    }                                                                          \
    BM_BLK_REGISTER_F(BlockedFixture, BM_BLK_IMPL_NAME(name, __VA_ARGS__))     \
        ->Name(BM_BLK_NAME(name, __VA_ARGS__))

void configure_benchmarks(bool fix_m, int64_t M, int64_t N);

int main(int argc, char **argv) {
    if (argc < 1)
        return 1;
    // Parse command-line arguments
    std::vector<char *> argvv{argv, argv + argc};
    bool fix_m = false;
    int64_t M = 8, N = 256;
    for (auto it = argvv.begin(); it != argvv.end();) {
        std::string_view arg = *it;
        if (arg == "--fix-m") {
            fix_m = true;
            argvv.erase(it);
        } else if (arg == "--fix-n") {
            fix_m = false;
            argvv.erase(it);
        } else if (std::string_view flag = "--m="; arg.starts_with(flag)) {
            M = std::stoi(std::string(arg.substr(flag.size())));
            argvv.erase(it);
        } else if (std::string_view flag = "--n="; arg.starts_with(flag)) {
            N = std::stoi(std::string(arg.substr(flag.size())));
            argvv.erase(it);
        } else {
            ++it;
        }
    }
    configure_benchmarks(fix_m, M, N);

    argc = static_cast<int>(argvv.size());
    ::benchmark::Initialize(&argc, argvv.data());
    if (::benchmark::ReportUnrecognizedArguments(argc, argvv.data()))
        return 1;
#if HYHOUND_WITH_OPENMP
    benchmark::AddCustomContext("OMP_NUM_THREADS",
                                std::to_string(omp_get_max_threads()));
#endif
    benchmark::AddCustomContext("hyhound_build_time", hyhound_build_time);
    benchmark::AddCustomContext("hyhound_commit_hash", hyhound_commit_hash);
#if defined(__AVX512F__)
    benchmark::AddCustomContext("arch", "avx512f");
#elif defined(__AVX2__)
    benchmark::AddCustomContext("arch", "avx2");
#elif defined(__AVX__)
    benchmark::AddCustomContext("arch", "avx");
#elif defined(__SSE3__)
    benchmark::AddCustomContext("arch", "sse3");
#endif
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
}

using namespace hyhound;
// clang-format off
BM_BLK_REGISTER_F(CholeskyFixture, full_factorize)->Name("full_factorize");
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 1, 32);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 4);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 8);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 12);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 16);
#if __AVX512F__
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 24);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 32);
#endif
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 12, 2);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 4, 12, 4);

BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 8);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 12);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 16);
#if __AVX512F__
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 24);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 32);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 16, 8);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 16, 12);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 16, 16);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 24, 2);
BENCHMARK_BLOCKED(hyh_update, update_cholesky, Downdate, 8, 24, 4);
#endif
// clang-format on

void configure_benchmarks(bool fix_m, int64_t M, int64_t N) {
    std::vector<int64_t> m_range{M}, n_range{N};
    if (fix_m) {
        n_range.clear();
        for (int64_t i = 4; i <= N; i *= 2)
            n_range.push_back(i);
    } else {
        m_range.clear();
        for (int64_t i = 1; i <= std::min<int64_t>(16, M); ++i)
            m_range.push_back(i);
        for (int64_t i = 32; i <= M; i *= 2)
            m_range.push_back(i);
    }
    for (auto *bm : benchmarks)
        bm->ArgNames({"m", "n"})
            ->ArgsProduct({m_range, n_range})
            ->MeasureProcessCPUTime()
            ->UseRealTime();
}
