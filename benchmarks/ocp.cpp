#include <benchmark/benchmark.h>
#include <guanaqo/openmp.h>
#include <hyhound-version.h>
#include <utility>

#include <hyhound/ocp/riccati.hpp>
#include <hyhound/ocp/schur.hpp>

using namespace hyhound;
using namespace hyhound::ocp;

auto generate_ocp() {
    using std::exp2;
    std::mt19937 rng{321};
    std::normal_distribution<real_t> nrml{0, 10};

    // OCPDataRiccati ocp{.N = 20, .nx = 240, .nu = 80, .ny = 240}; // large
    OCPDataRiccati ocp{.N = 20, .nx = 24, .nu = 8, .ny = 24}; // medium
    // OCPDataRiccati ocp{.N = 20, .nx = 6, .nu = 2, .ny = 6}; // small
    ocp.init_random(123);

    mat Σ = mat::Zero(ocp.ny, ocp.N + 1);
    std::ranges::generate(Σ.reshaped(), [&] { return exp2(nrml(rng)); });
    return std::pair{std::move(ocp), std::move(Σ)};
}

void bm_factor_riccati(benchmark::State &state) {
    auto [ocp, Σ] = generate_ocp();
    RiccatiFactor fac{.ocp = ocp};
    for (auto _ : state)
        factor(fac, Σ);
}

void bm_update_riccati(benchmark::State &state) {
    using std::exp2;
    std::mt19937 rng{54321};
    std::normal_distribution<real_t> nrml{0, 10};
    std::bernoulli_distribution bern{static_cast<double>(state.range(0)) / 100};

    auto [ocp, Σ] = generate_ocp();
    RiccatiFactor fac{.ocp = ocp};
    mat ΔΣ = Σ;
    for (auto _ : state) {
        state.PauseTiming();
        factor(fac, Σ);
        std::ranges::generate(ΔΣ.reshaped(), [&] {
            return bern(rng) ? exp2(nrml(rng)) : real_t{0};
        });
        state.ResumeTiming();
        update(fac, ΔΣ);
    }
}

void bm_factor_schur(benchmark::State &state) {
    auto [ocp, Σ] = generate_ocp();
    auto ocp_sch  = OCPDataSchur::from_riccati(ocp);
    SchurFactor factor_sch{.ocp = ocp_sch};
    for (auto _ : state)
        factor(factor_sch, Σ);
}

void bm_update_schur(benchmark::State &state) {
    using std::exp2;
    std::mt19937 rng{54321};
    std::normal_distribution<real_t> nrml{0, 10};
    std::bernoulli_distribution bern{static_cast<double>(state.range(0)) / 100};

    auto [ocp, Σ] = generate_ocp();
    auto ocp_sch  = OCPDataSchur::from_riccati(ocp);
    SchurFactor factor_sch{.ocp = ocp_sch};
    mat ΔΣ = Σ;

    for (auto _ : state) {
        state.PauseTiming();
        factor(factor_sch, Σ);
        std::ranges::generate(ΔΣ.reshaped(), [&] {
            return bern(rng) ? exp2(nrml(rng)) : real_t{0};
        });
        state.ResumeTiming();
        update(factor_sch, ΔΣ);
    }
}

BENCHMARK(bm_factor_riccati);
BENCHMARK(bm_update_riccati)->ArgNames({"prcnt"})->DenseRange(1, 20, 1);
BENCHMARK(bm_factor_schur);
BENCHMARK(bm_update_schur)->ArgNames({"prcnt"})->DenseRange(1, 20, 1);

int main(int argc, char **argv) {
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
        return 1;
#if GUANAQO_WITH_OPENMP
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
