#pragma once

#include <hyhound/ocp/riccati.hpp>

namespace hyhound::ocp {

struct HYHOUND_OCP_EXPORT OCPDataSchur {
    index_t N  = 31;
    index_t nx = 40, nu = 10, ny = 10;
    mat ABs  = mat::Zero(nx, (nx + nu) * N);
    mat CDs  = mat::Zero(ny, (nx + nu) * (N + 1));
    mat QSRs = mat::Zero(nx + nu, (nx + nu) * (N + 1));
    mat qrs  = mat::Zero(nx + nu, N + 1);
    mat es   = mat::Zero(nx, N + 1);

    auto F(index_t j) { return ABs.middleCols((nu + nx) * j, nu + nx); }
    auto G(index_t j) { return CDs.middleCols((nu + nx) * j, nu + nx); }
    auto B(index_t j) { return ABs.middleCols((nu + nx) * j + nx, nu); }
    auto A(index_t j) { return ABs.middleCols((nu + nx) * j, nx); }
    auto D(index_t j) { return CDs.middleCols((nu + nx) * j + nx, nu); }
    auto C(index_t j) { return CDs.middleCols((nu + nx) * j, nx); }
    auto H(index_t j) { return QSRs.middleCols((nu + nx) * j, nu + nx); }
    auto Q(index_t j) { return H(j).topLeftCorner(nx, nx); }
    auto R(index_t j) { return H(j).bottomRightCorner(nu, nu); }
    auto S(index_t j) { return H(j).bottomLeftCorner(nu, nx); }
    auto q(index_t j) { return qrs.col(j).topRows(nx); }
    auto r(index_t j) { return qrs.col(j).bottomRows(nu); }
    auto F(index_t j) const { return ABs.middleCols((nu + nx) * j, nu + nx); }
    auto G(index_t j) const { return CDs.middleCols((nu + nx) * j, nu + nx); }
    auto B(index_t j) const { return ABs.middleCols((nu + nx) * j + nx, nu); }
    auto A(index_t j) const { return ABs.middleCols((nu + nx) * j, nx); }
    auto D(index_t j) const { return CDs.middleCols((nu + nx) * j + nx, nu); }
    auto C(index_t j) const { return CDs.middleCols((nu + nx) * j, nx); }
    auto H(index_t j) const { return QSRs.middleCols((nu + nx) * j, nu + nx); }
    auto Q(index_t j) const { return H(j).topLeftCorner(nx, nx); }
    auto R(index_t j) const { return H(j).bottomRightCorner(nu, nu); }
    auto S(index_t j) const { return H(j).bottomLeftCorner(nu, nx); }
    auto q(index_t j) const { return qrs.col(j).topRows(nx); }
    auto r(index_t j) const { return qrs.col(j).bottomRows(nu); }

    static OCPDataSchur from_riccati(const OCPDataRiccati &ric) {
        OCPDataSchur schur{
            .N = ric.N, .nx = ric.nx, .nu = ric.nu, .ny = ric.ny};
        for (index_t i = 0; i < schur.N; ++i) {
            schur.A(i) = ric.A(i);
            schur.B(i) = ric.B(i);
            schur.C(i) = ric.C(i);
            schur.D(i) = ric.D(i);
            schur.Q(i) = ric.Q(i);
            schur.S(i) = ric.S(i).transpose();
            schur.R(i) = ric.R(i);
            schur.r(i) = ric.r(i);
            schur.q(i) = ric.q(i);
        }
        schur.C(schur.N) = ric.C(schur.N);
        schur.Q(schur.N) = ric.Q(schur.N);
        schur.q(schur.N) = ric.q(schur.N);
        schur.es         = ric.es;
        return schur;
    }
};

struct HYHOUND_OCP_EXPORT SchurFactor {
    const OCPDataSchur &ocp;
    mat LHs  = mat::Zero(ocp.nx + ocp.nu, (ocp.nx + ocp.nu) * (ocp.N + 1));
    mat LΨds = mat::Zero(ocp.nx, ocp.nx *(ocp.N + 1));
    mat LΨss = mat::Zero(ocp.nx, ocp.nx *ocp.N);
    mat Vs   = mat::Zero(ocp.nx, (ocp.nx + ocp.nu) * ocp.N);
    mat Wᵀs  = mat::Zero(ocp.nx + ocp.nu, ocp.nx *(ocp.N + 1));
    mat ΣG   = mat::Zero(ocp.ny, ocp.nx + ocp.nu);
    mat v    = mat::Zero(ocp.nx + ocp.nu, ocp.N + 1);
    mat λ    = mat::Zero(ocp.nx, ocp.N + 1);
    vec e̅    = vec::Zero(ocp.nx + ocp.nu);
    vec ẽ    = vec::Zero(ocp.nx + ocp.nu);
    vec ψ    = vec::Zero(2 * ocp.nx);

    auto LH(index_t j) {
        return LHs.middleCols((ocp.nx + ocp.nu) * j, ocp.nx + ocp.nu);
    }
    auto LHxx(index_t j) { return LH(j).topLeftCorner(ocp.nx, ocp.nx); }
    auto LΨd(index_t j) { return LΨds.middleCols(ocp.nx * j, ocp.nx); }
    auto LΨs(index_t j) { return LΨss.middleCols(ocp.nx * j, ocp.nx); }
    auto V(index_t j) {
        return Vs.middleCols((ocp.nx + ocp.nu) * j, ocp.nx + ocp.nu);
    }
    auto Wᵀ(index_t j) { return Wᵀs.middleCols(ocp.nx * j, ocp.nx); }
    auto x(index_t j) { return v.col(j).topRows(ocp.nx); }
    auto u(index_t j) { return v.col(j).bottomRows(ocp.nu); }
};

HYHOUND_OCP_EXPORT void factor(SchurFactor &factor, Eigen::Ref<const mat> Σ);
HYHOUND_OCP_EXPORT void update(SchurFactor &factor, Eigen::Ref<const mat> ΔΣ);
HYHOUND_OCP_EXPORT void solve(SchurFactor &factor);

} // namespace hyhound::ocp
