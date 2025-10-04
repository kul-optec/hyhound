import hyhound
import numpy as np
import numpy.linalg as la


def test_update_cholesky():
    rng = np.random.default_rng(seed=123)
    m, n = 7, 13
    L = np.tril(rng.uniform(-2, 2, (n, n)))
    A = rng.uniform(-1, 1, (n, m))
    L̃, _ = hyhound.update_cholesky(L, A)
    assert la.norm(L @ L.T + A @ A.T - L̃ @ L̃.T, "fro") < 1e-12
    L2, _ = hyhound.downdate_cholesky(L̃, A)
    assert la.norm(L @ L.T - L2 @ L2.T, "fro") < 1e-12


def test_update_cholesky_tall():
    rng = np.random.default_rng(seed=123)
    m, n, p = 7, 13, 43
    L = np.tril(rng.uniform(-2, 2, (p, n)))
    A = rng.uniform(-1, 1, (p, m))
    L̃, Ã = hyhound.update_cholesky(L, A)
    Ã[:n, :] = 0
    assert la.norm(L @ L.T + A @ A.T - L̃ @ L̃.T - Ã @ Ã.T, "fro") < 1e-12
    L2, A2 = hyhound.downdate_cholesky(L̃, A)
    A2[:n, :] = 0
    assert la.norm(L @ L.T - L2 @ L2.T - Ã @ Ã.T + A2 @ A2.T, "fro") < 1e-12


def test_update_cholesky_inplace():
    rng = np.random.default_rng(seed=123)
    m, n = 7, 13
    L = np.tril(rng.uniform(-2, 2, (n, n)))
    A = rng.uniform(-1, 1, (n, m))
    hyhound.update_cholesky_inplace(L̃ := L.copy(order="F"), A.copy(order="F"))
    assert la.norm(L @ L.T + A @ A.T - L̃ @ L̃.T, "fro") < 1e-12
    hyhound.downdate_cholesky_inplace(L̃, A.copy(order="F"))
    assert la.norm(L @ L.T - L̃ @ L̃.T, "fro") < 1e-12


def test_update_cholesky_diag():
    rng = np.random.default_rng(seed=123)
    m, n = 7, 13
    L = np.tril(rng.uniform(-2, 2, (n, n)))
    L += 10 * np.eye(n)
    A = rng.uniform(-1, 1, (n, m))
    d = rng.uniform(-1, 1, (m,))
    L̃, _ = hyhound.update_cholesky_diag(L, A, d)
    assert la.norm(L @ L.T + A @ np.diag(d) @ A.T - L̃ @ L̃.T, "fro") < 1e-12
    L2, _ = hyhound.update_cholesky_diag(L̃, A, -d)
    assert la.norm(L @ L.T - L2 @ L2.T, "fro") < 1e-12


def test_update_cholesky_sign():
    rng = np.random.default_rng(seed=123)
    m, n = 7, 13
    L = np.tril(rng.uniform(-2, 2, (n, n)))
    L += 10 * np.eye(n)
    A = rng.uniform(-1, 1, (n, m))
    d = rng.uniform(-1, 1, (m,))
    s = 0 * d  # np.copysign(np.zeros((m,)), d)
    e = np.copysign(np.ones((m,)), d)
    L̃, _ = hyhound.update_cholesky_sign(L, A, s)
    assert la.norm(L @ L.T + A @ np.diag(e) @ A.T - L̃ @ L̃.T, "fro") < 1e-12
    L2, _ = hyhound.update_cholesky_sign(L̃, A, -s)
    assert la.norm(L @ L.T - L2 @ L2.T, "fro") < 1e-12
