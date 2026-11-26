## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2026 by the adcc authors
##
## This file is part of adcc.
##
## adcc is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## adcc is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with adcc. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------
"""
Tests for the spin-pure block ("quadruple") approximation of the first-order
pphh_pphh block (order "1q", methods adc2xQ / adc3Q). The reference is the
full first-order block in which only couplings between doubles sharing the
same spatial orbitals {i, j, a, b} are kept, evaluated in numpy.
"""

import numpy as np
import pytest

import adcc
from adcc.adc_pp.matrix import block, shell_degeneracy_threshold
from adcc.AmplitudeVector import AmplitudeVector
from adcc.solver.explicit_symmetrisation import IndexSpinSymmetrisation

pytestmark = pytest.mark.skipif("pyscf" not in adcc.backends.available(), reason="pyscf not found.")

WATER = "O 0 0 0; H 0 0 1.8; H 1.7 0 -0.5"


def _scf(atom=WATER, basis="6-31g"):
    from pyscf import gto, scf

    mol = gto.M(atom=atom, unit="Bohr", basis=basis, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    return mf


def _antisym(x):
    x = 0.5 * (x - x.transpose(1, 0, 2, 3))
    return 0.5 * (x - x.transpose(0, 1, 3, 2))


class _NumpyReference:
    """
    Spin-orbital numpy version of the first order pphh_pphh block, optionally
    masked to the couplings within the same set of degenerate shells.
    """

    def __init__(self, hf):
        ms = hf.mospaces
        self.no, self.nv = ms.n_orbs("o1"), ms.n_orbs("v1")
        self.noa, self.nva = ms.n_orbs_alpha("o1"), ms.n_orbs_alpha("v1")
        self.fo, self.fv = hf.foo.to_ndarray(), hf.fvv.to_ndarray()
        self.oooo = hf.oooo.to_ndarray()
        self.vvvv = hf.vvvv.to_ndarray()
        self.ovov = hf.ovov.to_ndarray()
        so = self.shells(hf, "o1")
        sv = self.shells(hf, "v1")

        def pair_mask(s):
            eq = s[:, None] == s[None, :]
            return np.einsum("pr,qs->pqrs", eq, eq) | np.einsum("ps,qr->pqrs", eq, eq)

        self.m_oo = pair_mask(so)
        self.m_vv = pair_mask(sv)
        # ovov[k, b, j, c]: shell(k) == shell(j) and shell(b) == shell(c)
        self.m_ov = np.einsum("kj,bc->kbjc", so[:, None] == so[None, :], sv[:, None] == sv[None, :])
        spo, spv = np.arange(self.no) // self.noa, np.arange(self.nv) // self.nva
        self.ms0 = (
            spo[:, None, None, None] + spo[None, :, None, None]
            == spv[None, None, :, None] + spv[None, None, None, :]
        )

    @staticmethod
    def shells(hf, space):
        """Shell label of each spin orbital (degenerate within threshold)"""
        n_alpha = hf.mospaces.n_orbs_alpha(space)
        e = hf.orbital_energies(space).to_ndarray()
        assert np.allclose(e[:n_alpha], e[n_alpha:])
        labels = np.zeros(n_alpha, dtype=int)
        for k in range(1, n_alpha):
            if abs(e[k] - e[k - 1]) < shell_degeneracy_threshold:
                labels[k] = labels[k - 1]
            else:
                labels[k] = labels[k - 1] + 1
        return labels[np.arange(len(e)) % n_alpha]

    def apply(self, t, masked):
        m_oo, m_vv, m_ov = (self.m_oo, self.m_vv, self.m_ov) if masked else (1, 1, 1)
        return self.ms0 * (
            +2 * _antisym(np.einsum("ijac,bc->ijab", t, self.fv))
            - 2 * _antisym(np.einsum("ik,kjab->ijab", self.fo, t))
            - 4 * _antisym(np.einsum("ikac,kbjc->ijab", t, self.ovov * m_ov))
            + 0.5 * np.einsum("ijkl,klab->ijab", self.oooo * m_oo, t)
            + 0.5 * np.einsum("ijcd,abcd->ijab", t, self.vvvv * m_vv)
        )

    def random_doubles(self, seed=42):
        rng = np.random.default_rng(seed)
        t = rng.standard_normal((self.no, self.no, self.nv, self.nv)) * self.ms0
        return 4 * _antisym(t)


def _doubles_tensor(gs, data):
    """Doubles tensor without spin-block symmetry holding data."""
    tt = adcc.guess_zero(adcc.AdcMatrix("adc2x", gs), spin_block_symmetrisation="none").pphh
    tt.set_from_ndarray(data, 1e-12)
    return tt


@pytest.mark.parametrize("atom,basis", [(WATER, "6-31g"), ("N 0 0 0; N 0 0 2.07", "cc-pvdz")])
def test_quadruple_block_matches_masked_full_block(atom, basis):
    hf = adcc.ReferenceState(_scf(atom, basis))
    ref = _NumpyReference(hf)
    t = ref.random_doubles()
    gs = adcc.LazyMp(hf)
    tt = _doubles_tensor(gs, t)  # no spin symmetry: tests all spin blocks

    res = block(gs, ["pphh", "pphh"], "1q").apply(AmplitudeVector(pphh=tt))
    expected = ref.apply(t, masked=True)
    np.testing.assert_allclose(res.pphh.to_ndarray(), expected, atol=1e-11)
    # the diagonal alone is not sufficient
    diag = block(gs, ["pphh", "pphh"], "1d").apply(AmplitudeVector(pphh=tt))
    assert np.max(np.abs(diag.pphh.to_ndarray() - expected)) > 1e-2


@pytest.mark.parametrize("spin_block_symmetrisation", ["symmetric", "antisymmetric"])
def test_cvs_quadruple_block_matches_masked_full_block(spin_block_symmetrisation):
    mf = _scf()
    hf = adcc.ReferenceState(mf)
    hfc = adcc.ReferenceState(mf, core_orbitals=1)
    ref = _NumpyReference(hf)
    msc = hfc.mospaces
    n1, n1a = msc.n_orbs("o1"), msc.n_orbs_alpha("o1")
    nc, nca = msc.n_orbs("o2"), msc.n_orbs_alpha("o2")
    nv, nva = ref.nv, ref.nva
    # map CVS orbitals to the general occupied spin orbitals (core = lowest)
    g_o1 = np.array([(p // n1a) * ref.noa + nca + p % n1a for p in range(n1)])
    g_o2 = np.array([(p // nca) * ref.noa + p % nca for p in range(nc)])
    s1, s2, sv = (np.arange(n1) // n1a, np.arange(nc) // nca, np.arange(nv) // nva)
    ms0 = (
        s1[:, None, None, None] + s2[None, :, None, None]
        == sv[None, None, :, None] + sv[None, None, None, :]
    )

    rng = np.random.default_rng(1)
    tc = rng.standard_normal((n1, nc, nv, nv)) * ms0
    tc = 0.5 * (tc - tc.transpose(0, 1, 3, 2))
    flipped = tc
    for ax in range(4):
        flipped = np.roll(flipped, flipped.shape[ax] // 2, axis=ax)
    sign = 1 if spin_block_symmetrisation == "symmetric" else -1
    tc = 0.5 * (tc + sign * flipped)

    tf = np.zeros((ref.no, ref.no, nv, nv))
    tf[np.ix_(g_o1, g_o2, range(nv), range(nv))] = tc
    tf[np.ix_(g_o2, g_o1, range(nv), range(nv))] = -tc.transpose(1, 0, 2, 3)
    expected = ref.apply(tf, masked=True)[np.ix_(g_o1, g_o2, range(nv), range(nv))]

    gsc = adcc.LazyMp(hfc)
    tt = adcc.guess_zero(
        adcc.AdcMatrix("cvs-adc2x", gsc), spin_block_symmetrisation=spin_block_symmetrisation
    ).pphh
    tt.set_from_ndarray(tc, 1e-12)
    res = block(gsc, ["pphh", "pphh"], "1q", variant="cvs").apply(AmplitudeVector(pphh=tt))
    np.testing.assert_allclose(res.pphh.to_ndarray(), expected * ms0, atol=1e-12)


def test_quadruple_block_commutes_with_s2():
    """[M, S^2] = 0 on Ms = 0 doubles for the masked block, not for the diagonal"""
    hf = adcc.ReferenceState(_scf())
    ref = _NumpyReference(hf)

    def shift(x, up):
        # S_+: virtual beta -> alpha (+), occupied alpha -> beta (-)
        # S_-: virtual alpha -> beta (+), occupied beta -> alpha (-)
        out = np.zeros_like(x)
        for ax in range(4):
            n = x.shape[ax]
            occ = ax < 2
            lo, hi = slice(0, n // 2), slice(n // 2, n)
            src, dst = (hi, lo) if up != occ else (lo, hi)
            sl_src, sl_dst = [slice(None)] * 4, [slice(None)] * 4
            sl_src[ax], sl_dst[ax] = src, dst
            out[tuple(sl_dst)] += (-1.0 if occ else 1.0) * x[tuple(sl_src)]
        return out

    def s2(x):  # Ms = 0: S^2 = S_- S_+
        return shift(shift(x, True), False)

    gs = adcc.LazyMp(hf)
    t = ref.random_doubles()
    for order, commutes in [("1q", True), ("1d", False)]:
        bl = block(gs, ["pphh", "pphh"], order)

        def apply(x, bl=bl):
            tt = _doubles_tensor(gs, x)
            return bl.apply(AmplitudeVector(pphh=tt)).pphh.to_ndarray() * ref.ms0

        comm = apply(s2(t)) - s2(apply(t))
        rel = np.linalg.norm(comm) / np.linalg.norm(apply(t))
        assert (rel < 1e-12) == commutes, (order, rel)


@pytest.mark.parametrize("method", ["adc2xQ", "adc3Q"])
def test_quadruple_matrix_preserves_singlets(method):
    hf = adcc.ReferenceState(_scf())
    matrix = adcc.AdcMatrix(method, adcc.LazyMp(hf))
    guesses = adcc.guesses_singlet(matrix, n_guesses=4, block="pphh") + adcc.guesses_singlet(
        matrix, n_guesses=2, block="ph"
    )
    x = sum((0.3 + 0.1 * k) * g for k, g in enumerate(guesses))
    y = matrix @ x
    sym = IndexSpinSymmetrisation(matrix, enforce_spin_kind="singlet")
    ye = sym.symmetrise([y.copy()])[0]
    assert np.sqrt((ye - y).dot(ye - y)) < 1e-12


@pytest.mark.parametrize("kind", ["singlet", "triplet"])
@pytest.mark.parametrize("method,core_orbitals", [("adc2xQ", None), ("cvs-adc2xQ", 1)])
def test_quadruple_eigenpairs_converge(method, core_orbitals, kind):
    """Davidson must converge to true eigenpairs (spin maps propagate)."""
    kwargs = {f"n_{kind}s": 3}
    if core_orbitals:
        kwargs["core_orbitals"] = core_orbitals
    states = adcc.run_adc(_scf(), method=method, conv_tol=1e-9, output=None, **kwargs)
    for v, w in zip(states.excitation_vector, states.excitation_energy):
        r = states.matrix @ v - w * v
        assert np.sqrt(r.dot(r)) < 1e-7


@pytest.mark.parametrize(
    "atom,basis,degenerate_sets",
    [
        # Pi states of N2 (pairs)
        ("N 0 0 0; N 0 0 2.07", "cc-pvdz", [(0, 1), (3, 4)]),
        # T2 states of CH4 (triples)
        (
            "C 0 0 0; H 1.19 1.19 1.19; H -1.19 -1.19 1.19; H -1.19 1.19 -1.19; H 1.19 -1.19 -1.19",
            "6-31g",
            [(0, 1, 2), (3, 4, 5)],
        ),
        # slightly broken symmetry: orbital splittings ~1e-7 must still form shells
        (
            "C 0 0 0; H 1.19 1.19 1.190001; H -1.19 -1.19 1.19; H -1.19 1.19 -1.19; "
            + "H 1.19 -1.19 -1.19",
            "6-31g",
            [(0, 1, 2), (3, 4, 5)],
        ),
    ],
)
def test_quadruple_degenerate_states(atom, basis, degenerate_sets):
    """With degenerate orbitals, degenerate states have to stay degenerate"""
    states = adcc.run_adc(
        _scf(atom, basis), method="adc2xQ", n_singlets=6, conv_tol=1e-9, output=None
    )
    e = states.excitation_energy
    for group in degenerate_sets:
        # without shells the splitting would be ~1e-4 (meV)
        assert np.ptp(e[list(group)]) < 1e-6, (group, e)


def test_quadruple_unrestricted_not_implemented():
    from pyscf import gto, scf

    mol = gto.M(atom=WATER, unit="Bohr", basis="6-31g", verbose=0)
    mf = scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    with pytest.raises(NotImplementedError):
        adcc.run_adc(mf, method="adc2xQ", n_states=1, output=None)


@pytest.mark.parametrize("method", ["adc2xD", "adc2xQ", "adc3Q"])
def test_drop_vvvv(method):
    """vvvv is dropped by default, changes nothing and is not re-imported"""
    mf = _scf()
    hf_keep = adcc.ReferenceState(mf)
    ref = adcc.run_adc(
        hf_keep, method=method, n_singlets=3, conv_tol=1e-9, output=None, drop_vvvv=False
    )
    assert "v1v1v1v1" in hf_keep.cached_eri_blocks
    hf = adcc.ReferenceState(mf)
    states = adcc.run_adc(hf, method=method, n_singlets=3, conv_tol=1e-9, output=None)
    assert "v1v1v1v1" not in hf.cached_eri_blocks
    np.testing.assert_allclose(states.excitation_energy, ref.excitation_energy, atol=1e-10)


@pytest.mark.parametrize("method", ["adc2x", "adc3"])
def test_drop_vvvv_safeguard(method):
    """Methods which need vvvv in the matrix-vector products keep it"""
    hf = adcc.ReferenceState(_scf())
    adcc.run_adc(hf, method=method, n_singlets=1, output=None)
    assert "v1v1v1v1" in hf.cached_eri_blocks
    with pytest.raises(adcc.exceptions.InputError):
        adcc.run_adc(hf, method=method, n_singlets=1, output=None, drop_vvvv=True)
