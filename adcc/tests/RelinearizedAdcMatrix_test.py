#!/usr/bin/env python3
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
import io

import numpy as np
import pytest

import adcc
from adcc.AdcMatrix import AdcExtraTerm, AdcMatrix, RelinearizedAdcMatrix
from adcc.adc_pp.matrix import AdcBlock

from .testdata_cache import testdata_cache


def adc2_matrix(system):
    refstate = testdata_cache.refstate(system, case="gen")
    return adcc.AdcMatrix("adc2", adcc.LazyMp(refstate))


@pytest.mark.parametrize("system,tol", [("h2o_sto3g", 1e-6),
                                        ("h2o_def2tzvp", 1e-4)])
class TestRelinearizedAdc2:
    def test_finalized_states_match_exact(self, system, tol):
        matrix = adc2_matrix(system)
        exact = adcc.run_adc(matrix, n_singlets=3, conv_tol=1e-8,
                             output=None)
        raw = adcc.run_adc(matrix, n_singlets=3, conv_tol=1e-8, relin=True,
                           finalize_relinearized=False, output=None)
        final = adcc.run_adc(matrix, n_singlets=3, conv_tol=1e-8, relin=True,
                             output=None)

        raw_error = np.abs(raw.excitation_energy - exact.excitation_energy)
        final_error = np.abs(final.excitation_energy - exact.excitation_energy)
        assert np.all(final_error < tol)
        assert np.max(raw_error) > 10 * np.max(final_error)

        for v, ref in zip(final.excitation_vector, exact.excitation_vector):
            assert not any(getattr(v, b).needs_evaluation for b in v.blocks)
            assert v.dot(v) == pytest.approx(1.0)
            assert abs(v.dot(ref)) > 0.9999
            assert v.pphh.dot(v.pphh) == pytest.approx(ref.pphh.dot(ref.pphh),
                                                       rel=1e-2)

    def test_finalize_energy_is_rayleigh_quotient(self, system, tol):
        matrix = adc2_matrix(system)
        state = adcc.run_adc(matrix, n_singlets=2, conv_tol=1e-8, relin=True,
                             finalize_relinearized=False, output=None)
        rel = state.matrix
        for ampl, omega in zip(state.excitation_vector,
                               state.excitation_energy):
            energy, v, n_steps = rel.finalize(ampl, omega)
            rayleigh = v.dot(AdcMatrix.matvec(rel, v))  # real matrix
            assert energy == pytest.approx(rayleigh, abs=1e-10)
            assert 1 <= n_steps <= 5

    def test_raw_vectors_lack_eliminated_part(self, system, tol):
        matrix = adc2_matrix(system)
        state = adcc.run_adc(matrix, n_singlets=2, conv_tol=1e-8, relin=True,
                             finalize_relinearized=False, output=None)
        inactive = state.matrix._mask_inactive_tensor
        for v in state.excitation_vector:
            eliminated = v.pphh * inactive
            assert eliminated.dot(eliminated) == pytest.approx(0.0, abs=1e-24)


@pytest.mark.parametrize("system,kind,n_states,guess_energies,n_reference,tol", [
    ("cn_sto3g", "any", 4, 1.0, 6, 2e-3),
    ("h2o_def2tzvp", "singlet", 12, 1.5, 12, 1e-3),
])
def test_states_sorted_after_finalize(system, kind, n_states, guess_energies,
                                      n_reference, tol):
    # A guess energy far above the states makes their corrections differ
    # enough to swap close-lying states before sorting (the energies are
    # correspondingly less accurate).
    matrix = adc2_matrix(system)
    exact = adcc.run_adc(matrix, n_states=n_reference, kind=kind,
                         conv_tol=1e-8, output=None)
    output = io.StringIO()
    final = adcc.run_adc(matrix, n_states=n_states, kind=kind, conv_tol=1e-8,
                         relin=True, guess_energies=guess_energies,
                         output=output)
    assert "State order changed" in output.getvalue()
    assert np.all(np.diff(final.excitation_energy) >= 0)
    assert final.excitation_energy == pytest.approx(
        exact.excitation_energy[:n_states], abs=tol)


@pytest.mark.parametrize("case", ["gen", "fc-fv"])
def test_spin_flip(case):
    refstate = testdata_cache.refstate("hf_631g", case=case)
    matrix = adcc.AdcMatrix("adc2", adcc.LazyMp(refstate))
    exact = adcc.run_adc(matrix, n_spin_flip=4, conv_tol=1e-8, output=None)
    final = adcc.run_adc(matrix, n_spin_flip=2, conv_tol=1e-8, relin=True,
                         output=None)
    assert final.excitation_energy == pytest.approx(
        exact.excitation_energy[:2], abs=1e-4)


def test_guesses_require_guess_energies():
    matrix = adc2_matrix("h2o_sto3g")
    exact = adcc.run_adc(matrix, n_singlets=2, conv_tol=1e-8, output=None)
    with pytest.raises(adcc.InputError):
        adcc.run_adc(matrix, n_singlets=2, relin=True,
                     guesses=exact.excitation_vector, output=None)


@pytest.mark.parametrize("system,guess_energies", [("h2o_sto3g", 0.5),
                                                   ("h2o_def2tzvp", 0.35)])
def test_user_guesses_with_eliminated_components(system, guess_energies):
    # Exact ADC(2) vectors have components in the eliminated doubles; unless
    # projected away, these converge to spurious zero eigenvalues.
    matrix = adc2_matrix(system)
    exact = adcc.run_adc(matrix, n_singlets=3, conv_tol=1e-8, output=None)
    guesses = [v.copy() for v in exact.excitation_vector]
    final = adcc.run_adc(matrix, n_singlets=3, conv_tol=1e-8, relin=True,
                         guesses=guesses, guess_energies=guess_energies,
                         output=None)
    assert final.excitation_energy == pytest.approx(exact.excitation_energy,
                                                    abs=1e-4)


def test_states_above_lowest_eliminated_configuration_kept_raw():
    # Without explicit configurations, the doubles-dominated exact states 9
    # and 10 (just above the lowest eliminated D0) are not part of the fold;
    # its 9th and 10th eigenvalues lie far above and must be kept as they are.
    matrix = adc2_matrix("h2o_sto3g")
    kwargs = dict(n_singlets=10, n_guesses=10, conv_tol=1e-8, relin=True,
                  guess_energies=0.5, explicit_width=0.0, output=None)
    exact = adcc.run_adc(matrix, n_singlets=10, n_guesses=10, conv_tol=1e-8,
                         output=None)
    raw = adcc.run_adc(matrix, finalize_relinearized=False, **kwargs)
    with pytest.warns(UserWarning, match="not described by the fold"):
        final = adcc.run_adc(matrix, **kwargs)
    pole = final.matrix._lowest_eliminated_d0
    above = raw.excitation_energy >= pole
    assert np.count_nonzero(above) == 2
    assert final.excitation_energy[above] == pytest.approx(
        raw.excitation_energy[above])
    assert final.excitation_energy[~above] == pytest.approx(
        exact.excitation_energy[~above], abs=1e-3)


def relinearized_s1(guess_energies):
    """x_P, b and raw omega of the water/STO-3G S1 without explicit
    configurations, for a guess energy far above it (omega < root)."""
    matrix = adc2_matrix("h2o_sto3g")
    state = adcc.run_adc(matrix, n_singlets=1, conv_tol=1e-8, relin=True,
                         guess_energies=guess_energies, explicit_width=0.0,
                         finalize_relinearized=False, output=None)
    rel = state.matrix
    x_p = rel._active_masked(state.excitation_vector[0])
    b = rel._coupling_into_eliminated(x_p)
    return rel, state.excitation_vector[0], x_p, b, state.excitation_energy[0]


@pytest.mark.parametrize("guess_energies", [1.0, 1.5])
def test_newton_bracketed_by_lowest_pole(guess_energies):
    rel, _, x_p, b, omega = relinearized_s1(guess_energies)
    root, _, _ = rel._self_consistent_energy(x_p, b, omega, tol=1e-12)
    assert omega < root

    # A (fictitious) pole just above the root clamps the start to the left
    # of the root, so the first Newton step overshoots the pole and must be
    # replaced by bisection, without changing the root.
    rel._lowest_eliminated_d0 = root + 1e-7
    output = io.StringIO()
    root_bracketed, _, _ = rel._self_consistent_energy(x_p, b, omega,
                                                       tol=1e-12,
                                                       output=output)
    assert "(bisection)" in output.getvalue()
    assert root_bracketed == pytest.approx(root, abs=1e-10)


def test_raw_pair_kept_if_iteration_ends_at_pole():
    # A (fictitious) pole between omega and the root, where b does not
    # couple: no root below it, so the iteration ends at the pole.
    rel, ampl, x_p, b, omega = relinearized_s1(1.5)
    root, _, _ = rel._self_consistent_energy(x_p, b, omega, tol=1e-12)
    rel._lowest_eliminated_d0 = 0.5 * (omega + root)
    with pytest.warns(UserWarning, match="ends at the lowest eliminated"):
        energy, vector, n_steps = rel.finalize(ampl, omega)
    assert energy == omega and vector is ampl and n_steps == 0


@pytest.mark.parametrize("case,method,n_singlets", [
    ("gen", "adc2", 12),     # more states than singles guesses
    ("gen", "adc2", 6),      # fewer guesses than requested by default
    ("fv", "adc2", 3),
    ("cvs", "cvs-adc2", 2),
    ("cvs", "cvs-adc2", 3),
])
def test_few_guesses(case, method, n_singlets):
    refstate = testdata_cache.refstate("h2o_sto3g", case=case)
    matrix = adcc.AdcMatrix(method, adcc.LazyMp(refstate))
    exact = adcc.run_adc(matrix, n_singlets=n_singlets, conv_tol=1e-8,
                         output=None)
    final = adcc.run_adc(matrix, n_singlets=n_singlets, conv_tol=1e-8,
                         relin=True, output=None)
    assert final.excitation_energy == pytest.approx(exact.excitation_energy,
                                                    abs=1e-4)


def test_fold_finite_at_explicit_diagonal():
    # omega hitting the diagonal of an explicit configuration exactly (e.g. a
    # pure explicit doubles state) must not give 0/0 in the fold.
    matrix = adc2_matrix("h2o_sto3g")
    rel = RelinearizedAdcMatrix(matrix, 0.5, explicit_width=2.0)
    active = rel._mask_active_tensor.to_ndarray() > 0.5
    omega = rel.diagonal_0.pphh.to_ndarray()[active][0]  # numpy scalar

    v = adcc.guess_zero(rel)
    v.set_random()
    x_q = rel.inverter.apply(rel._coupling_into_eliminated(v), omega)
    assert np.all(np.isfinite(x_q.pphh.to_ndarray()))


@pytest.mark.parametrize("system,kind,n_states,n_reference", [
    ("h2o_sto3g", "triplet", 3, 3),
    ("h2o_def2tzvp", "triplet", 3, 3),
    # The full solver needs more states to find the lowest 4 for CN
    ("cn_sto3g", "any", 4, 6),
    ("cn_ccpvdz", "any", 3, 3),
])
def test_other_kinds(system, kind, n_states, n_reference):
    matrix = adc2_matrix(system)
    exact = adcc.run_adc(matrix, n_states=n_reference, kind=kind,
                         conv_tol=1e-8, output=None)
    final = adcc.run_adc(matrix, n_states=n_states, kind=kind, conv_tol=1e-8,
                         relin=True, output=None)
    assert final.excitation_energy == pytest.approx(
        exact.excitation_energy[:n_states], abs=1e-4)


@pytest.mark.parametrize("system,case", [
    ("h2o_sto3g", "gen"), ("h2o_sto3g", "fc"), ("h2o_sto3g", "fv"),
    ("h2o_sto3g", "cvs"), ("cn_sto3g", "gen"), ("cn_sto3g", "fc-fv"),
])
@pytest.mark.parametrize("omega,width", [(0.5, 1.0), (0.8, 0.0), (1.0, 1.5),
                                         ([0.3, 0.5], 0.5)])
def test_masks_match_dense_selection(system, case, omega, width):
    method = "cvs-adc2" if "cvs" in case else "adc2"
    refstate = testdata_cache.refstate(system, case=case)
    matrix = adcc.AdcMatrix(method, adcc.LazyMp(refstate))
    rel = RelinearizedAdcMatrix(matrix, omega, explicit_width=width)

    diagonal = rel.diagonal_0.pphh.to_ndarray()
    active = rel._mask_active_tensor.to_ndarray() > 0.5
    inactive = rel._mask_inactive_tensor.to_ndarray() > 0.5
    np.testing.assert_array_equal(active,
                                  diagonal <= rel.omega_screen + width)
    np.testing.assert_array_equal(inactive, ~active)

    # The pole is the lowest eliminated element (conservatively including
    # those vanishing by antisymmetry, e.g. i == j).
    eliminated = diagonal[inactive]
    if eliminated.size:
        assert rel._lowest_eliminated_d0 == pytest.approx(eliminated.min())
    else:
        assert rel._lowest_eliminated_d0 == np.inf


def test_ignored_parameters_warn():
    matrix = adc2_matrix("h2o_sto3g")
    with pytest.warns(UserWarning, match="relin is not set"):
        adcc.run_adc(matrix, n_singlets=1, guess_energies=0.5, output=None)
    rel = RelinearizedAdcMatrix(matrix, 0.5)
    with pytest.warns(UserWarning, match="already a RelinearizedAdcMatrix"):
        adcc.run_adc(rel, n_singlets=1, explicit_width=0.2, output=None)
    # order1_width is not used (nor validated) for a 0th order block
    RelinearizedAdcMatrix(matrix, 0.5, order1_width=0.0)


def test_explicit_guess_count_reduced_with_warning():
    # Water/STO-3G has 10 singlet singles and, with the default window, no
    # explicit doubles to offer as guesses.
    matrix = adc2_matrix("h2o_sto3g")
    with pytest.warns(UserWarning, match="Only 10 guesses available"):
        adcc.run_adc(matrix, n_singlets=3, n_guesses=12, relin=True,
                     output=None)


def test_missing_doubles_guesses_filled_with_singles():
    from adcc.workflow import obtain_guesses_by_inspection

    matrix = adc2_matrix("h2o_sto3g")
    state = adcc.run_adc(matrix, n_singlets=3, conv_tol=1e-8, relin=True,
                         output=None)
    with pytest.warns(UserWarning, match="Only 0 doubles guesses"):
        guesses = obtain_guesses_by_inspection(state.matrix, 6, "singlet",
                                               n_guesses_doubles=2,
                                               min_guesses=3)
    assert len(guesses) == 6


def test_extra_terms():
    matrix = adc2_matrix("h2o_sto3g")
    ones = matrix.diagonal().ones_like()

    def shift_ph(hf, mp, intermediates):
        return AdcBlock(lambda v: adcc.AmplitudeVector(ph=0.1 * v.ph),
                        adcc.AmplitudeVector(ph=0.1 * ones.ph))
    term = AdcExtraTerm(matrix, {"ph_ph": shift_ph})

    # Terms on the input matrix would be dropped by the rebuild.
    with pytest.raises(NotImplementedError):
        RelinearizedAdcMatrix(matrix + term, 0.5)
    rel = RelinearizedAdcMatrix(matrix, 0.5)
    with pytest.raises(NotImplementedError):
        rel + term
    with pytest.raises(NotImplementedError):
        term + rel

    # A term on the complementary block would be ignored by the fold.
    def shift_pphh(hf, mp, intermediates):
        return AdcBlock(lambda v: adcc.AmplitudeVector(pphh=0.1 * v.pphh),
                        adcc.AmplitudeVector(pphh=0.1 * ones.pphh))
    with pytest.raises(NotImplementedError):
        rel += AdcExtraTerm(matrix, {"pphh_pphh": shift_pphh})

    # In place (as run_adc adds environment terms) the term is applied.
    exact = adcc.run_adc(matrix + term, n_singlets=2, conv_tol=1e-8,
                         output=None)
    rel += term
    final = adcc.run_adc(rel, n_singlets=2, conv_tol=1e-8, output=None)
    assert final.excitation_energy == pytest.approx(exact.excitation_energy,
                                                    abs=1e-5)


@pytest.mark.parametrize("kwargs", [
    dict(omega_guess=[]), dict(omega_guess=[0.5, np.nan]),
    dict(omega_guess=0.5, explicit_width=-0.1),
    dict(omega_guess=0.5, explicit_width=np.inf),
])
def test_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        RelinearizedAdcMatrix(adc2_matrix("h2o_sto3g"), **kwargs)


@pytest.mark.parametrize("method", ["adc2x", "adc3"])
def test_non_trivial_complementary_block_not_implemented(method):
    refstate = testdata_cache.refstate("h2o_sto3g", case="gen")
    matrix = adcc.AdcMatrix(method, adcc.LazyMp(refstate))
    with pytest.raises(NotImplementedError):
        RelinearizedAdcMatrix(matrix, 0.5)
