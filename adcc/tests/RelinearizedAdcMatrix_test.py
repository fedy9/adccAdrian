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
import numpy as np
import pytest

import adcc
from adcc.AdcMatrix import AdcMatrix, RelinearizedAdcMatrix

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


@pytest.mark.parametrize("method", ["adc2x", "adc3"])
def test_non_trivial_complementary_block_not_implemented(method):
    refstate = testdata_cache.refstate("h2o_sto3g", case="gen")
    matrix = adcc.AdcMatrix(method, adcc.LazyMp(refstate))
    with pytest.raises(NotImplementedError):
        RelinearizedAdcMatrix(matrix, 0.5)
