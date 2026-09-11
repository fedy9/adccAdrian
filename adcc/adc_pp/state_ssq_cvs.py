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
"""
<S^2> for CVS-ADC excited states.

No CVS-specific formula exists for the 2-particle difference density matrix
(see `state_diffdm_2p.py`, which has no "cvs-*" dispatch entries). Instead of
deriving one, a CVS excitation vector is *embedded* into the excitation space
of a "flat" reference state built on the same orbitals but without
core-valence separation (`ReferenceState._flat_reference_state`): the
excluded (non-core) configurations are set to zero. Since a CVS state is,
by construction, exactly the state with zero amplitude on those
configurations, contracting this embedded vector with the ordinary
(non-CVS) ISR density-matrix formulas and the ordinary ssq_1p/ssq_2p
operators gives the same <S^2> as a (hypothetical) CVS-specific formula
would -- for any state which is dominantly described by single-core-hole
configurations. Cross-checks against independently-computed unrestricted
CVS-ADC(2)/(2)x results agree to ~1e-12 for such states, see
`test_state_ssq_restricted_matches_unrestricted_cvs`. Pure-doubles satellite
states, which are not spin-pure even in the non-CVS case, are the one known
exception.

Currently only CVS-ADC(2)-shaped excitation vectors (ph: 'o2v1',
pphh: 'o1o2v1v1', i.e. a single core hole) are supported.
"""
import numpy as np

from ..AmplitudeVector import AmplitudeVector
from ..Intermediates import Intermediates
from ..NParticleOperator import product_trace
from .state_diffdm import DISPATCH as DIFFDM_1P_DISPATCH
from .state_diffdm_2p import DISPATCH as DIFFDM_2P_DISPATCH


def _orbital_index_maps(mospaces_cvs, mospaces_flat):
    """
    Positions of the CVS valence-occupied (o1) and core-occupied (o2)
    orbitals within the merged occupied ("o1") space of the flat reference.
    """
    raw_o1 = list(mospaces_cvs.occupied_orbitals)
    raw_o2 = list(mospaces_cvs.core_orbitals)
    raw_o1_flat = list(mospaces_flat.occupied_orbitals)
    idx_o1_in_flat = [raw_o1_flat.index(r) for r in raw_o1]
    idx_o2_in_flat = [raw_o1_flat.index(r) for r in raw_o2]
    return idx_o1_in_flat, idx_o2_in_flat


def _embed_cvs_amplitude(excitation_vector, idx_o1_in_flat, idx_o2_in_flat,
                         template):
    """
    Embed a CVS-ADC(2) excitation vector into the flat excitation space
    defined by `template` (a zero AmplitudeVector on the flat reference),
    with zero coefficients on all excluded (non-core) configurations.
    """
    ev = excitation_vector
    if list(ev.ph.subspaces) != ["o2", "v1"] \
       or list(ev.pphh.subspaces) != ["o1", "o2", "v1", "v1"]:
        raise NotImplementedError(
            "<S^2> for CVS-ADC states is only implemented for CVS-ADC(2)"
            "-shaped excitation vectors (ph: o2 x v1, pphh: o1 x o2 x v1 x "
            f"v1); got ph subspaces {list(ev.ph.subspaces)} and pphh "
            f"subspaces {list(ev.pphh.subspaces)}."
        )
    n_v = template.ph.shape[1]
    ph_nd = ev.ph.to_ndarray()
    pphh_nd = ev.pphh.to_ndarray()

    ph_arr = np.zeros(template.ph.shape)
    ph_arr[np.ix_(idx_o2_in_flat, range(n_v))] = ph_nd

    # The flat pphh tensor is antisymmetric in its first two (merged
    # occupied) indices; the CVS pphh tensor lives entirely in the o1 x o2
    # block and has no such symmetry (o1 and o2 are different subspaces
    # there), so the antisymmetric partner block has to be filled in by hand.
    pphh_arr = np.zeros(template.pphh.shape)
    ix = np.ix_(idx_o1_in_flat, idx_o2_in_flat, range(n_v), range(n_v))
    pphh_arr[ix] = pphh_nd
    ix_swap = np.ix_(idx_o2_in_flat, idx_o1_in_flat, range(n_v), range(n_v))
    pphh_arr[ix_swap] = -pphh_nd.transpose(1, 0, 2, 3)

    ph_t = template.ph.copy()
    pphh_t = template.pphh.copy()
    ph_t.set_from_ndarray(ph_arr, 1e-14)
    pphh_t.set_from_ndarray(pphh_arr, 1e-14)
    return AmplitudeVector(ph=ph_t, pphh=pphh_t)


def cvs_state_ssq(ground_state, excitation_vector, property_method) -> float:
    """
    Compute <S^2> of a single CVS-ADC excited state.

    Parameters
    ----------
    ground_state : LazyMp
        The (CVS) ground state upon which the excitation was based.
    excitation_vector : AmplitudeVector
        The CVS excitation amplitude vector of the state.
    property_method : IsrMethod
        The (cvs-prefixed) ISR method used for property calculations,
        e.g. as returned by `ElectronicStates.property_method`.
    """
    flat_ground_state = ground_state._flat_ground_state
    flat_refstate = flat_ground_state.reference_state
    template = ground_state._flat_amplitude_template

    idx_o1_in_flat, idx_o2_in_flat = _orbital_index_maps(
        ground_state.mospaces, flat_ground_state.mospaces
    )
    amplitude = _embed_cvs_amplitude(
        excitation_vector, idx_o1_in_flat, idx_o2_in_flat, template
    )

    level_str = property_method.name.replace("cvs-", "")
    if level_str not in DIFFDM_2P_DISPATCH:
        raise NotImplementedError(
            f"<S^2> for CVS-ADC states is not implemented at the "
            f"'{property_method.name}' property level."
        )
    intermediates = Intermediates(flat_ground_state)
    ddm_1p = DIFFDM_1P_DISPATCH[level_str](flat_ground_state, amplitude,
                                           intermediates)
    ddm_2p = DIFFDM_2P_DISPATCH[level_str](flat_ground_state, amplitude,
                                           intermediates)

    ssq_1p_op = flat_refstate.operators.ssq_1p
    ssq_2p_op = flat_refstate.operators.ssq_2p
    gs_ssq = flat_ground_state.ssq(property_method.level.to_int())

    ssq_1p = product_trace(ssq_1p_op, ddm_1p)
    ssq_2p = product_trace(ssq_2p_op, ddm_2p)
    return ssq_1p + ssq_2p + gs_ssq
