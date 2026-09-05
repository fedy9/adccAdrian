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

from .functions import evaluate
from .AmplitudeVector import AmplitudeVector


class ComplementaryBlockInverter:
    """
    Approximate application of (omega - A_QQ)^{-1} to a vector living in a
    single excitation space (e.g. "pphh"), where the self-coupling A_QQ of
    that space is split as A_QQ = D + V. D is the diagonal, 0th order
    (bare orbital energy) part of the block -- diagonal and spin-blind by
    construction, so splitting it off this way keeps the resulting
    approximation spin-pure. V is everything else (see
    adc_pp.matrix.block_pphh_pphh_1_v for the concrete example).

    The space is split into an arbitrary number of disjoint index windows
    (a boolean mask each), each treated via a truncated Neumann series in
    V of its own order:

        A^{-1} v ~= sum_{k=0}^{order} (D_shifted^{-1} V)^k D_shifted^{-1} v

    with D_shifted = omega - D. `order=0` reduces to a plain divide by
    D_shifted for that window. There is never any coupling introduced
    between different windows, even when several are treated
    approximately: V's output is re-masked back onto the originating
    window after every application.

    Parameters
    ----------
    space : str
        Name of the excitation space this inverter acts on (e.g. "pphh").
    diagonal_0 : AmplitudeVector
        The 0th order (bare) diagonal of the space_space block, as
        returned by e.g. adc_pp.matrix.diagonal_pphh_pphh_0. Only the
        `space` block of this AmplitudeVector is used.
    v_apply : callable
        AmplitudeVector -> AmplitudeVector, the "fluctuation" part V of
        the space_space block (i.e. the actual block used minus its own
        0th order part). Only invoked for windows with `order > 0`.
    windows : list of (numpy.ndarray of bool, int)
        List of (mask, order) pairs. `mask` selects the part of the
        `space` tensor treated at the given Neumann `order` (>= 0). The
        masks must be pairwise disjoint; they need not cover the full
        space (uncovered elements are left at zero).
    """

    def __init__(self, space, diagonal_0, v_apply, windows):
        for mask, order in windows:
            if order < 0:
                raise ValueError("The Neumann expansion order must be >= 0.")
        for i, (mask_i, _) in enumerate(windows):
            for mask_j, _ in windows[i + 1:]:
                if np.any(mask_i & mask_j):
                    raise ValueError("The windows passed to ComplementaryBlockInverter "
                                     "must be pairwise disjoint.")
        self.space = space
        self.diagonal_0 = diagonal_0
        self.v_apply = v_apply
        self.windows = windows

    def _masked(self, ampl, mask):
        """
        Return a new single-block (`space`-only) AmplitudeVector holding
        the `space` block of `ampl`, restricted to `mask`. Any other
        blocks `ampl` might carry are dropped.
        """
        tensor = getattr(ampl, self.space)
        arr = tensor.to_ndarray()
        out = tensor.zeros_like()
        out.set_from_ndarray(arr * mask, 1e-14)
        return AmplitudeVector(**{self.space: out})

    def apply(self, numerator, omega):
        """
        Compute the approximate action of (omega - A_QQ)^{-1} on `numerator`
        (an AmplitudeVector holding only the `space` block).
        """
        d_shifted = evaluate(omega - self.diagonal_0)

        tensor = getattr(numerator, self.space)
        res = AmplitudeVector(**{self.space: tensor.zeros_like()})
        for mask, order in self.windows:
            term = evaluate(self._masked(numerator, mask) / d_shifted)
            window_res = term.copy()
            for _ in range(order):
                vterm = self._masked(self.v_apply(term), mask)
                term = evaluate(vterm / d_shifted)
                window_res += term
            res += window_res
        return res


__all__ = ["ComplementaryBlockInverter"]
