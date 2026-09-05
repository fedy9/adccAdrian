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

    The space is split into two disjoint index windows (a boolean mask
    each); no coupling between the two windows is ever introduced, even
    when both are treated approximately:

    - `mask_expansion`: treated via a truncated Neumann series in V,
        A^{-1} v ~= sum_{k=0}^{order} (D_shifted^{-1} V)^k D_shifted^{-1} v
      with D_shifted = omega - D.
    - `mask_zeroth`: treated at 0th order only, i.e. a plain divide by
      D_shifted (equivalent to `order=0` for that part of the space).

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
        0th order part). Only invoked if `order > 0` and `mask_expansion`
        is not empty.
    mask_expansion : numpy.ndarray of bool
        Boolean mask (shape of the `space` tensor) selecting the part of
        the space treated via the order-`order` Neumann expansion.
    mask_zeroth : numpy.ndarray of bool
        Boolean mask selecting the part of the space treated at 0th
        order only. Must be disjoint from `mask_expansion`.
    order : int, optional
        Neumann expansion order used for `mask_expansion` (default 1).
        `order=0` makes the expansion part behave exactly like the
        zeroth-order part (no V ever applied) -- used automatically by
        RelinearizedAdcMatrix when the space_space block is 0th order to
        begin with, in which case V is identically zero anyway.
    """

    def __init__(self, space, diagonal_0, v_apply, mask_expansion,
                mask_zeroth, order=1):
        if order < 0:
            raise ValueError("The Neumann expansion order must be >= 0.")
        if np.any(mask_expansion & mask_zeroth):
            raise ValueError("mask_expansion and mask_zeroth must be disjoint.")
        self.space = space
        self.diagonal_0 = diagonal_0
        self.v_apply = v_apply
        self.mask_expansion = mask_expansion
        self.mask_zeroth = mask_zeroth
        self.order = order

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

        term = evaluate(self._masked(numerator, self.mask_expansion) / d_shifted)
        res = term.copy()
        for _ in range(self.order):
            vterm = self._masked(self.v_apply(term), self.mask_expansion)
            term = evaluate(vterm / d_shifted)
            res += term

        res += evaluate(self._masked(numerator, self.mask_zeroth) / d_shifted)
        return res


__all__ = ["ComplementaryBlockInverter"]
