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
from .functions import evaluate
from .AmplitudeVector import AmplitudeVector


class ComplementaryBlockInverter:
    """
    Approximate application of `(omega - A_QQ)^{-1}` to a vector of a
    single excitation space (e.g. "pphh"), with `A_QQ = D + V` split into
    the spin-blind 0th order diagonal `D` and the remainder `V`. The space
    is split into disjoint windows, each treated by a truncated Neumann
    series of its own order,

        A^{-1} v ~= sum_{k=0}^{order} (D_shifted^{-1} V)^k D_shifted^{-1} v,

    with `D_shifted = omega - D` (order 0: a plain division). Windows
    never couple to each other: the output of `V` is masked back onto the
    originating window.

    Parameters
    ----------
    space : str
        Name of the excitation space (e.g. "pphh").
    diagonal_0 : AmplitudeVector
        The 0th order diagonal of the space_space block; only its `space`
        block is used.
    v_apply : callable
        AmplitudeVector -> AmplitudeVector, the remainder `V`. Only called
        for windows with `order > 0`.
    windows : list of (Tensor, int)
        Pairwise disjoint `(mask, order)` pairs. `mask` is a 0.0/1.0
        `Tensor` with the symmetry of the diagonal (not an amplitude
        tensor), e.g. from `_mask_from_selection` in AdcMatrix.py. The
        masks need not cover the whole space; uncovered elements give
        zero.
    """

    def __init__(self, space, diagonal_0, v_apply, windows):
        for mask, order in windows:
            if order < 0:
                raise ValueError("The Neumann expansion order must be >= 0.")
        for i, (mask_i, _) in enumerate(windows):
            for mask_j, _ in windows[i + 1:]:
                overlap = float(mask_i.dot(mask_j))
                if abs(overlap) > 1e-8:
                    raise ValueError("The windows passed to ComplementaryBlockInverter "
                                     "must be pairwise disjoint.")
        self.space = space
        self.diagonal_0 = diagonal_0
        self.v_apply = v_apply
        self.windows = [(evaluate(mask), order) for mask, order in windows]

    def _masked(self, ampl, iwindow):
        """`space` block of `ampl` restricted to window `iwindow`, as a
        single-block AmplitudeVector."""
        tensor = getattr(ampl, self.space)
        mask, _ = self.windows[iwindow]
        out = evaluate(tensor * mask)
        return AmplitudeVector(**{self.space: out})

    def apply(self, numerator, omega):
        """Approximate `(omega - A_QQ)^{-1} numerator` for a `space`-only
        AmplitudeVector `numerator`."""
        d_shifted = evaluate(omega - self.diagonal_0)

        tensor = getattr(numerator, self.space)
        res = AmplitudeVector(**{self.space: tensor.zeros_like()})
        for iwindow, (mask, order) in enumerate(self.windows):
            term = evaluate(self._masked(numerator, iwindow) / d_shifted)
            window_res = term.copy()
            for _ in range(order):
                vterm = self._masked(self.v_apply(term), iwindow)
                term = evaluate(vterm / d_shifted)
                window_res += term
            res += window_res
        return res


__all__ = ["ComplementaryBlockInverter"]
