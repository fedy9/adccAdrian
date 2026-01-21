#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the adcc authors
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

from .guesses_from_diagonal import guesses_from_diagonal
from adcc.exceptions import InputError
import warnings


def determine_spin_change(adc_type: str, kind: str, is_alpha: bool):
    if adc_type == "pp":
        if kind == "spin_flip":
            spin_change = -1
        else:
            spin_change = 0
    elif adc_type == "ip":
        if kind == "spin_flip":
            spin_change = 0
        else:
            spin_change = +0.5 - int(is_alpha)
    elif adc_type == "ea":
        if kind == "spin_flip":
            spin_change = 0
        else:
            spin_change = -0.5 + int(is_alpha)
    return spin_change


class Guesses:
    def __init__(self, guesses: list, kind: str,
                 adc_type: str, is_alpha: bool):
        if guesses is None:
            self.guesses = []
        else:
            self.guesses = guesses
        self.n_guesses = len(self.guesses)
        self.kind = kind
        self.is_alpha = is_alpha

        symmetrisation = {
            "singlet": "symmetric",
            "doublet": "none",
            "triplet": "antisymmetric",
            "spin_flip":"none",
            "any": "none"
        }
        self.spin_block_symmetrisation = symmetrisation.get(kind, 'none')

        if guesses is not None:
            self.spin_change = None
        else:
            self.spin_change = determine_spin_change(adc_type, kind, is_alpha)


    def estimate_n_guesses(self, matrix, n_states, n_guesses_per_state=2,
                           singles_only=True):
        """
        Implementation of a basic heuristic to find a good number of guess
        vectors to be searched for using the find_guesses function.
        Internal function called from run_adc.

        matrix             ADC matrix
        n_states           Number of states to be computed
        singles_only       Try to stay withing the singles excitation space
                        with the number of guess vectors.
        n_guesses_per_state  Number of guesses to search for for each state
        """
        # Try to use at least 4 or twice the number of states
        # to be computed as guesses
        n_guesses = n_guesses_per_state * max(2, n_states)

        if singles_only:
            # Compute the maximal number of sensible singles block guesses.
            # This is roughly the number of occupied alpha orbitals
            # times the number of virtual alpha orbitals
            #
            # If the system is core valence separated, then only the
            # core electrons count as "occupied".
            mospaces = matrix.mospaces
            sp_occ = "o2" if matrix.is_core_valence_separated else "o1"
            n_virt_a = mospaces.n_orbs_alpha("v1")
            n_occ_a = mospaces.n_orbs_alpha(sp_occ)
            estimate = n_occ_a * n_virt_a
            if matrix.method.level < 2 and matrix.method.adc_type != "pp":
                # Adjustment for IP- and EA-ADC(0/1) calculations
                estimate = n_occ_a if matrix.method.adc_type == "ip" else n_virt_a
            n_guesses = min(n_guesses, estimate)

        # Adjust if we overshoot the maximal number of sensible singles block
        # guesses, but make sure we get at least n_states guesses
        self.n_guesses = max(n_states, n_guesses)


    def obtain_guesses_by_inspection(self, matrix, n_guesses_doubles=None,
                                     **kwargs):
        """
        Obtain guesses by inspecting the diagonal matrix elements.
        If n_guesses_doubles is not None, this number is always adhered to.
        Otherwise the number of doubles guesses is adjusted to fill up whatever
        the singles guesses cannot provide to reach n_guesses.

        matrix      The matrix for which guesses are to be constructed
        is_alpha    Is the detached/attached electron alpha spin for the
                    respective IP-/EA-ADC calculation.
        kwargs      Any other argument understood by guesses_from_diagonal.
        """
        if n_guesses_doubles is not None and len(matrix.axis_blocks) < 2:
            raise InputError("n_guesses_doubles > 0 is only sensible if the "
                            "ADC method has a doubles block (i.e. it is *not* "
                            "ADC(0), ADC(1) or a variant thereof.")

        # Determine number of singles guesses to request
        if n_guesses_doubles is None:
            n_guesses_doubles = 0
        else:
            if len(matrix.axis_blocks) < 2:
                warnings.warn("Ignoring n_guesses_doubles parameter, since no "
                              "double guesses exist.")

        n_guesses_singles = self.n_guesses - n_guesses_doubles

        self.guesses += guesses_from_diagonal(
            matrix, n_guesses_singles, matrix.axis_blocks[0], self.kind,
            self.is_alpha, self.spin_change, self.spin_block_symmetrisation,
            **kwargs)

        # Determine number of doubles guesses to request if not
        # explicitly specified
        n_guesses_doubles = self.n_guesses - len(self.guesses)
        if n_guesses_doubles > 0:
            self.guesses += guesses_from_diagonal(
                matrix, n_guesses_doubles, matrix.axis_blocks[1], self.kind,
                self.is_alpha, self.spin_change,
                self.spin_block_symmetrisation, **kwargs)

        if len(self.guesses) < self.n_guesses:
            raise InputError("Less guesses found than requested: {} found, "
                            "{} requested".format(
                                len(self.guesses), self.n_guesses))
        self.n_guesses = len(self.guesses)