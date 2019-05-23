#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the adcc authors
##
## This file is part of adcc.
##
## adcc is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## adcc is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with adcc. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------
import numpy as np
import scipy.constants

from adcc.timings import Timer


class SolverStateBase:
    def __init__(self, matrix):
        self.method = matrix.method       # The adc method which was used to
        #                                   obtain these states
        self.ground_state = matrix.ground_state  # The (MP) ground state upon
        #                                   which this solver state is based
        self.reference_state = matrix.ground_state.reference_state
        self.matrix = matrix
        self.eigenvalues = None           # Current eigenvalues
        self.eigenvectors = None          # Current eigenvectors
        self.converged = False            # Flag whether iteration is converged
        self.n_iter = 0                   # Number of iterations
        self.n_applies = 0                # Number of applies
        self.timer = Timer()              # Construct a new timer

    def describe(self):
        text = ""
        toeV = scipy.constants.value("Hartree energy in eV")

        if hasattr(self, "kind") and self.kind \
           and self.kind not in [" ", ""]:
            kind = self.kind + " "
        else:
            kind = ""

        if kind == "spin_flip" and hasattr(self, "spin_change") \
           and self.spin_change and self.spin_change != -1:
            spin_change = "(ΔMS={:+2d})".format(self.spin_change)
        else:
            spin_change = ""

        if self.converged:
            conv = "converged"
        else:
            conv = "NOT CONVERGED"

        text += "+" + 53 * "-" + "+\n"
        head = "| {0:15s}  {1:>34s} |\n"
        delim = ",  " if kind else ""
        text += head.format(self.method.name,
                            kind + spin_change + delim + conv)
        text += "+" + 53 * "-" + "+\n"

        # TODO Certain methods such as ADC(0), ADC(1) do not have
        #      a doubles part and it does not really make sense to
        #      display it here.

        # TODO Add dominant amplitudes
        body = "| {0:2d} {1:13.7g} {2:13.7g} {3:9.4g} {4:9.4g}  |\n"
        text += "|  #        excitation energy       |v1|^2    |v2|^2  |\n"
        text += "|          (au)           (eV)                        |\n"
        for i, vec in enumerate(self.eigenvectors):
            v1_norm = vec["s"].dot(vec["s"])
            if "d" in vec.blocks:
                v2_norm = vec["d"].dot(vec["d"])
            else:
                v2_norm = 0
            text += body.format(i, self.eigenvalues[i],
                                self.eigenvalues[i] * toeV,
                                v1_norm, v2_norm)
        text += "+" + 53 * "-" + "+"
        return text

    def _repr_pretty_(self, pp, cycle):
        if cycle:
            n_states = 0
            if self.eigenvalues is not None:
                n_states = len(self.eigenvalues)
            pp.text("SolverState(" + self.method.name + ", n_states=" +
                    str(n_states) + ")")
        else:
            pp.text(self.describe())
