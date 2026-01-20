#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2020 by the adcc authors
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


def __getattr__(attr):
    """
    Return a multi-dimensional block string like 'o1o2v1v1'
    when requesting the attribute 'ocvv'.
    """
    if any(c not in ["o", "v", "c"] for c in attr):
        raise AttributeError
    mapping = {"o": "o1", "v": "v1", "c": "o2"}
    return "".join(mapping[c] for c in attr)


def get_block_name(spaces, order, variant, has_core_occupied_space):
    """
    Assembles name of block function with variant, spaces, and order.
    Does some sanity checks.
    """
    if isinstance(variant, str):
        variant = [variant]
    elif variant is None:
        variant = []

    if has_core_occupied_space and "cvs" not in variant:
        raise ValueError("Cannot run a general (non-core-valence approximated) "
                         "ADC method on top of a ground state with a "
                         "core-valence separation.")
    if not has_core_occupied_space and "cvs" in variant:
        raise ValueError("Cannot run a core-valence approximated ADC method on "
                         "top of a ground state without a "
                         "core-valence separation.")

    return "_".join(["block"] + variant + spaces + [str(order)])
