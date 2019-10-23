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
import os
import warnings

from .available_backends import available, first_available, have_backend

import h5py

__all__ = ["import_scf_results", "run_hf", "have_backend", "available"]


def import_scf_results(res):
    """
    Import an scf_result from an SCF program. Tries to be smart
    and guess what the host program was and how to import it.
    """
    if have_backend("pyscf"):
        from . import pyscf as backend_pyscf
        from pyscf import scf

        if isinstance(res, scf.hf.SCF):
            return backend_pyscf.import_scf(res)

    if have_backend("molsturm"):
        from . import molsturm as backend_molsturm
        from molsturm.State import State

        if isinstance(res, State):
            return backend_molsturm.import_scf(res)

    if have_backend("veloxchem"):
        import veloxchem as vlx

        from . import veloxchem as backend_veloxchem

        if isinstance(res, vlx.scfrestdriver.ScfRestrictedDriver):
            return backend_veloxchem.import_scf(res)

    if have_backend("psi4"):
        import psi4
        from . import psi4 as backend_psi4

        if isinstance(res, psi4.core.HF):
            return backend_psi4.import_scf(res)

    from libadcc import HartreeFockSolution_i
    if isinstance(res, HartreeFockSolution_i):
        return res

    if isinstance(res, (dict, h5py.File)):
        from adcc.DataHfProvider import DataHfProvider
        return DataHfProvider(res)

    if isinstance(res, str):
        if os.path.isfile(res) and (res.endswith(".h5")
                                    or res.endswith(".hdf5")):
            return import_scf_results(h5py.File(res, "r"))
        else:
            raise ValueError("Unrecognised path or file extension: {}"
                             "".format(res))

    # Note: Add more backends here

    raise NotImplementedError("No means to import an SCF result of "
                              "type " + str(type(res)) + " implemented.")


def run_hf(backend=None, xyz=None, basis="sto-3g", charge=0, multiplicity=1,
           conv_tol=1e-12, conv_tol_grad=1e-8, max_iter=150):
    """
        Run a HF calculation with a specified backend, molecule, and SCF
        parameters

        backend:        name of the backend (pyscf, psi4, or veloxchem)
        xyz:            string with coordinates in Bohr
        basis:          basis set name
        charge:         charge of the molecule
        multiplicity:   spin multiplicity 2S + 1
        conv_tol:       energy convergence tolerance
        conv_tol_grad:  convergence tolerance of the electronic gradient
        max_iter:       maximum number of SCF iterations
    """

    if not backend:
        backend = first_available()
        warnings.warn("No backend specified. Using {}.".format(backend))

    if not have_backend(backend):
        raise ValueError("Backend {} not found.".format(backend))
    if backend == "psi4":
        from . import psi4 as backend_hf
    elif backend == "pyscf":
        from . import pyscf as backend_hf
    elif backend == "veloxchem":
        from . import veloxchem as backend_hf
    elif backend == "molsturm":
        from . import molsturm as backend_hf
    else:
        raise NotImplementedError("No run_hf function implemented for backend "
                                  "{}.".format(backend))

    return backend_hf.run_hf(
        xyz, basis, charge, multiplicity, conv_tol,
        conv_tol_grad, max_iter
    )
