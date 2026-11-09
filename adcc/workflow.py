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
import sys
import warnings

import numpy as np

from libadcc import ReferenceState

from . import solver
from .guess import (guesses_any, guesses_singlet, guesses_spin_flip,
                    guesses_triplet)
from .LazyMp import LazyMp
from .AdcMatrix import (AdcMatrix, AdcMatrixlike, AdcExtraTerm,
                        RelinearizedAdcMatrix)
from .AdcMethod import AdcMethod, IsrMethod
from .exceptions import InputError
from .functions import evaluate
from .ExcitedStates import ExcitedStates
from .ReferenceState import ReferenceState as adcc_ReferenceState
from .solver.lanczos import lanczos
from .solver.davidson import jacobi_davidson
from .solver.explicit_symmetrisation import (IndexSpinSymmetrisation,
                                             IndexSymmetrisation)

__all__ = ["run_adc"]


def run_adc(data_or_matrix, n_states=None, kind="any", conv_tol=None,
            eigensolver=None, guesses=None, n_guesses=None, isr_order=None,
            n_guesses_doubles=None, output=sys.stdout, core_orbitals=None,
            frozen_core=None, frozen_virtual=None, method=None,
            n_singlets=None, n_triplets=None, n_spin_flip=None,
            environment=None, relin=False, finalize_relinearized=True,
            guess_energies=None, explicit_width=None, order1_width=None,
            **solverargs):
    """Run an ADC calculation.

    Main entry point to run an ADC calculation. The reference to build the ADC
    calculation upon is supplied using the `data_or_matrix` argument.
    `adcc` is pretty flexible here. Possible options include:

        a. Hartree-Fock data from a host program, e.g. a molsturm SCF
           state, a pyscf SCF object or any class implementing the
           :py:class:`adcc.HartreeFockProvider` interface. From this data all
           objects mentioned in (b) to (d) will be implicitly created and will
           become available in the returned state.
        b. A :py:class:`adcc.ReferenceState` object
        c. A :py:class:`adcc.LazyMp` object
        d. A :py:class:`adcc.AdcMatrix` object

    Parameters
    ----------
    data_or_matrix
        Data containing the SCF reference
    n_states : int, optional
    kind : str, optional
    n_singlets : int, optional
    n_triplets : int, optional
    n_spin_flip : int, optional
        Specify the number and kind of states to be computed. Possible values
        for kind are "singlet", "triplet", "spin_flip" and "any", which is
        the default. For unrestricted references clamping spin-pure
        singlets/triplets is currently not possible and kind has to remain as
        "any". For restricted references `kind="singlets"` or `kind="triplets"`
        may be employed to enforce a particular excited states manifold.
        Specifying `n_singlets` is equivalent to setting `kind="singlet"` and
        `n_states=5`. Similarly for `n_triplets` and `n_spin_flip`.
        `n_spin_flip` is only valid for unrestricted references.

    conv_tol : float, optional
        Convergence tolerance to employ in the iterative solver for obtaining
        the ADC vectors (default: `1e-6` or 10 * SCF tolerance,
        whatever is larger)

    eigensolver : str, optional
        The eigensolver algorithm to use.

    n_guesses : int, optional
        Total number of guesses to compute. By default only guesses derived from
        the singles block of the ADC matrix are employed. See
        `n_guesses_doubles` for alternatives. If no number is given here
        `n_guesses = min(4, 2 * number of excited states to compute)`
        or a smaller number if the number of excitation is estimated to be less
        than the outcome of above formula.

    n_guesses_doubles : int, optional
        Number of guesses to derive from the doubles block. By default none
        unless n_guesses as explicitly given or automatically determined is
        larger than the number of singles guesses, which can be possibly found.

    guesses : list, optional
        Provide the guess vectors to be employed for the ADC run. Takes
        preference over `n_guesses` and `n_guesses_doubles`, such that these
        parameters are ignored. For a relinearized matrix, they are projected
        onto its explicit configurations and orthonormalized.

    isr_order: int or str, optional
        Order of the ISR expansion used for property calculations.
        Can be an integer (e.g. 3 for ISR(3)) or a string (e.g. "isr1s")

    output : stream, optional
        Python stream to which output will be written. If `None` all output
        is disabled.

    core_orbitals : int or list or tuple, optional
        The orbitals to be put into the core-occupied space. For ways to
        define the core orbitals see the description in
        :py:class:`adcc.ReferenceState`.
        Required if core-valence separation is applied and the input data is
        given as data from the host program (i.e. option (a) discussed above)

    frozen_core : int or list or tuple, optional
        The orbitals to select as frozen core orbitals (i.e. inactive occupied
        orbitals for both the MP and ADC methods performed). For ways to define
        these see the description in :py:class:`adcc.ReferenceState`.

    frozen_virtual : int or list or tuple, optional
        The orbitals to select as frozen virtual orbitals (i.e. inactive
        virtuals for both the MP and ADC methods performed). For ways to define
        these see the description in :py:class:`adcc.ReferenceState`.

    environment : bool or list or dict, optional
        The keywords to specify how coupling to an environment model,
        e.g. PE, is treated. For details see :ref:`environment`.

    relin : bool, optional
        If `True`, solve a :class:`adcc.RelinearizedAdcMatrix` built from
        the ADC matrix instead. Without `guess_energies`, `omega_guess` is
        estimated from the diagonal energies of the first `n_states`
        internal guess vectors. Default `False`. Has no effect if
        `data_or_matrix` already is a `RelinearizedAdcMatrix`.

    finalize_relinearized : bool, optional
        Only relevant when solving a :class:`adcc.RelinearizedAdcMatrix`
        (via `relin=True`, or because `data_or_matrix` already was one).

        - `True` (the default): each converged state is finalized by
          :meth:`RelinearizedAdcMatrix.finalize`: its energy becomes the
          self-consistent energy (Newton, no extra matrix-vector
          products) and its vector the full-space vector unfolded at that
          energy. The raw eigenvalue and the Newton iterations of each
          state are printed to `output`.
        - `False`: the raw solver eigenpairs are kept (energies at the
          fixed fold energy, vectors without the eliminated
          configurations).

        Has no effect for an ordinary (non-relinearized) matrix.

    guess_energies : float or list, optional
        Only with `relin=True`: target energy or energies used as
        `omega_guess` instead of the automatic estimate. Required if
        `guesses` are given.

    explicit_width, order1_width : float, optional
        Only with `relin=True`: forwarded to `RelinearizedAdcMatrix`.

    Other parameters
    ----------------
    max_subspace : int, optional
        Maximal subspace size
    max_iter : int, optional
        Maximal number of iterations
    max_subspace_iter : int, optional
        Maximum number of iterations for diagonalizing the subspace matrix
        during the Davidson solver procedure

    Returns
    -------
    ExcitedStates
        An :class:`adcc.ExcitedStates` object containing the
        :class:`adcc.AdcMatrix`, the :class:`adcc.LazyMp` ground state and the
        :class:`adcc.ReferenceState` as well as computed eigenpairs.

    Examples
    --------

    Run an ADC(2) calculation on top of a `pyscf` RHF reference of
    hydrogen flouride.

    >>> from pyscf import gto, scf
    ... mol = gto.mole.M(atom="H 0 0 0; F 0 0 1.1", basis="sto-3g")
    ... mf = scf.RHF(mol)
    ... mf.conv_tol_grad = 1e-8
    ... mf.kernel()
    ...
    ... state = adcc.run_adc(mf, method="adc2", n_singlets=3)

    The same thing can also be achieved using the `adcc.adcN` family of
    short-hands (see e.g. :py:func:`adcc.adc2`, :py:func:`adcc.cvs_adc2x`):

    >>> state = adcc.adc2(mf, n_singlets=3)

    Run a CVS-ADC(3) calculation of O2 with one core-occupied orbital

    >>> from pyscf import gto, scf
    ... mol = gto.mole.M(atom="O 0 0 0; O 0 0 1.2", basis="sto-3g")
    ... mf = scf.RHF(mol)
    ... mf.conv_tol_grad = 1e-8
    ... mf.kernel()
    ...
    ... state = adcc.cvs_adc3(mf, core_orbitals=1, n_singlets=3)
    """
    was_already_relin = isinstance(data_or_matrix, RelinearizedAdcMatrix)
    if relin and not was_already_relin and guesses is not None \
       and guess_energies is None:
        raise InputError("With relin=True, guess_energies must be given "
                         "together with guesses.")
    matrix = construct_adcmatrix(
        data_or_matrix, core_orbitals=core_orbitals, frozen_core=frozen_core,
        frozen_virtual=frozen_virtual, method=method, relin=relin, 
        omega_guess=guess_energies, explicit_width=explicit_width,
        order1_width=order1_width)

    n_states, kind = validate_state_parameters(
        matrix.reference_state, n_states=n_states, n_singlets=n_singlets,
        n_triplets=n_triplets, n_spin_flip=n_spin_flip, kind=kind)

    # Determine spin change during excitation. If guesses is not None,
    # i.e. user-provided, we cannot guarantee for obtaining a particular
    # spin_change in case of a spin_flip calculation.
    spin_change = None
    if kind == "spin_flip" and guesses is None:
        spin_change = -1

    # Select solver to run
    if eigensolver is None:
        eigensolver = "davidson"

    # Setup environment coupling terms and energy corrections. Safe to do
    # before the omega_guess update below: update_omega_guess() only
    # touches the tier partition, never `self.blocks`, so an
    # environment-coupling term added here is preserved regardless of
    # ordering.
    ret = setup_environment(matrix, environment)
    env_matrix_term, env_energy_corrections = ret
    # add terms to matrix
    if env_matrix_term:
        matrix += env_matrix_term

    if isinstance(matrix, RelinearizedAdcMatrix) and not was_already_relin:
        # Built with a placeholder omega_guess (see construct_adcmatrix);
        # a RelinearizedAdcMatrix passed in by the caller is left as is.
        if guess_energies is None:
            omega_guess = estimate_omega_guess(matrix, n_states, kind)
            if output is not None:
                print("No guess energies provided, taking the average of "
                      "internal guesses.", file=output)
            matrix.update_omega_guess(omega_guess)
    if isinstance(matrix, RelinearizedAdcMatrix) and output is not None:
        occupancy = matrix.window_occupancy
        print(f"RelinearizedAdcMatrix(omega_fixed={matrix.omega_fixed:.4f}, "
              f"{matrix.complementary_space}): " + ", ".join(
                  f"{100 * w['fraction']:.2f}% {name}"
                  for name, w in occupancy.items() if name != "total"
              ), file=output)

    if isinstance(matrix, RelinearizedAdcMatrix) and guesses is not None:
        guesses = _project_guesses_to_explicit(matrix, guesses, n_states)

    property_method = None
    if isr_order is not None:
        property_method = matrix.method.as_method_at_level(IsrMethod, isr_order)

    diagres = diagonalise_adcmatrix(
        matrix, n_states, kind, guesses=guesses, n_guesses=n_guesses,
        n_guesses_doubles=n_guesses_doubles, conv_tol=conv_tol, output=output,
        eigensolver=eigensolver, **solverargs)

    _finalize_relinearized_states(matrix, diagres, finalize_relinearized,
                                  output=output)

    exstates = ExcitedStates(diagres, property_method=property_method)
    exstates.kind = kind
    exstates.spin_change = spin_change

    # add environment corrections to excited states
    exstates += env_energy_corrections
    return exstates


#
# Individual steps
#
def construct_adcmatrix(data_or_matrix, core_orbitals=None, frozen_core=None,
                        frozen_virtual=None, method=None, relin=False,
                        omega_guess=None, explicit_width=None, order1_width=None):
    """
    Use the provided data or AdcMatrix object to check consistency of the
    other passed parameters and construct the AdcMatrix object representing
    the problem to be solved. If `relin`, and the result is not already a
    RelinearizedAdcMatrix, wrap it in one, with a placeholder
    `omega_guess=0.0` that run_adc updates once the target energies are
    known. Internal function called from run_adc.
    """
    if not isinstance(data_or_matrix, AdcMatrixlike) and method is None:
        raise InputError("method needs to be explicitly provided unless "
                         "data_or_matrix is an AdcMatrixlike.")
    if method is not None and not isinstance(method, AdcMethod):
        try:
            method = AdcMethod(method)
        except ValueError as e:
            raise InputError(str(e))  # In case the method is unknown

    if not isinstance(data_or_matrix, (ReferenceState, AdcMatrixlike, LazyMp)):
        if method.is_core_valence_separated and core_orbitals is None:
            raise InputError("If core-valence separation approximation is "
                             "applied then the number of core orbitals needs "
                             "to be specified via the parameter "
                             "core_orbitals.")
        try:
            refstate = adcc_ReferenceState(data_or_matrix,
                                           core_orbitals=core_orbitals,
                                           frozen_core=frozen_core,
                                           frozen_virtual=frozen_virtual)
        except ValueError as e:
            raise InputError(str(e))  # In case of an issue with the spaces
        data_or_matrix = refstate
    elif core_orbitals is not None:
        mospaces = data_or_matrix.mospaces
        warnings.warn("Ignored core_orbitals parameter because data_or_matrix"
                      " is a ReferenceState, a LazyMp or an AdcMatrixlike object "
                      " (which has a value of core_orbitals={})."
                      "".format(mospaces.n_orbs_alpha("o2")))
    elif frozen_core is not None:
        mospaces = data_or_matrix.mospaces
        warnings.warn("Ignored frozen_core parameter because data_or_matrix"
                      " is a ReferenceState, a LazyMp or an AdcMatrixlike object "
                      " (which has a value of frozen_core={})."
                      "".format(mospaces.n_orbs_alpha("o3")))
    elif frozen_virtual is not None:
        mospaces = data_or_matrix.mospaces
        warnings.warn("Ignored frozen_virtual parameter because data_or_matrix"
                      " is a ReferenceState, a LazyMp or an AdcMatrixlike object "
                      " (which has a value of frozen_virtual={})."
                      "".format(mospaces.n_orbs_alpha("v2")))

    # Make AdcMatrix (if not done)
    if isinstance(data_or_matrix, (ReferenceState, LazyMp)):
        try:
            matrix = AdcMatrix(method, data_or_matrix)
        except ValueError as e:
            # In case of an issue with CVS <-> chosen spaces
            raise InputError(str(e))
        if relin:
            if omega_guess is None:
                omega_guess = 0.0
            matrix = RelinearizedAdcMatrix(
                matrix, omega_guess=omega_guess, explicit_width=explicit_width,
                order1_width=order1_width)
        return matrix
    elif method is not None and method != data_or_matrix.method:
        warnings.warn("Ignored method parameter because data_or_matrix is an"
                      " AdcMatrixlike, which implicitly sets the method")
    if isinstance(data_or_matrix, AdcMatrixlike):
        if relin and not isinstance(data_or_matrix, RelinearizedAdcMatrix):
            if omega_guess is None:
                omega_guess = 0.0
            return RelinearizedAdcMatrix(
                data_or_matrix, omega_guess=omega_guess,
                explicit_width=explicit_width, order1_width=order1_width)
        return data_or_matrix


def validate_state_parameters(reference_state, n_states=None, n_singlets=None,
                              n_triplets=None, n_spin_flip=None, kind="any"):
    """
    Check the passed state parameters for consistency with itself and with
    the passed reference and normalise them. In the end return the number of
    states and the corresponding kind parameter selected.
    Internal function called from run_adc.
    """
    if sum(nst is not None for nst in [n_states, n_singlets,
                                       n_triplets, n_spin_flip]) > 1:
        raise InputError("One may only specify one out of n_states, "
                         "n_singlets, n_triplets and n_spin_flip")

    if n_singlets is not None:
        if not reference_state.restricted:
            raise InputError("The n_singlets parameter may only be employed "
                             "for restricted references")
        if kind not in ["singlet", "any"]:
            raise InputError(f"Kind parameter {kind} not compatible "
                             "with n_singlets > 0")
        kind = "singlet"
        n_states = n_singlets
    if n_triplets is not None:
        if not reference_state.restricted:
            raise InputError("The n_triplets parameter may only be employed "
                             "for restricted references")
        if kind not in ["triplet", "any"]:
            raise InputError(f"Kind parameter {kind} not compatible "
                             "with n_triplets > 0")
        kind = "triplet"
        n_states = n_triplets
    if n_spin_flip is not None:
        if reference_state.restricted:
            raise InputError("The n_spin_flip parameter may only be employed "
                             "for unrestricted references")
        if kind not in ["spin_flip", "any"]:
            raise InputError(f"Kind parameter {kind} not compatible "
                             "with n_spin_flip > 0")
        kind = "spin_flip"
        n_states = n_spin_flip

    # Check if there are states to be computed
    if n_states is None or n_states == 0:
        raise InputError("No excited states to be computed. Specify at least "
                         "one of n_states, n_singlets, n_triplets, "
                         "or n_spin_flip")
    if n_states < 0:
        raise InputError("n_states needs to be positive")

    if kind not in ["any", "spin_flip", "singlet", "triplet"]:
        raise InputError("The kind parameter may only take the values 'any', "
                         "'singlet', 'triplet' or 'spin_flip'")
    if kind in ["singlet", "triplet"] and not reference_state.restricted:
        raise InputError("kind==singlet and kind==triplet are only valid for "
                         "ADC calculations in combination with a restricted "
                         "ground state.")
    if kind in ["spin_flip"] and reference_state.restricted:
        raise InputError("kind==spin_flip is only valid for "
                         "ADC calculations in combination with an unrestricted "
                         "ground state.")
    return n_states, kind


def diagonalise_adcmatrix(matrix, n_states, kind, eigensolver="davidson",
                          guesses=None, n_guesses=None, n_guesses_doubles=None,
                          conv_tol=None, output=sys.stdout, **solverargs):
    """
    This function seeks appropriate guesses and afterwards proceeds to
    diagonalise the ADC matrix using the specified eigensolver.
    Internal function called from run_adc.
    """
    reference_state = matrix.reference_state

    # Determine default ADC convergence tolerance
    if conv_tol is None:
        conv_tol = max(10 * reference_state.conv_tol, 1e-6)
    if reference_state.conv_tol > conv_tol:
        raise InputError(
            "Convergence tolerance of SCF results "
            f"(== {reference_state.conv_tol}) needs to be lower than ADC "
            f"convergence tolerance parameter conv_tol (== {conv_tol})."
        )

    # Determine explicit_symmetrisation
    explicit_symmetrisation = IndexSymmetrisation
    if kind in ["singlet", "triplet"]:
        explicit_symmetrisation = IndexSpinSymmetrisation(
            matrix, enforce_spin_kind=kind
        )

    # Set some solver-specific parameters
    if eigensolver == "davidson":
        n_guesses_per_state = 2
        callback = setup_solver_printing(
            "Jacobi-Davidson", matrix, kind, solver.davidson.default_print,
            output=output)
        run_eigensolver = jacobi_davidson
    elif eigensolver == "lanczos":
        n_guesses_per_state = 1
        callback = setup_solver_printing(
            "Lanczos", matrix, kind, solver.lanczos.default_print,
            output=output)
        run_eigensolver = lanczos
    else:
        raise InputError(f"Solver {eigensolver} unknown, try 'davidson'.")

    # Obtain or check guesses
    if guesses is None:
        if n_guesses is None:
            # restrict to the number of available singles guesses if no doubles
            # are available
            n_guesses = estimate_n_guesses(
                matrix=matrix, n_states=n_states,
                singles_only=("pphh" not in matrix.axis_blocks),
                n_guesses_per_state=n_guesses_per_state
            )
        # A relinearized matrix offers only singles and explicit doubles
        # as guesses, possibly fewer than requested; n_states suffice.
        min_guesses = (n_states if isinstance(matrix, RelinearizedAdcMatrix)
                       else None)
        guesses = obtain_guesses_by_inspection(matrix, n_guesses, kind,
                                               n_guesses_doubles,
                                               min_guesses=min_guesses)
    else:
        if len(guesses) < n_states:
            raise InputError("Less guesses provided via guesses (== {}) "
                             "than states to be computed (== {})"
                             "".format(len(guesses), n_states))
        if n_guesses is not None:
            warnings.warn("Ignoring n_guesses parameter, since guesses are "
                          "explicitly provided.")
        if n_guesses_doubles is not None:
            warnings.warn("Ignoring n_guesses_doubles parameter, since guesses "
                          "are explicitly provided.")

    solverargs.setdefault("which", "SA")
    return run_eigensolver(matrix, guesses, n_ep=n_states, conv_tol=conv_tol,
                           callback=callback,
                           explicit_symmetrisation=explicit_symmetrisation,
                           **solverargs)


def estimate_omega_guess(matrix, n_states, kind):
    """
    Automatic `omega_guess` for `relin=True`: the diagonal energies of
    `n_states` guesses, obtained as the solver itself would, but from the
    actual diagonal of the RelinearizedAdcMatrix `matrix`, not from its
    relinearized one (which hides the eliminated configurations, at this
    point those of the placeholder `omega_guess`). Doubles fill up if
    there are too few singles guesses. Internal function called from
    run_adc.
    """
    actual = AdcMatrix(matrix.method, matrix.ground_state,
                       block_orders=matrix.block_orders,
                       intermediates=matrix.intermediates,
                       diagonal_precomputed=AdcMatrix.diagonal(matrix))
    guesses = obtain_guesses_by_inspection(actual, n_states, kind)
    diagonal = actual.diagonal()
    return np.array([g.dot(diagonal * g) for g in guesses])


def _project_guesses_to_explicit(matrix, guesses, n_states):
    """
    Restrict user `guesses` to the explicit configurations of the
    RelinearizedAdcMatrix `matrix` and orthonormalize them (Gram-Schmidt,
    two passes). `matvec` never populates the eliminated configurations,
    so a subspace free of them stays free of them; components there would
    otherwise give spurious zero eigenvalues. Guesses (almost) entirely
    in the eliminated space are dropped. Internal function called from
    run_adc.
    """
    projected = []
    for guess in guesses:
        v = matrix._active_masked(guess)
        for _ in range(2):
            for w in projected:
                v = v - w.dot(v) * w
        norm = np.sqrt(v.dot(v))
        if norm > 1e-6 * np.sqrt(guess.dot(guess)):
            projected.append(evaluate(v / norm))
    if len(projected) < n_states:
        raise InputError(f"Only {len(projected)} of the guesses remain "
                         f"linearly independent after projection onto the "
                         f"explicit configurations of the relinearized "
                         f"matrix, but {n_states} states are requested.")
    return projected


def _finalize_relinearized_states(matrix, diagres, finalize=True,
                                  output=sys.stdout):
    """
    If `matrix` is a RelinearizedAdcMatrix and `finalize` is set, replace
    each eigenpair in `diagres` in-place by its self-consistent energy and
    the full-space vector unfolded at it (see
    `RelinearizedAdcMatrix.finalize`), printing the Newton iterations of
    each state, and sort the states by their final energies. Internal
    function called from run_adc.
    """
    if not isinstance(matrix, RelinearizedAdcMatrix) or not finalize:
        return

    raw_energies = np.array(diagres.eigenvalues)
    if output is not None:
        print(f"Relinearized states: self-consistent energies by Newton "
              f"(omega_fixed = {matrix.omega_fixed:.10f}), vectors unfolded "
              f"to the full space", file=output)
    vectors, energies = [], []
    for i, (ampl, omega) in enumerate(zip(diagres.eigenvectors, raw_energies)):
        if output is not None:
            print(f"  State {i + 1:3d}   lambda = {omega:.10f}", file=output)
        energy, v, _ = matrix.finalize(ampl, omega, output=output)
        vectors.append(v)
        energies.append(energy)

    # The corrections differ per state, so close-lying states can swap.
    order = np.argsort(energies, kind="stable")
    if output is not None and np.any(order != np.arange(len(order))):
        print("  State order changed by the finalization; states are "
              "sorted by their final energies.", file=output)
    diagres.eigenvectors = [vectors[i] for i in order]
    diagres.eigenvalues = np.array(energies)[order]
    for attr in ("residual_norms", "residuals_converged"):
        values = getattr(diagres, attr, None)
        if values is not None:
            setattr(diagres, attr, np.asarray(values)[order])
    if getattr(diagres, "residuals", None) is not None:
        diagres.residuals = [diagres.residuals[i] for i in order]


def estimate_n_guesses(matrix, n_states, singles_only=True,
                       n_guesses_per_state=2):
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
        n_guesses = min(n_guesses, n_occ_a * n_virt_a)

    # Adjust if we overshoot the maximal number of sensible singles block
    # guesses, but make sure we get at least n_states guesses
    return max(n_states, n_guesses)


def obtain_guesses_by_inspection(matrix, n_guesses, kind, n_guesses_doubles=None,
                                 min_guesses=None):
    """
    Obtain guesses by inspecting the diagonal matrix elements.
    If n_guesses_doubles is not None, this is number is always adhered to.
    Otherwise the number of doubles guesses is adjusted to fill up whatever
    the singles guesses cannot provide to reach n_guesses.
    If min_guesses is given, fewer than n_guesses (but at least
    min_guesses) guesses are accepted.
    Internal function called from run_adc.
    """
    if n_guesses_doubles is not None and n_guesses_doubles > 0 \
       and "pphh" not in matrix.axis_blocks:
        raise InputError("n_guesses_doubles > 0 is only sensible if the ADC "
                         "method has a doubles block (i.e. it is *not* ADC(0), "
                         "ADC(1) or a variant thereof.")

    # Determine guess function
    guess_function = {"any": guesses_any, "singlet": guesses_singlet,
                      "triplet": guesses_triplet,
                      "spin_flip": guesses_spin_flip}[kind]

    # Determine number of singles guesses to request
    n_guess_singles = n_guesses
    if n_guesses_doubles is not None:
        n_guess_singles = n_guesses - n_guesses_doubles
    singles_guesses = guess_function(matrix, n_guess_singles, block="ph")

    doubles_guesses = []
    if "pphh" in matrix.axis_blocks:
        # Determine number of doubles guesses to request if not
        # explicitly specified
        if n_guesses_doubles is None:
            n_guesses_doubles = n_guesses - len(singles_guesses)
        if n_guesses_doubles > 0:
            doubles_guesses = guess_function(matrix, n_guesses_doubles,
                                             block="pphh")

    total_guesses = singles_guesses + doubles_guesses
    n_required = n_guesses if min_guesses is None else min_guesses
    if len(total_guesses) < n_required:
        hint = ""
        if isinstance(matrix, RelinearizedAdcMatrix):
            hint = (". Only singles and explicit doubles of the relinearized "
                    "matrix qualify; increase explicit_width")
        raise InputError("Less guesses found than requested: {} found, "
                         "{} requested{}".format(len(total_guesses), n_required,
                                                 hint))
    return total_guesses


def setup_solver_printing(solmethod_name, matrix, kind, default_print,
                          output=None):
    """
    Setup default printing for solvers. Internal function called from run_adc.
    """
    kstr = " "
    if kind != "any":
        kstr = " " + kind
    method_name = f"{matrix}"
    if hasattr(matrix, "method"):
        method_name = matrix.method.name

    if output is not None:
        print(f"Starting {method_name}{kstr} {solmethod_name} ...",
              file=output)

        def inner_callback(state, identifier):
            default_print(state, identifier, output)
        return inner_callback


def setup_environment(matrix, environment):
    """
    Setup environment matrix terms and/or energy corrections.
    Internal function called from run_adc.
    """
    valid_envs = ["ptss", "ptlr", "linear_response"]
    hf = matrix.reference_state
    if hf.environment and environment is None:
        raise InputError(
            "Environment found in reference state, but no environment"
            " configuration specified. Please select from the following"
            f" schemes: {valid_envs} or set to False."
        )
    elif environment and not hf.environment:
        raise InputError(
            "Environment specified, but no environment"
            " was found in reference state."
        )
    elif not hf.environment:
        environment = {}

    convertor = {
        bool: lambda value: {"ptss": True, "ptlr": True} if value else {},
        list: lambda value: {k: True for k in value},
        str: lambda value: {value: True},
        dict: lambda value: value,
    }
    conversion = convertor.get(type(environment), None)
    if conversion is None:
        raise TypeError("Cannot convert environment parameter of type"
                        f"'{type(environment)}' to dict.")
    environment = conversion(environment)

    if any(env not in valid_envs for env in environment):
        raise InputError("Invalid key specified for environment."
                         f" Valid keys are '{valid_envs}'.")

    env_matrix_term = None
    energy_corrections = []

    forbidden_combinations = [
        ["ptlr", "linear_response"],
    ]
    for fbc in forbidden_combinations:
        if all(environment.get(k, False) for k in fbc):
            raise InputError("Combination of environment schemes"
                             f" '{fbc}' not allowed. Check the"
                             " adcc documentation for more details.")

    for pt in ["ptss", "ptlr"]:
        if not environment.get(pt, False):
            continue
        hf_corr = hf.excitation_energy_corrections
        eec_key = f"{hf.environment}_{pt}_correction"
        if eec_key not in hf_corr:
            raise ValueError(f"{pt} correction requested, but could not find"
                             f" the needed function {eec_key} in"
                             f" reference state from backend {hf.backend}.")
        energy_corrections.append(hf_corr[eec_key])
    if environment.get("linear_response", False):
        from adcc.adc_pp import environment as adcpp_env
        block_key = f"block_ph_ph_0_{hf.environment}"
        if not hasattr(adcpp_env, block_key):
            raise NotImplementedError("Matrix term for linear response coupling"
                                      f" with solvent {hf.environment}"
                                      " not implemented.")
        block_fun = getattr(adcpp_env, block_key)
        env_matrix_term = AdcExtraTerm(matrix, {'ph_ph': block_fun})

    return env_matrix_term, energy_corrections
