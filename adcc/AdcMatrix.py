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
import itertools
import warnings
import numpy as np

import libadcc

from .LazyMp import LazyMp
from .adc_pp import matrix as ppmatrix
from .timings import Timer, timed_member_call
from .AdcMethod import AdcMethod, Method, AdcType
from .functions import ones_like, evaluate
from .Intermediates import Intermediates
from .AmplitudeVector import AmplitudeVector
from .ComplementaryBlockInverter import ComplementaryBlockInverter


class AdcExtraTerm:
    def __init__(self, matrix, blocks):
        """Initialise an AdcExtraTerm.
        This class can be used to add customs terms
        to an existing :py:class:`AdcMatrix`

        Parameters
        ----------
        matrix : AdcMatrix
            The matrix for which the extra term
            should be created.
        blocks : dict
            A dictionary where the key labels the matrix block
            and the item denotes a callable to construct
            an :py:class:`AdcBlock`
        """
        self.ground_state = matrix.ground_state
        self.reference_state = matrix.reference_state
        self.intermediates = matrix.intermediates
        self.blocks = {}
        if not isinstance(blocks, dict):
            raise TypeError("blocks needs to be a dict.")
        for space in blocks:
            block_fun = blocks[space]
            if not callable(block_fun):
                raise TypeError("Items in additional_blocks must be callable.")
            block = block_fun(
                self.reference_state, self.ground_state, self.intermediates
            )
            self.blocks[space] = block


class AdcMatrixlike:
    """
    Base class marker for all objects like ADC matrices.
    """

    _special_block_orders = {
        "adc2x": {"ph_ph": 2, "ph_pphh": 1, "pphh_ph": 1, "pphh_pphh": 1},
        "isr1s": {"ph_ph": 1, "ph_pphh": None, "pphh_ph": None, "pphh_pphh": None},
        "isr2d": {"ph_ph": 2, "ph_pphh": 1, "pphh_ph": 1, "pphh_pphh": 0},
        "isr3d": {"ph_ph": 3, "ph_pphh": 2, "pphh_ph": 2, "pphh_pphh": 1},
    }

    @classmethod
    def _default_block_orders(cls, method: Method,
                              bandwidth: int) -> dict[str, int]:
        """
        Determines the default block orders for the given adc method.

        Parameters
        ----------
        method: Method
            The method to generate default block orders for.
        bandwidth: int
            The number of coupling blocks in the ADC/ISR matrix with
            non-vanishing zeroth-order contributions, e.g.,
            0 for the secular matrix and 1 for the 1-particle ISR matrix.
        """
        # check if we have a special method like adc2x
        # I guess base_method should also contain the adc_type prefix so
        # we don't need to separate different adc_types
        block_orders = cls._special_block_orders.get(
            method.base_method.name, None
        )
        if block_orders is not None:
            return block_orders.copy()
        # otherwise assume that we have a "normal" PP/IP/...-ADC(n) method
        # - determine which spaces are available in the ADC(n) matrix
        #   starting from the given minimal space
        min_space = {
            AdcType.PP: "ph"
        }.get(method.adc_type, None)
        if min_space is None:
            raise ValueError(f"Unknown adc type {method.adc_type.to_str()} for "
                             f"method {method.name}. Can not determine default "
                             "block orders.")
        assert bandwidth >= 0
        n_spaces = ((method.level.to_int() + bandwidth) // 2) + 1
        spaces = [
            "p" * i + min_space + "h" * i for i in range(0, n_spaces)
        ]
        # exploit the fact that the spaces are sorted from small to high:
        # If we walk the adc matrix in any direction we always have to subtract 1!
        # Therefore, we can determine the order according to the position of the
        # spaces: maxorder - 0 for singles, maxorder - 2 for doubles and
        # maxorder - 3 for the doubles/triples coupling
        ret = {}
        for ((i1, bra), (i2, ket)) in \
                itertools.product(enumerate(spaces), repeat=2):
            order = method.level.to_int() - i1 - i2
            # For ISR matrices allow missing diagonal blocks.
            order = None if order < 0 else order
            ret[f"{bra}_{ket}"] = order
        return ret

    @classmethod
    def _validate_block_orders(cls, block_orders: dict[str, int],
                               method: Method,
                               allow_missing_diagonal_blocks: bool = False) -> None:
        """
        Validates that the given block_orders form a valid adc/isr matrix for the
        given adc(isr) method.

        Parameters
        ----------
        block_orders: dict[str, int]
            The block orders to validate. Block orders should be of the form
            {'ph_ph': 2, 'ph_pphh': 1, ...}
        method: Method
            The adc/isr method/adc type (PP-ADC/ISR, ...) for which to validate
            the block_orders.
        allow_missing_diagonal_blocks: bool, optional
            If set, couplings between missing diagonal blocks are allowed, e.g.,
            {'ph_ph': 1, 'ph_pphh': 0, 'pphh_ph': 0}
            will be valid although there is a coupling between the 'ph_ph' block
            and the missing 'pphh_pphh' block. (default: False)
        """
        for block, order in block_orders.items():
            if order is None:
                continue
            assert order >= 0
            # ensure that the block is valid for the given adc type
            bra, ket = block.split("_")
            if not cls._is_valid_space(bra, method) or \
                    not cls._is_valid_space(ket, method):
                raise ValueError(f"Invalid block {block} for a "
                                 f"{method.adc_type.to_str()} ADC matrix.")
            if bra == ket:  # done for diagonal blocks
                continue
            # ensure that the matrix is symmetric
            inv_block = f"{ket}_{bra}"
            inv_order = block_orders.get(inv_block, None)
            if inv_order is None or inv_order != order:
                raise ValueError(f"{block} and {inv_block} should always have "
                                 "the same order.")
            if allow_missing_diagonal_blocks:
                continue
            # ensure that we have no coupling between missing diagonal blocks
            bra_diag = f"{bra}_{bra}"
            ket_diag = f"{ket}_{ket}"
            bra_diag_order = block_orders.get(bra_diag, None)
            ket_diag_order = block_orders.get(ket_diag, None)
            if bra_diag_order is None or ket_diag_order is None:
                raise ValueError(f"Can only have a couling block {block} "
                                 f"if both diagonal blocks {bra_diag} and "
                                 f"{ket_diag} are in the matrix too.")

    @classmethod
    def _is_valid_space(cls, space: str, method: Method) -> bool:
        """
        Checks whether the given space ('ph' for instance) is valid for the given
        adc method. Thereby we only verify that the space matches the adc_type of
        method!
        """
        n_particle, n_hole = space.count("p"), space.count("h")
        # ensure that the space is of the form pp...hh...
        if ("p" * n_particle + "h" * n_hole) != space:
            return False
        # depending on the adc type n_particle and n_hole have to
        # be equal or differ e.g. by +-1 (IP/EA)
        if method.adc_type is AdcType.PP:
            return n_particle == n_hole
        raise ValueError(f"Unknown adc type {method.adc_type.to_str()} for method "
                         f"{method.name}. Can not validate space.")


class AdcMatrix(AdcMatrixlike):

    def __init__(self, method, hf_or_mp, block_orders=None, intermediates=None,
                 diagonal_precomputed=None):
        """
        Initialise an ADC matrix.

        Parameters
        ----------
        method : str or AdcMethod
            Method to use.
        hf_or_mp : adcc.ReferenceState or adcc.LazyMp
            HF reference or MP ground state
        block_orders : optional
            The order of perturbation theory to employ for each matrix block.
            If not set, defaults according to the selected ADC method are chosen.
        intermediates : adcc.Intermediates or NoneType
            Allows to pass intermediates to re-use to this class.
        diagonal_precomputed: adcc.AmplitudeVector
            Allows to pass a pre-computed diagonal, for internal use only.
        """
        if isinstance(hf_or_mp, (libadcc.ReferenceState,
                                 libadcc.HartreeFockSolution_i)):
            hf_or_mp = LazyMp(hf_or_mp)
        if not isinstance(hf_or_mp, LazyMp):
            raise TypeError("hf_or_mp is not a valid object. It needs to be "
                            "either a LazyMp, a ReferenceState or a "
                            "HartreeFockSolution_i.")

        if not isinstance(method, AdcMethod):
            method = AdcMethod(method)

        if diagonal_precomputed:
            if not isinstance(diagonal_precomputed, AmplitudeVector):
                raise TypeError("diagonal_precomputed needs to be"
                                " an AmplitudeVector.")
            if diagonal_precomputed.needs_evaluation:
                raise ValueError("diagonal_precomputed must already"
                                 " be evaluated.")

        self.timer = Timer()
        self.method = method
        self.ground_state = hf_or_mp
        self.reference_state = hf_or_mp.reference_state
        self.mospaces = hf_or_mp.reference_state.mospaces
        self.is_core_valence_separated = method.is_core_valence_separated
        self.ndim = 2
        self.extra_terms = []

        self.intermediates = intermediates
        if self.intermediates is None:
            self.intermediates = Intermediates(self.ground_state)

        self.block_orders = self._default_block_orders(
            self.method, bandwidth=0
        )
        if block_orders is not None:
            self.block_orders.update(block_orders)
        self._validate_block_orders(
            block_orders=self.block_orders, method=self.method,
            allow_missing_diagonal_blocks=False
        )

        # Build the blocks and diagonals
        with self.timer.record("build"):
            variant = None
            if self.is_core_valence_separated:
                variant = "cvs"
            blocks = {
                block: ppmatrix.block(self.ground_state, block.split("_"),
                                      order=order, intermediates=self.intermediates,
                                      variant=variant)
                for block, order in self.block_orders.items() if order is not None
            }
            self.blocks = {bl: blocks[bl].apply for bl in blocks}
            if diagonal_precomputed:
                self._diagonal: AmplitudeVector = diagonal_precomputed
            else:
                self._diagonal: AmplitudeVector = sum(
                    bl.diagonal for bl in blocks.values() if bl.diagonal
                )
                self._diagonal.evaluate()
            self._init_space_data(self._diagonal)

    def __iadd__(self, other):
        """In-place addition of an :py:class:`AdcExtraTerm`

        Parameters
        ----------
        other : AdcExtraTerm
            the extra term to be added
        """
        if not isinstance(other, AdcExtraTerm):
            return NotImplemented
        if not all(k in self.blocks for k in other.blocks):
            raise ValueError("Can only add to blocks of"
                             " AdcMatrix that already exist.")
        for sp in other.blocks:
            orig_app = self.blocks[sp]
            other_app = other.blocks[sp].apply

            def patched_apply(ampl, original=orig_app, other=other_app):
                return sum(app(ampl) for app in (original, other))
            self.blocks[sp] = patched_apply
        other_diagonal = sum(bl.diagonal for bl in other.blocks.values()
                             if bl.diagonal)
        self._diagonal = self._diagonal + other_diagonal
        self._diagonal.evaluate()
        self.extra_terms.append(other)
        return self

    def __add__(self, other):
        """Addition of an :py:class:`AdcExtraTerm`, creating
        a copy of self and adding the term to the new matrix

        Parameters
        ----------
        other : AdcExtraTerm
            the extra term to be added

        Returns
        -------
        AdcMatrix
            a copy of the AdcMatrix with the extra term added
        """
        if not isinstance(other, AdcExtraTerm):
            return NotImplemented
        ret = AdcMatrix(self.method, self.ground_state,
                        block_orders=self.block_orders,
                        intermediates=self.intermediates,
                        diagonal_precomputed=self.diagonal())
        ret += other
        return ret

    def __radd__(self, other):
        return self.__add__(other)

    def _init_space_data(self, diagonal):
        """Update the cached data regarding the spaces of the ADC matrix"""
        self.axis_spaces = {}
        self.axis_lengths = {}
        for block in diagonal.blocks:
            self.axis_spaces[block] = getattr(diagonal, block).subspaces
            self.axis_lengths[block] = np.prod([
                self.mospaces.n_orbs(sp) for sp in self.axis_spaces[block]
            ])
        self.shape = (sum(self.axis_lengths.values()),
                      sum(self.axis_lengths.values()))

    def __repr__(self):
        ret = f"AdcMatrix({self.method.name}, "
        for b, o in self.block_orders.items():
            ret += f"{b}={o}, "
        return ret + ")"

    def __len__(self):
        return self.shape[0]

    @property
    def axis_blocks(self):
        """
        Return the blocks used along one of the axes of the ADC matrix
        (e.g. ['ph', 'pphh']).
        """
        # sort the keys by length to ensure that we always get
        # [singles, doubles, ...]
        return sorted(self.axis_spaces, key=len)

    def diagonal(self):
        """Return the diagonal of the ADC matrix"""
        return self._diagonal

    def block_apply(self, block, tensor):
        """
        Compute the application of a block of the ADC matrix
        with another AmplitudeVector or Tensor. Non-matching blocks
        in the AmplitudeVector will be ignored.
        """
        if not isinstance(tensor, libadcc.Tensor):
            raise TypeError("tensor should be an adcc.Tensor")

        with self.timer.record(f"apply/{block}"):
            outblock, inblock = block.split("_")
            ampl = AmplitudeVector(**{inblock: tensor})
            ret = self.blocks[block](ampl)
            return getattr(ret, outblock)

    @timed_member_call()
    def matvec(self, v):
        """
        Compute the matrix-vector product of the ADC matrix
        with an excitation amplitude and return the result.
        """
        return sum(block(v) for block in self.blocks.values())

    def rmatvec(self, v):
        # ADC matrix is symmetric
        return self.matvec(v)

    def __matmul__(self, other):
        if isinstance(other, AmplitudeVector):
            return self.matvec(other)
        if isinstance(other, list):
            if all(isinstance(elem, AmplitudeVector) for elem in other):
                return [self.matvec(ov) for ov in other]
        return NotImplemented

    def block_view(self, block):
        """
        Return a view into the AdcMatrix that represents a single
        block of the matrix. Currently only diagonal blocks are supported.
        """
        b1, b2 = block.split("_")
        if b1 != b2:
            raise NotImplementedError("Off-diagonal block views not yet "
                                      "implemented.")
            # TODO For off-diagonal blocks we probably need a different
            #      data structure as the AdcMatrix class as these block
            #      are inherently different than an AdcMatrix (no Hermiticity
            #      for example) and basically they only need to support some
            #      form of matrix-vector product and some statistics like
            #      spaces and sizes etc.
        block_orders = {bl: None for bl in self.block_orders.keys()}
        block_orders[block] = self.block_orders[block]
        return AdcMatrix(self.method, self.ground_state,
                         block_orders=block_orders,
                         intermediates=self.intermediates)

    def construct_symmetrisation_for_blocks(self):
        """
        Construct the symmetrisation functions, which need to be
        applied to relevant blocks of an AmplitudeVector in order
        to symmetrise it to the right symmetry in order to be used
        with the various matrix-vector-products of this function.

        Most importantly the returned functions antisymmetrise
        the occupied and virtual parts of the doubles parts
        if this is sensible for the method behind this adcmatrix.

        Returns a dictionary block identifier -> function
        """
        ret = {}
        if self.is_core_valence_separated:
            # CVS doubles part is antisymmetric wrt. (i,K,a,b) <-> (i,K,b,a)
            ret["pphh"] = lambda v: v.antisymmetrise([(2, 3)])
        else:
            def symmetrise_generic_adc_doubles(invec):
                # doubles part is antisymmetric wrt. (i,j,a,b) <-> (i,j,b,a)
                # doubles part is antisymmetric wrt. (i,j,a,b) <-> (j,i,a,b)
                scratch = invec.antisymmetrise([(0, 1)]).antisymmetrise([(2, 3)])
                # doubles part is symmetric wrt. (i,j,a,b) <-> (j,i,b,a)
                return scratch.symmetrise([(0, 1), (2, 3)])
            ret["pphh"] = symmetrise_generic_adc_doubles

            def symmetrise_generic_adc_triples(invec):
                # triples part is antisymmetric wrt. permutations of (i,j,k)
                # and wrt. permutations of (a,b,c)
                scratch = (
                    invec.antisymmetrise([(0, 1, 2)]).antisymmetrise([(3, 4, 5)])
                )
                # triples part is symmetric wrt. permutations of (i,j,k) and
                # permutations of (a,b,c)
                # NOTE: The doubles fix the (numerical) symmetry with a single
                # symmetrise call. This is not possible for triples, since the
                # following symmetrise call only covers 6 of the 18
                # even permutations generated by the 2 antisymmetrise calls above:
                # ijkabc + ikjacb + jikbac + jkibca + kijcab + kjicba
                return scratch.symmetrise([(0, 1, 2), (3, 4, 5)])
            ret["ppphhh"] = symmetrise_generic_adc_triples
        return ret

    def dense_basis(self, axis_blocks=None, ordering="adcc"):
        """
        Return the list of indices and their values
        of the dense basis representation

        ordering: adcc, spin, spatial
        """
        ret = []
        if axis_blocks is None:
            axis_blocks = self.axis_blocks
        if not isinstance(axis_blocks, list):
            axis_blocks = [axis_blocks]

        # Define function to impose the order in the basis
        if ordering == "adcc":
            def reduce_index(n_orbsa, idx):
                return idx, idx
        elif ordering == "spin":
            def reduce_index(n_orbsa, idx):
                is_beta = [idx[i] >= n_orbsa[i] for i in range(len(idx))]
                spatial = [idx[i] - n_orbsa[i] if is_beta[i] else idx[i]
                           for i in range(len(idx))]
                # Sort first by spin, then by spatial
                return (is_beta, spatial)
        elif ordering == "spatial":
            def reduce_index(n_orbsa, idx):
                is_beta = [idx[i] >= n_orbsa[i] for i in range(len(idx))]
                spatial = [idx[i] - n_orbsa[i] if is_beta[i] else idx[i]
                           for i in range(len(idx))]
                # Sort first by spatial, then by spin
                return (spatial, is_beta)

        if "ph" in axis_blocks:
            ret_s = []
            sp_s = self.axis_spaces["ph"]
            n_orbs_s = [self.mospaces.n_orbs(sp) for sp in sp_s]
            n_orbsa_s = [self.mospaces.n_orbs_alpha(sp) for sp in sp_s]
            for i in range(n_orbs_s[0]):
                for a in range(n_orbs_s[1]):
                    ret_s.append([((i, a), 1)])

            def sortfctn(x):
                return min(reduce_index(n_orbsa_s, idx) for idx, factor in x)
            ret_s.sort(key=sortfctn)
            ret_s.sort(key=sortfctn)
            ret.extend(ret_s)

        if "pphh" in axis_blocks:
            ret_d = []
            sp_d = self.axis_spaces["pphh"]
            n_orbsa_d = [self.mospaces.n_orbs_alpha(sp) for sp in sp_d]

            if sp_d[0] == sp_d[1] and sp_d[2] == sp_d[3]:
                nso = self.mospaces.n_orbs(sp_d[0])
                nsv = self.mospaces.n_orbs(sp_d[2])
                ret_d.extend([[((i, j, a, b), +1 / 2),
                               ((j, i, a, b), -1 / 2),
                               ((i, j, b, a), -1 / 2),
                               ((j, i, b, a), +1 / 2)]
                              for i in range(nso) for j in range(i)
                              for a in range(nsv) for b in range(a)])
            elif sp_d[2] == sp_d[3]:
                nso = self.mospaces.n_orbs(sp_d[0])
                nsc = self.mospaces.n_orbs(sp_d[1])
                nsv = self.mospaces.n_orbs(sp_d[2])
                ret_d.extend([[((i, j, a, b), +1 / np.sqrt(2)),
                               ((i, j, b, a), -1 / np.sqrt(2))]
                              for i in range(nso) for j in range(nsc)
                              for a in range(nsv) for b in range(a)])
            else:
                nso = self.mospaces.n_orbs(sp_d[0])
                nsc = self.mospaces.n_orbs(sp_d[1])
                nsv = self.mospaces.n_orbs(sp_d[2])
                nsw = self.mospaces.n_orbs(sp_d[3])
                ret_d.append([((i, j, b, a), 1)
                              for i in range(nso) for j in range(nsc)
                              for a in range(nsv) for b in range(nsw)])

            def sortfctn(x):
                return min(reduce_index(n_orbsa_d, idx) for idx, factor in x)
            ret_d.sort(key=sortfctn)
            ret_d.sort(key=sortfctn)
            ret.extend(ret_d)

        if any(b not in ("ph", "pphh") for b in self.axis_blocks):
            raise NotImplementedError("Blocks other than ph and pphh "
                                      "not implemented")
        return ret

    def to_ndarray(self, out=None):
        """
        Return the ADC matrix object as a dense numpy array. Converts the sparse
        internal representation of the ADC matrix to a dense matrix and return
        as a numpy array.

        Notes
        -----

        This method is only intended to be used for debugging and
        visualisation purposes as it involves computing a large amount of
        matrix-vector products and the returned array consumes a considerable
        amount of memory.

        The resulting matrix has no spin symmetry imposed, which means that
        its eigenspectrum may contain non-physical excitations (e.g. with linear
        combinations of α->β and α->α components in the excitation vector).

        This function has not been sufficiently tested to be considered stable.
        """
        # TODO Update to ph / pphh
        # TODO Still uses deprecated functions
        import tqdm

        from adcc import guess_zero

        # Get zero amplitude of the appropriate symmetry
        # (TODO: Only true for C1, where there is only a single irrep)
        ampl_zero = guess_zero(self)
        assert self.mospaces.point_group == "C1"

        # Build the shape of the returned array
        # Since the basis of the doubles block is not the unit vectors
        # this *not* equal to the shape of the AdcMatrix object
        basis = {b: self.dense_basis(b) for b in self.axis_blocks}
        mat_len = sum(len(basis[b]) for b in basis)

        if out is None:
            out = np.zeros((mat_len, mat_len))
        else:
            if out.shape != (mat_len, mat_len):
                raise ValueError("Output array has shape ({0:}, {1:}), but "
                                 "shape ({2:}, {2:}) is required."
                                 "".format(*out.shape, mat_len))
            out[:] = 0  # Zero all data in out.

        # Check for the cases actually implemented
        if any(b not in ("ph", "pphh") for b in self.axis_blocks):
            raise NotImplementedError("Blocks other than ph and pphh "
                                      "not implemented")
        if "ph" not in self.axis_blocks:
            raise NotImplementedError("Block 'ph' needs to be present")

        # Extract singles-singles block (contiguous)
        assert "ph" in self.axis_blocks
        n_orbs_ph = [self.mospaces.n_orbs(sp) for sp in self.axis_spaces["ph"]]
        n_ph = np.prod(n_orbs_ph)
        assert len(basis["ph"]) == n_ph
        view_ss = out[:n_ph, :n_ph].reshape(*n_orbs_ph, *n_orbs_ph)
        for i in range(n_orbs_ph[0]):
            for a in range(n_orbs_ph[1]):
                ampl = ampl_zero.copy()
                ampl.ph[i, a] = 1
                view_ss[:, :, i, a] = (self @ ampl).ph.to_ndarray()

        # Extract singles-doubles and doubles-doubles block
        if "pphh" in self.axis_blocks:
            assert self.axis_blocks == ["ph", "pphh"]
            view_sd = out[:n_ph, n_ph:].reshape(*n_orbs_ph, len(basis["pphh"]))
            view_dd = out[n_ph:, n_ph:]
            for j, bas1 in tqdm.tqdm(enumerate(basis["pphh"]),
                                     total=len(basis["pphh"])):
                ampl = ampl_zero.copy()
                for idx, val in bas1:
                    ampl.pphh[idx] = val
                ret_ampl = self @ ampl
                view_sd[:, :, j] = ret_ampl.ph.to_ndarray()

                for i, bas2 in enumerate(basis["pphh"]):
                    view_dd[i, j] = sum(val * ret_ampl.pphh[idx]
                                        for idx, val in bas2)

            out[n_ph:, :n_ph] = np.transpose(out[:n_ph, n_ph:])
        return out


class AdcMatrixShifted(AdcMatrix):
    def __init__(self, matrix, shift=0.0):
        """
        Initialise a shifted ADC matrix. Applying this class to a vector ``v``
        represents an efficient version of ``matrix @ v + shift * v``.

        Parameters
        ----------
        matrix : AdcMatrix
            Matrix which is shifted
        shift : float
            Value by which to shift the matrix
        """
        super().__init__(matrix.method, matrix.ground_state,
                         block_orders=matrix.block_orders,
                         intermediates=matrix.intermediates)
        self.shift = shift

    def matvec(self, in_ampl):
        out = super().matvec(in_ampl)
        out = out + self.shift * in_ampl
        return out

    def to_ndarray(self, out=None):
        super().to_ndarray(self, out)
        out = out + self.shift * np.eye(*out.shape)
        return out

    def block_apply(self, block, in_vec):
        ret = super().block_apply(block, in_vec)
        inblock, outblock = block.split("_")
        if inblock == outblock:
            ret += self.shift * in_vec
        return ret

    def diagonal(self):
        out = super().diagonal()
        out = out + self.shift  # Shift the diagonal
        return out

    def block_view(self, block):
        raise NotImplementedError("Block-view not yet implemented for "
                                  "shifted ADC matrices.")
        # TODO The way to implement this is to ask the inner matrix to
        #      a block_view and then wrap that in an AdcMatrixShifted.


class AdcMatrixProjected(AdcMatrix):
    def __init__(self, matrix, excitation_blocks, core_orbitals=None,
                 outer_virtuals=None):
        """
        Initialise a projected ADC matrix, i.e. represents the expression
        ``P @ M @ P`` where ``P`` is a projector onto a subset of
        ``excitation_blocks``.

        The ``excitation_blocks`` are defined by partitioning the ``o1`` occupied
        and ``v1`` virtual space of the ``matrix.mospaces`` into a core-occupied
        ``c``, valence-occupied ``o``, inner-virtual ``v`` and outer-virtual ``w``.
        This matrix will only keep selected blocks in the amplitudes non-zero, which
        are selected in the ``excitation_blocks`` list
        (e.g. ``["cv", "ccvv", "ocvv"]``).

        For details on the option how to select the spaces, see the documentation
        in :py:`adcc.ReferenceState.__init__` (``outer_virtuals`` follows the same
        rules as ``frozen_virtuals``).

        Parameters
        ----------
        matrix : AdcMatrix
            Matrix which is projected
        excitation_blocks : list
            Excitation blocks to keep in the Amplitudes.
        core_orbitals : int or list or tuple, optional
            The orbitals to be put into the ``c`` core space.
        outer_virtuals : int or list or tuple, optional
            The virtuals to be put into the ``w`` outer-virtual space.
        """
        from .projection import Projector, SubspacePartitioning

        for sp in excitation_blocks:
            if not any(len(sp) == len(ax_sp)
                       for ax_sp in matrix.axis_spaces.values()):
                raise ValueError(f"Invalid partition block {sp}.")

        super().__init__(matrix.method, matrix.ground_state,
                         block_orders=matrix.block_orders,
                         intermediates=matrix.intermediates)
        partitioning = SubspacePartitioning(matrix.mospaces, core_orbitals,
                                            outer_virtuals)

        projectors = {}
        for block in matrix.axis_spaces.keys():
            block_partitions = [sp for sp in excitation_blocks
                                if len(sp) == len(matrix.axis_spaces[block])]
            projectors[block] = Projector(matrix.axis_spaces[block],
                                          partitioning, block_partitions)
        self.projectors = projectors

    def apply_projection(self, in_ampl):
        return AmplitudeVector(**{
            block: self.projectors[block] @ in_ampl[block]
            for block in in_ampl.keys()
        })

    def matvec(self, in_ampl):
        in_proj = self.apply_projection(in_ampl)
        out = super().matvec(in_proj)
        return self.apply_projection(out)

    def block_apply(self, block, in_vec):
        inblock, outblock = block.split("_")
        in_proj = self.projectors[inblock].apply(in_vec)
        ret = super().block_apply(block, in_proj)
        return self.projectors[outblock].apply(ret)

    def diagonal(self):
        blocks = {}
        for (block, diagblock) in super().diagonal().items():
            # On the diagonal don't set the ignored amplitudes to zero,
            # but instead set them to a very large value to (a) avoid these being
            # selected for the initial guess and (b) ensure the preconditioning
            # naturally reduces components along this direction.
            P = self.projectors[block]
            one = ones_like(diagblock)
            blocks[block] = P @ diagblock + 100000 * (one - P @ one)
        return AmplitudeVector(**blocks)

    def block_view(self, block):
        raise NotImplementedError("Block-view not yet implemented for "
                                  "projected ADC matrices.")
        # TODO The way to implement this is to ask the inner matrix to
        #      a block_view and then wrap that in an AdcMatrixProjected.


def _select_within_cutoff_via_search(tensor, cutoff):
    """
    Return `(index_tuple, value)` for every canonical element of `tensor`
    with `|value| <= cutoff`, via a growing `Tensor.select_n_absmin`
    search (no densification). Sequential and not parallelised, so only
    the fallback of `_select_within_cutoff` (e.g. for CVS).
    """
    if cutoff <= 0 or tensor.size == 0:
        return []
    n = min(10, tensor.size)
    res = tensor.select_n_absmin(n)
    while (not res or max(abs(v) for _, v in res) <= cutoff) and n < tensor.size:
        n = max(n + 1, min(tensor.size, 2 * n))
        res = tensor.select_n_absmin(n)
    return [(tuple(idx), val) for idx, val in res if abs(val) <= cutoff]


def _select_within_cutoff_via_energies(tensor, reference_state, axis_spaces,
                                       omega_screen, cutoff):
    """
    Fast path for `_select_within_cutoff`: the bare diagonal is a sum of
    orbital energies (e.g. `D0_ijab = eps_a + eps_b - eps_i - eps_j`), so
    all canonical elements with `|D0 - omega_screen| <= cutoff` are found
    by vectorised NumPy over the occupied and virtual energy sums
    (`numpy.searchsorted` on the sorted virtual sums).

    Index tuples are enumerated *with* replacement, since repeated indices
    are distinct elements of the symmetric-declared diagonal tensor.
    `Tensor.is_allowed` is checked for every candidate. Returns `None` (use
    `_select_within_cutoff_via_search`) unless all occupied axes share one
    MO space and all virtual axes another (not the case e.g. for CVS).
    """
    occ_labels = [label for label in axis_spaces if label.startswith("o")]
    virt_labels = [label for label in axis_spaces if label.startswith("v")]
    if len(occ_labels) + len(virt_labels) != len(axis_spaces):
        return None
    if len(set(occ_labels)) != 1 or len(set(virt_labels)) != 1:
        return None
    r_occ, r_virt = len(occ_labels), len(virt_labels)

    eps_o = reference_state.orbital_energies(occ_labels[0]).to_ndarray()
    eps_v = reference_state.orbital_energies(virt_labels[0]).to_ndarray()
    if eps_o.size == 0 or eps_v.size == 0:
        return []

    occ_idx = np.array(list(
        itertools.combinations_with_replacement(range(eps_o.size), r_occ)
    ))
    virt_idx = np.array(list(
        itertools.combinations_with_replacement(range(eps_v.size), r_virt)
    ))
    occ_sums = eps_o[occ_idx].sum(axis=1)
    virt_sums = eps_v[virt_idx].sum(axis=1)
    order = np.argsort(virt_sums)
    virt_sums_sorted = virt_sums[order]

    lo = occ_sums + omega_screen - cutoff
    hi = occ_sums + omega_screen + cutoff
    left = np.searchsorted(virt_sums_sorted, lo, side="left")
    right = np.searchsorted(virt_sums_sorted, hi, side="right")

    selection = []
    for oi in range(occ_idx.shape[0]):
        occ_tuple = tuple(int(x) for x in occ_idx[oi])
        for pos in range(int(left[oi]), int(right[oi])):
            vi = order[pos]
            idx = occ_tuple + tuple(int(x) for x in virt_idx[vi])
            if not tensor.is_allowed(idx):
                continue
            selection.append((idx, float(virt_sums[vi] - occ_sums[oi])))
    return selection


def _select_within_cutoff(template, omega_screen, cutoff,
                          reference_state=None, axis_spaces=None):
    """
    Return `(index_tuple, value)` for every canonical element of the
    unshifted diagonal `template` with `|value - omega_screen| <= cutoff`,
    without densifying `template`. Uses
    `_select_within_cutoff_via_energies` if possible, else
    `_select_within_cutoff_via_search` on the shifted tensor.
    """
    if cutoff <= 0 or template.size == 0:
        return []
    if reference_state is not None and axis_spaces is not None:
        fast = _select_within_cutoff_via_energies(
            template, reference_state, axis_spaces, omega_screen, cutoff
        )
        if fast is not None:
            return fast
    shifted = evaluate(template - omega_screen)
    # shift back, so both backends return raw D0 values
    return [(idx, v + omega_screen)
           for idx, v in _select_within_cutoff_via_search(shifted, cutoff)]


def _mask_from_selection(template, selection):
    """
    Build a 0.0/1.0 `Tensor` with the shape and symmetry of `template`,
    set to 1.0 at the indices in `selection` (and, via `__setitem__`, at
    all symmetry-equivalent elements). `selection` can contain many
    symmetry-equivalent entries, so each is read first and only written
    if not yet set, since `__setitem__` is far more expensive than
    `__getitem__`.
    """
    mask = template.zeros_like()
    for idx, _ in selection:
        if mask[idx] > 0.5:
            continue
        mask[idx] = 1.0
    return evaluate(mask)


class RelinearizedAdcMatrix(AdcMatrix):
    """
    Relinearized ("folded") ADC matrix: the configurations of one
    excitation space (the "complementary" space, by default the highest
    one, e.g. "pphh" for ADC(2)) are split by their bare diagonal `D0`:

    - explicit: `D0 <= omega_screen + explicit_width`, kept as unknowns
      and coupled via the real, unapproximated complementary block;
    - eliminated: all others, folded into the remaining spaces via
      `(omega_fixed - D)^{-1}`, with `omega_fixed = mean(omega_guess)` and
      `omega_screen = max(omega_guess)`.

    The window is one-sided: configurations below the targeted states are
    always explicit, so every eliminated `D0` lies above them and the fold
    has no poles there. Since `omega_fixed` is fixed, the resulting
    operator is linear and can be diagonalised with the ordinary
    Davidson/Lanczos solvers; its eigenvalues are exact only for
    `omega_fixed` equal to the true eigenvalue. Use `finalize` after
    solving (as `run_adc` does by default) to get the self-consistent
    energy of each state and its full-space eigenvector.

    `matvec` never populates the eliminated configurations, so guess
    vectors must not have components there (they span a zero eigenspace).
    `run_adc` projects user guesses accordingly; when calling a solver
    directly, pass guesses restricted to the explicit configurations.

    Only a 0th order complementary block (e.g. the default space of
    ADC(2) and ADC(4)) is supported, for which `A_QQ = D` is diagonal and
    the elimination is exact at a given energy. Code paths for higher
    orders (an additional 1st order tier treated by a Neumann series in
    `V = A_QQ - D`, and `include_coupling`) are kept dormant for a future
    ADC(3) extension.

    Parameters
    ----------
    matrix : AdcMatrix
        The full ADC matrix to relinearize: a plain `AdcMatrix` without
        extra terms (add those to the relinearized matrix via `+=`).
    omega_guess : float or Sequence[float]
        Approximate energy/energies of the targeted state(s). The fold is
        evaluated at their mean (`omega_fixed`), the explicit window is
        measured from their maximum (`omega_screen`), so that it covers
        the highest targeted state.
    explicit_width : float, optional
        Width of the explicit window in Ha above `omega_screen`. Default
        `1.0`, a conservative choice: for ADC(2) the self-consistent
        energy is already accurate without any explicit configurations.
    order1_width : float, optional
        Dormant (higher orders only): width in Ha of the 1st order tier
        on top of `explicit_width`. Default `4.0`.
    complementary_space : str, optional
        Excitation space to relinearize. Defaults to the highest
        excitation class present in `matrix` (`matrix.axis_blocks[-1]`).
    include_coupling : bool, optional
        Dormant (higher orders only): whether the explicit tier couples
        directly to the eliminated tiers via the complementary block.
    """

    _DEFAULT_EXPLICIT_WIDTH = 1.0
    _DEFAULT_ORDER1_WIDTH = 4.0

    def __init__(self, matrix, omega_guess, explicit_width=None,
                order1_width=None, complementary_space=None,
                include_coupling=True):
        # The blocks are rebuilt from the method, which would silently drop
        # anything else the input matrix carries.
        if type(matrix) is not AdcMatrix or matrix.extra_terms:
            raise NotImplementedError(
                "RelinearizedAdcMatrix can only be built from a plain "
                "AdcMatrix without extra terms. Add extra terms to the "
                "RelinearizedAdcMatrix afterwards via +=."
            )
        super().__init__(matrix.method, matrix.ground_state,
                         block_orders=matrix.block_orders,
                         intermediates=matrix.intermediates)

        if complementary_space is None:
            complementary_space = self.axis_blocks[-1]
        if complementary_space not in self.axis_blocks:
            raise ValueError(f"Invalid complementary_space {complementary_space}: "
                             f"not one of {self.axis_blocks}.")
        if len(self.axis_blocks) < 2:
            raise ValueError("RelinearizedAdcMatrix needs at least one "
                             "excitation space besides the complementary one.")
        self.complementary_space = complementary_space
        self.include_coupling = include_coupling

        if explicit_width is None:
            explicit_width = self._DEFAULT_EXPLICIT_WIDTH
        if order1_width is None:
            order1_width = self._DEFAULT_ORDER1_WIDTH
        if explicit_width < 0:
            raise ValueError("Need explicit_width >= 0.")
        if order1_width <= 0:
            raise ValueError("Need order1_width > 0.")
        self.explicit_width = explicit_width
        self.order1_width = order1_width

        space_space = f"{complementary_space}_{complementary_space}"
        order_used = self.block_orders.get(space_space, None)
        if order_used is None:
            raise ValueError(f"The {space_space} block is not part of this "
                             "ADC matrix.")

        # Independent of omega_guess, so not redone by update_omega_guess
        variant = "cvs" if self.is_core_valence_separated else None
        diagonal_0_block = ppmatrix.block(
            self.ground_state, [complementary_space, complementary_space],
            order=0, intermediates=self.intermediates, variant=variant
        )
        self.diagonal_0 = diagonal_0_block.diagonal

        self.is_trivial = (order_used == 0)
        # The non-trivial (V != 0) code paths below and throughout this
        # class (order-1 tier, Neumann series, include_coupling) are kept
        # dormant for a future ADC(3) extension.
        if not self.is_trivial:
            raise NotImplementedError(
                f"RelinearizedAdcMatrix is only implemented for a 0th order "
                f"{space_space} block (e.g. the default complementary space "
                f"of ADC(2)), not for order {order_used} ({self.method.name})."
            )
        if self.is_trivial:
            self._v_apply = None
        elif (complementary_space == "pphh" and order_used == 1
             and not self.is_core_valence_separated):
            v_block = ppmatrix.block_pphh_pphh_1_v(
                self.reference_state, self.ground_state, self.intermediates
            )
            self._v_apply = v_block.apply
        else:
            # Generic (less efficient, but always correct) fallback:
            # V = (block at its actual order) - (block at 0th order)
            order_n_apply = self.blocks[space_space]
            order_0_apply = diagonal_0_block.apply

            def v_apply(ampl):
                return evaluate(order_n_apply(ampl) - order_0_apply(ampl))
            self._v_apply = v_apply

        self.update_omega_guess(omega_guess)

    def update_omega_guess(self, omega_guess):
        """
        (Re-)derive `omega_fixed`/`omega_screen` from `omega_guess` and
        recompute the partition (masks, `window_occupancy`, `inverter`).
        Everything else, e.g. environment terms added to `self.blocks`
        after construction, is kept. This allows building the matrix with
        a placeholder `omega_guess` before the target energies are known.
        """
        omega_guesses = np.atleast_1d(np.asarray(omega_guess, dtype=float))
        self.omega_fixed = float(np.mean(omega_guesses))
        self.omega_screen = float(np.max(omega_guesses))

        # One-sided windows: every configuration with D0 <= top is
        # selected. As D0 > 0 (aufbau reference), [0, top] is selected as
        # the symmetric window of radius top / 2 around top / 2. Only the
        # (small) selected tiers are searched; the eliminated 0th order
        # tier is the complement, so the diagonal is never densified.
        # For a 0th order block there is no 1st order tier (V = 0).
        space = self.complementary_space
        template = getattr(self.diagonal_0, space)
        axis_spaces = self.axis_spaces[space]
        explicit_top = self.omega_screen + self.explicit_width
        if self.is_trivial:
            active_sel = _select_within_cutoff(
                template, explicit_top / 2, explicit_top / 2,
                reference_state=self.reference_state, axis_spaces=axis_spaces
            )
            order1_sel = []
        else:
            order1_top = explicit_top + self.order1_width
            within_order1_cutoff = _select_within_cutoff(
                template, order1_top / 2, order1_top / 2,
                reference_state=self.reference_state, axis_spaces=axis_spaces
            )
            active_sel = [(idx, v) for idx, v in within_order1_cutoff
                          if v <= explicit_top]
            order1_sel = [(idx, v) for idx, v in within_order1_cutoff
                          if v > explicit_top]

        active_tensor = _mask_from_selection(template, active_sel)
        self._mask_active_tensor = active_tensor
        self._mask_inactive_tensor = evaluate(1 - active_tensor)

        elimination_windows = []
        if order1_sel:
            order1_tensor = _mask_from_selection(template, order1_sel)
            elimination_windows.append((order1_tensor, 1))
            order0_tensor = evaluate(1 - active_tensor - order1_tensor)
        else:
            order0_tensor = self._mask_inactive_tensor
        elimination_windows.append((order0_tensor, 0))

        n_total = template.size
        n_explicit = int(round(float(active_tensor.dot(active_tensor))))
        self.window_occupancy = {
            "total": n_total,
            "explicit": {"count": n_explicit, "fraction": n_explicit / n_total},
        }
        for mask, order in elimination_windows:
            n = int(round(float(mask.dot(mask))))
            self.window_occupancy[f"order-{order}"] = {
                "count": n, "fraction": n / n_total
            }

        self.inverter = ComplementaryBlockInverter(
            self.complementary_space, self.diagonal_0, self._v_apply,
            elimination_windows
        )
        self._lowest_eliminated_d0 = self._find_lowest_eliminated_d0()

    def _find_lowest_eliminated_d0(self):
        """Lowest bare diagonal `D0` of the eliminated configurations
        (`inf` if there are none), i.e. the lowest pole of the fold. Found
        by searching windows of growing width above the explicit window,
        which contains every lower configuration."""
        occupancy = self.window_occupancy
        if occupancy["explicit"]["count"] == occupancy["total"]:
            return np.inf
        space = self.complementary_space
        explicit_top = self.omega_screen + self.explicit_width
        width = 0.1
        while width < 1e4:
            found = [v for _, v in _select_within_cutoff(
                getattr(self.diagonal_0, space), explicit_top + width / 2,
                width / 2, reference_state=self.reference_state,
                axis_spaces=self.axis_spaces[space]
            ) if v > explicit_top]
            if found:
                return min(found)
            width *= 2
        return np.inf

    def _active_masked(self, v):
        """Copy of `v` with its `complementary_space` block restricted to
        the explicit configurations (eliminated amplitudes are always
        recomputed, never taken from an input vector)."""
        v_active = v.copy()
        space = self.complementary_space
        active_tensor = getattr(v_active, space)
        v_active[space] = evaluate(active_tensor * self._mask_active_tensor)
        return v_active

    def _build_numerator(self, v_active):
        """For `matvec`: the numerator `b` driving the eliminated
        configurations (as `_coupling_into_eliminated`) and the result of
        all blocks on the active-masked `v_active`, with its
        `complementary_space` block restricted to the explicit
        configurations."""
        space = self.complementary_space
        space_space = f"{space}_{space}"

        # All blocks except the complementary self-coupling
        result = v_active.zeros_like()
        for name, block_fn in self.blocks.items():
            if name != space_space:
                result += block_fn(v_active)
        external = evaluate(getattr(result, space))

        # Real complementary block on the explicit part: exact explicit
        # self-coupling, plus explicit->eliminated (if include_coupling)
        self_coupling = evaluate(
            getattr(self.blocks[space_space](v_active), space)
        )

        # Only the explicit configurations are unknowns
        result[space] = evaluate(
            self._mask_active_tensor * (external + self_coupling)
        )

        if self.include_coupling:
            numerator_tensor = evaluate(
                self._mask_inactive_tensor * (external + self_coupling)
            )
        else:
            numerator_tensor = evaluate(self._mask_inactive_tensor * external)
        numerator = AmplitudeVector(**{space: numerator_tensor})

        return numerator, result

    def __add__(self, other):
        # AdcMatrix.__add__ would return a plain AdcMatrix, silently
        # dropping the relinearization.
        if not isinstance(other, AdcExtraTerm):
            return NotImplemented
        raise NotImplementedError("Adding to a RelinearizedAdcMatrix is only "
                                  "supported in place, via +=.")

    def matvec(self, v):
        space = self.complementary_space
        space_space = f"{space}_{space}"

        v_active = self._active_masked(v)
        numerator, result = self._build_numerator(v_active)
        v_compl = self.inverter.apply(numerator, self.omega_fixed)

        # Fold the eliminated part back into the other excitation spaces.
        for name, block_fn in self.blocks.items():
            b1, b2 = name.split("_")
            if b2 == space and name != space_space:
                result += block_fn(v_compl)

        if self.include_coupling:
            fed_back = evaluate(
                getattr(self.blocks[space_space](v_compl), space)
            )
            result_space = getattr(result, space)
            result[space] = evaluate(
                result_space + fed_back * self._mask_active_tensor
            )

        return result

    def finalize(self, ampl, omega, output=None):
        """
        Finalize a converged relinearized eigenpair `(ampl, omega)`:
        determine its self-consistent energy (see
        `_self_consistent_energy`) and unfold the eliminated configurations
        at that energy (see `_unfold`), giving the full-space eigenvector.
        For a 0th order complementary block, the returned energy is exactly
        the Rayleigh quotient of the returned vector with the real,
        unapproximated matrix. No matrix-vector product with the full
        matrix is needed: besides elementwise operations on the
        complementary space, the only cost is one apply of the blocks
        coupling into it (see `_coupling_into_eliminated`).

        Parameters
        ----------
        ampl : AmplitudeVector
            A converged eigenvector of this matrix (e.g. an entry of
            `state.eigenvectors`). Only its explicit-window entries in
            `complementary_space` are used; anything outside the explicit
            window is discarded and recomputed.
        omega : float
            The converged (raw) eigenvalue belonging to `ampl`.
        output : file-like, optional
            If given, the Newton iterations are printed to it.

        Returns
        -------
        tuple(float, AmplitudeVector, int)
            The self-consistent energy, the full, unit-normalized
            eigenvector and the number of Newton steps taken. For a state
            not described by the fold (see `_self_consistent_energy`), the
            raw `omega` and `ampl` with zero steps, and a warning.
        """
        x_p = self._active_masked(ampl)
        b = self._coupling_into_eliminated(x_p)
        energy, x_q, n_steps = self._self_consistent_energy(x_p, b, omega,
                                                            output=output)
        if x_q is None:
            return omega, ampl, 0
        return energy, self._unfold(x_p, x_q), n_steps

    def _coupling_into_eliminated(self, v_active):
        """
        The numerator `b = M_QP x_P` driving the eliminated configurations,
        for the already active-masked `v_active`: only the blocks coupling
        into `complementary_space` are applied (the complementary block
        itself only if `include_coupling` and it is not 0th order -- a
        0th order block is diagonal, so it never couples explicit into
        eliminated configurations), and the result is restricted to the
        eliminated configurations. Equal to the numerator of
        `_build_numerator`, without everything else `matvec` needs.
        """
        space = self.complementary_space
        b = 0
        for name, block_fn in self.blocks.items():
            out_space, in_space = name.split("_")
            if out_space != space:
                continue
            if in_space == space and (self.is_trivial
                                      or not self.include_coupling):
                continue
            b = b + getattr(block_fn(v_active), space)
        return AmplitudeVector(**{
            space: evaluate(self._mask_inactive_tensor * b)
        })

    def _self_consistent_energy(self, v_active, numerator, omega, tol=1e-6,
                                maxiter=50, output=None):
        """
        Self-consistent energy of the converged relinearized eigenpair
        `(v_active, omega)`: the root of

            w = a + h(w),   h(w) = <b|(w - D)^{-1}|b> / <x_P|x_P>,

        i.e. of `w = <x_P|M(w)|x_P> / <x_P|x_P>` for the folded effective
        matrix `M(w)` (the self-consistency of Loewdin partitioning), with
        the explicit vector `x_P = v_active` held fixed (no eigenvector
        relaxation). `b` is `numerator`, and `a = omega - h(omega_fixed)`
        is the `w`-independent part `<x_P|M_PP|x_P> / <x_P|x_P>`. For a
        0th order complementary block, this root is exactly the Rayleigh
        quotient of the real matrix with the vector unfolded at it.

        Solved by Newton's method, started from the second-order solution
        around `omega_fixed`,

            omega_fixed + d / (1 + n) + s * d**2 / (1 + n)**3

        with `d = omega - omega_fixed`, `n = <x_Q|x_Q>` and
        `s = <x_Q|(omega_fixed - D)^{-1}|x_Q>`, `x_Q` being the eliminated
        part at `omega_fixed` (per unit-normalized `x_P`). Each step costs
        one divide by `w - D` and two dot products over the eliminated
        space. Iterated until a step is below `tol` (Ha), typically 2
        steps; as `g' >= 1`, `|w - root| <= |g(w)|`.

        Below the lowest eliminated `D` (the lowest pole), `g(w) = w - a -
        h(w)` increases from -inf to +inf and is convex, so it has exactly
        one root there. Newton steps that would leave the bracket known to
        contain it are replaced by bisection. If `omega` itself lies at or
        above the lowest pole, the state is not described by the fold
        (its near-degenerate configurations were eliminated), and the raw
        `omega` is kept with a warning; likewise if the iteration ends at
        the pole (a pole without coupling, `b = 0` there).

        If `output` is given, the starting value and every Newton step are
        printed to it, each with its change from the previous value (the
        starting value's from `omega`).

        Returns
        -------
        tuple(float, AmplitudeVector, int)
            The self-consistent energy; the eliminated part
            `x_Q = (w - D)^{-1} b` from the last Newton step, evaluated
            at an energy within `tol` of the returned one (so no further
            evaluation is needed to unfold); and the number of Newton
            steps taken. `(omega, None, 0)` if the raw `omega` is kept.
        """
        pole = self._lowest_eliminated_d0
        if omega >= pole:
            warnings.warn(f"Relinearized eigenvalue {omega:.6f} lies at or "
                          f"above the lowest eliminated configuration "
                          f"(D = {pole:.6f}), so the state is not described "
                          f"by the fold and is kept unchanged. Increase "
                          f"explicit_width or the omega guesses.")
            return omega, None, 0

        omega_fixed = self.omega_fixed
        norm2 = v_active.dot(v_active)
        x_q = self.inverter.apply(numerator, omega_fixed)
        n = x_q.dot(x_q) / norm2
        s = x_q.dot(x_q / evaluate(omega_fixed - self.diagonal_0)) / norm2
        a = omega - numerator.dot(x_q) / norm2
        d = omega - omega_fixed
        w = float(omega_fixed + d / (1 + n) + s * d**2 / (1 + n)**3)
        if w >= pole:
            w = 0.5 * (omega + pole)
        if output is not None:
            print(f"    2nd order   omega = {w:.10f}   "
                  f"change = {w - omega:+.4e}", file=output)

        lower, upper = -np.inf, pole  # bracket of the root
        for n_steps in range(1, maxiter + 1):
            x_q = self.inverter.apply(numerator, w)
            g = w - a - numerator.dot(x_q) / norm2
            w_new = float(w - g / (1 + x_q.dot(x_q) / norm2))
            # From the right (g > 0), convexity keeps Newton right of the
            # root; from the left, a step can overshoot the pole.
            bisected = False
            if g > 0:
                upper = w
            else:
                lower = w
                if w_new >= upper:
                    w_new = 0.5 * (lower + upper)
                    bisected = True
            step, w = w - w_new, w_new
            if output is not None:
                print(f"    Newton {n_steps:3d}  omega = {w:.10f}   "
                      f"change = {-step:+.4e}"
                      + ("   (bisection)" if bisected else ""), file=output)
            if abs(step) < tol:
                break
        else:
            warnings.warn(f"Self-consistent energy not converged within "
                          f"{maxiter} Newton steps (last step {step:.2e}).")

        if pole - w < tol:
            warnings.warn(f"Self-consistent energy of the relinearized "
                          f"eigenvalue {omega:.6f} ends at the lowest "
                          f"eliminated configuration (D = {pole:.6f}), so the "
                          f"state is kept unchanged. Increase explicit_width "
                          f"or the omega guesses.")
            return omega, None, 0
        return w, x_q, n_steps

    def _unfold(self, v_active, x_q):
        """Full, unit-normalized eigenvector: the explicit part `v_active`
        with the eliminated configurations' amplitudes `x_q` (from
        `_self_consistent_energy`) added back into `complementary_space`."""
        space = self.complementary_space
        # native Tensor addition, no densification
        full = v_active.copy()
        full[space] = evaluate(getattr(full, space) + getattr(x_q, space))
        return full / np.sqrt(full.dot(full))

    def diagonal(self):
        # Eliminated configurations are no unknowns: a large diagonal
        # keeps guesses and preconditioning away from them. Built from the
        # mask tensors, without densification.
        out = super().diagonal()
        space = self.complementary_space
        tensor = getattr(out, space)
        new_tensor = evaluate(
            tensor * self._mask_active_tensor + 1e5 * self._mask_inactive_tensor
        )
        return AmplitudeVector(**{
            **{b: t for b, t in out.items() if b != space},
            space: new_tensor,
        })
