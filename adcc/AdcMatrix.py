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


class RelinearizedAdcMatrix(AdcMatrix):
    """
    Relinearized ("folded"/"downfolded") ADC matrix: the self-coupling of
    one excitation space (the "complementary" space, e.g. "pphh") is
    partitioned into exactly three tiers, by proximity to `omega_guess`:

    - an "explicit" tier (closest to `omega_guess`), kept as genuine
      explicit unknowns, coupled to the rest of the matrix via the real,
      unapproximated complementary block;
    - a "1st order" tier, eliminated ("folded" into the remaining spaces)
      via one Neumann correction to the inverse of `(omega_fixed - A_QQ)`;
    - a "0th order" tier (everything else), eliminated via a plain
      diagonal divide (no Neumann correction).

    In both eliminated tiers, `A_QQ = D + V` is split into the diagonal,
    0th order (bare orbital energy) part `D` and the remainder `V`. `D` is
    spin-blind by construction (it cannot distinguish e.g. an "aaaa" from
    an "abab" spin-block, since orbital energies don't depend on spin for
    a restricted reference), so this split keeps the resulting
    approximation spin-pure, unlike using the diagonal of the full
    (correlated) block would.

    Tier membership is decided by `|D0 - omega_screen| <= cutoff`, where
    `omega_screen` is `max(omega_guess)`; the explicit cutoff is
    `explicit_width`, and the 1st-order cutoff is `explicit_width +
    order1_width` (`order1_width` is a *span* added on top of
    `explicit_width`, not itself an absolute cutoff -- so widening
    `explicit_width` alone keeps the same-sized 1st-order buffer right
    beyond it, rather than requiring `order1_width` to be kept in sync by
    hand). `pphh`-type diagonal energies always sit above any reasonable
    target state's energy, so both cutoffs are effectively one-sided,
    not a symmetric window.

    There is never any direct coupling between the two eliminated tiers.
    Coupling between the explicit tier and the eliminated tiers can be
    included or dropped via `include_coupling`.

    Since this relinearization fixes `omega_fixed` (`mean(omega_guess)`)
    once and for all, the resulting operator is linear and can be
    diagonalised with the ordinary Davidson/Lanczos solvers. The price is
    that the result may only be exact for `omega_fixed` equal to the true
    eigenvalue; for other choices it is an approximation whose quality
    improves as `omega_fixed` approaches the true eigenvalue and/or as the
    tier widths are widened to treat more configurations explicitly/at
    1st order. If the complementary block is 0th order to begin with
    (true, in particular, for the default complementary space -- i.e. the
    highest excitation class present -- of any even-level ADC(n) method,
    e.g. ADC(2) or ADC(4)), `V` is identically zero and the elimination is
    exact regardless of `omega_guess` or the tier widths.

    Use `finalize_vector` after solving to get a variationally better
    (Ritz value) estimate of the eigenvalue and a rigorous a posteriori
    residual-norm error bound, evaluated against the real, unapproximated
    matrix -- this is the intended normal usage pattern, and the default
    tier widths below were calibrated assuming it is used.

    Parameters
    ----------
    matrix : AdcMatrix
        The full ADC matrix to relinearize.
    omega_guess : float or Sequence[float]
        Approximate energy/energies of the state(s) being targeted. A
        single float for one state; a sequence when targeting several at
        once (e.g. with `n_ep > 1`). `omega_fixed` (used in the actual
        fold, `1/(omega_fixed - D)`) is set to the mean; the tier widths
        are measured from the maximum instead (`omega_screen`), since a
        configuration can be safely far from the mean while still being
        close to resonance with whichever targeted state sits at the top
        of the range -- no single shared `omega_fixed` can be equally
        good for every state, but anchoring tier membership at the
        highest one avoids under-treating configurations that matter
        specifically for it.
    explicit_width : float, optional
        Explicit-tier cutoff in Ha, measured from `max(omega_guess)`; see
        above. Default `1.0`.
    order1_width : float, optional
        1st-order-tier cutoff, as a *span* in Ha added on top of
        `explicit_width` (i.e. the actual cutoff is `explicit_width +
        order1_width`), not an absolute cutoff of its own; see above.
        Default `4.0`, i.e. `explicit_width=1.0` and `order1_width=4.0`
        together reproduce a 1st-order cutoff of `5.0` Ha. Both defaults
        were calibrated by joint worst-case testing across a 10-molecule
        ADC(3)/cc-pVDZ benchmark, targeting 5 states per molecule at once
        (not just the lowest) and 0.05 eV accuracy, using the Ritz value
        from `finalize_vector` (not the raw relinearized eigenvalue) as
        the accuracy criterion, since that is the intended normal usage
        pattern. Confirmed to hold up (worst case still comfortably under
        that target) when re-tested at cc-pVTZ on the hardest cc-pVDZ
        cases. The explicit tier barely engages at all under this default
        (well under 0.1% of configurations in every tested case) --
        essentially all of the accuracy comes from the 1st-order tier
        together with the Ritz correction, not from any expensive
        exact-diagonalization treatment. This is a starting point based
        on limited testing (small molecules only), not a rigorously tuned
        choice.
    complementary_space : str, optional
        Excitation space to relinearize. Defaults to the highest
        excitation class present in `matrix` (`matrix.axis_blocks[-1]`).
    include_coupling : bool, optional
        Whether the explicit tier directly couples to the eliminated
        tiers via the real complementary block (True, default) or
        whether this coupling is neglected (False), i.e. the explicit
        tier only communicates with the eliminated tiers indirectly, via
        the other excitation spaces.
    """

    # See `explicit_width`/`order1_width` above for how these were
    # calibrated. order1_width is a *span* on top of explicit_width (see
    # update_omega_guess), so this default reproduces the same absolute
    # 1st-order/0th-order boundary (1.0 + 4.0 == 5.0 Ha) the original,
    # non-relative calibration found.
    _DEFAULT_EXPLICIT_WIDTH = 1.0
    _DEFAULT_ORDER1_WIDTH = 4.0

    def __init__(self, matrix, omega_guess, explicit_width=None,
                order1_width=None, complementary_space=None,
                include_coupling=True):
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
        # order1_width is a *span* added on top of explicit_width, not an
        # absolute cutoff -- so widening explicit_width alone always keeps
        # the same-sized 1st-order buffer right beyond it, rather than
        # requiring order1_width to be manually kept in sync (or erroring
        # out if it is not).
        self.order1_width = order1_width

        space_space = f"{complementary_space}_{complementary_space}"
        order_used = self.block_orders.get(space_space, None)
        if order_used is None:
            raise ValueError(f"The {space_space} block is not part of this "
                             "ADC matrix.")

        # Everything from here down (the bare 0th order diagonal and the
        # fluctuation V) depends only on the method/complementary_space,
        # never on omega_guess -- computed once here and reused by every
        # later update_omega_guess() call, so a better omega_guess
        # becoming available later (e.g. once real guess vectors exist)
        # never needs any of this recomputed.
        variant = "cvs" if self.is_core_valence_separated else None
        diagonal_0_block = ppmatrix.block(
            self.ground_state, [complementary_space, complementary_space],
            order=0, intermediates=self.intermediates, variant=variant
        )
        self.diagonal_0 = diagonal_0_block.diagonal
        self._diag0_arr = getattr(self.diagonal_0, complementary_space).to_ndarray()

        self.is_trivial = (order_used == 0)
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
        recompute the tier partition (`mask_active`, `window_occupancy`,
        `inverter`) accordingly. Everything else -- in particular
        `self.blocks`, so any environment-coupling term added via `+=`
        after construction -- is left untouched, unlike constructing a
        new `RelinearizedAdcMatrix` from scratch would.

        This is what lets construction be split into two stages: build
        with a throwaway placeholder `omega_guess` (e.g. `0.0`) before
        the real target energies are known -- e.g. before initial guess
        vectors have been obtained -- then call this once they are.
        """
        omega_guesses = np.atleast_1d(np.asarray(omega_guess, dtype=float))
        self.omega_fixed = float(np.mean(omega_guesses))
        self.omega_screen = float(np.max(omega_guesses))

        # Configurations are screened by proximity to omega_screen, not
        # to the raw diagonal: a configuration is "hard" because it is
        # close to resonance with the highest targeted state, not merely
        # because its bare 0th order energy happens to be small in an
        # absolute sense.
        screen_arr = np.abs(self._diag0_arr - self.omega_screen)
        self.mask_active = screen_arr <= self.explicit_width
        # order1_width is a span added on top of explicit_width, not an
        # absolute cutoff -- see its class-attribute default comment.
        order1_cutoff = self.explicit_width + self.order1_width
        mask_order1 = (~self.mask_active) & (screen_arr <= order1_cutoff)
        mask_order0 = ~(self.mask_active | mask_order1)

        if self.is_trivial:
            # V is identically zero, so the 1st order Neumann correction
            # is mathematically a no-op (and v_apply is None, since it is
            # never needed) -- merge the order-1 tier into order-0.
            mask_order0 = mask_order0 | mask_order1
            mask_order1 = np.zeros_like(mask_order1)

        elimination_windows = []
        if np.any(mask_order1):
            elimination_windows.append((mask_order1, 1))
        if np.any(mask_order0):
            elimination_windows.append((mask_order0, 0))

        n_total = self.mask_active.size
        n_explicit = int(self.mask_active.sum())
        self.window_occupancy = {
            "total": n_total,
            "explicit": {"count": n_explicit, "fraction": n_explicit / n_total},
        }
        for mask, order in elimination_windows:
            n = int(mask.sum())
            self.window_occupancy[f"order-{order}"] = {
                "count": n, "fraction": n / n_total
            }

        descr = [f"{n_explicit}/{n_total} ({n_explicit/n_total*100:.2f}%) explicit"]
        descr += [f"{int(mask.sum())}/{n_total} ({int(mask.sum())/n_total*100:.2f}%) order-{order}"
                 for mask, order in elimination_windows]
        print(f"RelinearizedAdcMatrix(omega_fixed={self.omega_screen:.4f},"
              f"{self.complementary_space}): ")
        print(", ".join(descr))

        self.inverter = ComplementaryBlockInverter(
            self.complementary_space, self.diagonal_0, self._v_apply,
            elimination_windows
        )

    def _active_masked(self, v):
        """Return a copy of `v` whose `complementary_space` block is
        restricted to the explicit (mask_active) configurations -- never
        trust whatever the eliminated windows of an input vector happen to
        hold, they are always recomputed fresh, never stored."""
        v_active = v.copy()
        active_tensor = getattr(v_active, self.complementary_space)
        active_tensor.set_from_ndarray(
            active_tensor.to_ndarray() * self.mask_active, 1e-14
        )
        return v_active

    def _build_numerator(self, v_active):
        """Build the "external" numerator driving the eliminated windows --
        the coupling of the other excitation spaces into `complementary_space`
        plus, if `include_coupling`, the active->inactive part of the real
        complementary block -- evaluated at the already active-masked
        `v_active`, and restricted to the non-explicit (~mask_active)
        configurations. Also returns the `complementary_space`-coupling-only
        `result` AmplitudeVector (everything matvec() needs beyond the
        numerator itself: the other spaces' contributions, with `space`
        already reduced to the explicit window's own exact self-coupling).
        """
        space = self.complementary_space
        space_space = f"{space}_{space}"

        # Apply every block except the complementary space's own
        # self-coupling. Afterwards, result's `space` block holds exactly
        # the coupling of the *other* spaces into `space` (evaluated at
        # the active/window-1 amplitudes) -- this is the "external" part
        # of the numerator driving the eliminated windows.
        result = v_active.zeros_like()
        for name, block_fn in self.blocks.items():
            if name != space_space:
                result += block_fn(v_active)

        numerator_arr = getattr(result, space).to_ndarray().copy()

        # The real, unapproximated complementary block applied to the
        # active part: its window-1 output is window-1's own (exact)
        # self-coupling; its window-2/3 output is the direct
        # active->inactive coupling (only used if include_coupling).
        self_coupling_arr = getattr(
            self.blocks[space_space](v_active), space
        ).to_ndarray()

        # The `space` output only ever represents window-1 (explicit)
        # unknowns -- windows 2+3 are not independent degrees of freedom
        # and must never appear in the returned vector, even transiently
        # via the "external" contribution computed above.
        result_space = getattr(result, space)
        result_space.set_from_ndarray(
            result_space.to_ndarray() * self.mask_active
            + self_coupling_arr * self.mask_active,
            1e-14
        )

        if self.include_coupling:
            numerator_arr += self_coupling_arr * (~self.mask_active)
        numerator_arr *= (~self.mask_active)
        numerator = AmplitudeVector(
            **{space: getattr(v_active, space).zeros_like()}
        )
        getattr(numerator, space).set_from_ndarray(numerator_arr, 1e-14)

        return numerator, result

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
            fed_back = getattr(
                self.blocks[space_space](v_compl), space
            ).to_ndarray()
            result_space = getattr(result, space)
            result_space.set_from_ndarray(
                result_space.to_ndarray() + fed_back * self.mask_active, 1e-14
            )

        return result

    def finalize_vector(self, ampl, omega):
        """
        Finalize a converged relinearized eigenvector: reconstruct the
        full vector, evaluate it against the real, unapproximated ADC
        matrix, and return the vector together with its Ritz value and
        residual norm with respect to that real matrix.

        The eliminated windows' amplitudes are evaluated directly (via the
        same fold used internally by `matvec`, but without contracting
        them back into the explicit space) at the actual converged
        eigenvalue `omega` -- not `omega_fixed`, which is no longer needed
        for accuracy once the true eigenvalue is known -- and added into
        the corresponding (previously implicit) entries of the vector.
        The result is renormalized to unit norm, since adding those
        previously-implicit components changes the vector's norm.

        The returned Ritz value is a variationally more accurate estimate
        of the true eigenvalue than `omega` itself: the Rayleigh quotient
        of a vector is accurate to *second* order in the vector's error,
        while `omega` (an eigenvalue of the approximate, relinearized
        operator) is only first-order accurate. The residual norm gives a
        rigorous a posteriori error bound for free: since the ADC matrix
        is Hermitian, some true eigenvalue is guaranteed to lie within
        `residual_norm` of the returned Ritz value -- without ever
        needing a separately converged reference calculation.

        This method deliberately stops here and applies no further
        correction (e.g. a preconditioned residual update) to the
        returned vector. If a fully converged solution of the real,
        unapproximated matrix is needed, the right way to get it is a
        separate call to a Davidson/Lanczos solver on that real matrix,
        using the returned `full_vector` as its initial guess -- this
        reuses the solver's own, already-validated convergence and
        orthogonalization logic instead of duplicating a bespoke
        correction step here, and the excellent warm start from this
        method should mean very few further iterations are needed:

        >>> full, ritz, resnorm = rel.finalize_vector(ampl, omega)
        >>> state = jacobi_davidson(matrix, [full], n_ep=1, ...)  # doctest: +SKIP

        (`matrix` being the real `AdcMatrix` that `rel` relinearizes, not
        `rel` itself.)

        Parameters
        ----------
        ampl : AmplitudeVector
            A converged eigenvector of this matrix (e.g. an entry of
            `state.eigenvectors`). Only its explicit-window entries in
            `complementary_space` are used; anything outside the explicit
            window is discarded and recomputed from scratch.
        omega : float
            The converged eigenvalue belonging to `ampl`.

        Returns
        -------
        tuple(AmplitudeVector, float, float)
            `(full_vector, ritz_value, residual_norm)`: the full,
            unit-normalized eigenvector (including the eliminated
            windows' amplitudes); its Rayleigh quotient against the real
            (unfolded) ADC matrix; and the norm of the corresponding
            residual against that real matrix.
        """
        space = self.complementary_space

        v_active = self._active_masked(ampl)
        numerator, _ = self._build_numerator(v_active)
        v_compl = self.inverter.apply(numerator, omega)

        full = v_active.copy()
        full_space = getattr(full, space)
        full_space.set_from_ndarray(
            full_space.to_ndarray() + getattr(v_compl, space).to_ndarray(),
            1e-14
        )
        full = full / np.sqrt(full.dot(full))

        # The real, unapproximated matrix-vector product -- self.blocks
        # holds the true blocks (RelinearizedAdcMatrix is constructed with
        # the same block_orders as the matrix it relinearizes), so the
        # base AdcMatrix.matvec (bypassing this class's own override)
        # gives exactly that, with no separate matrix reference needed.
        a_full = super().matvec(full)
        ritz_value = full.dot(a_full)
        residual = a_full - ritz_value * full
        residual_norm = np.sqrt(residual.dot(residual))

        return full, ritz_value, residual_norm

    def diagonal(self):
        # The eliminated windows are never independent unknowns; report a
        # large diagonal there so guesses/preconditioning never favour
        # them (mirrors AdcMatrixProjected's treatment of ignored blocks).
        out = super().diagonal()
        space = self.complementary_space
        tensor = getattr(out, space)
        arr = tensor.to_ndarray()
        arr = np.where(self.mask_active, arr, 1e5)
        new_tensor = tensor.zeros_like()
        new_tensor.set_from_ndarray(arr, 1e-14)
        return AmplitudeVector(**{
            **{b: t for b, t in out.items() if b != space},
            space: new_tensor,
        })
