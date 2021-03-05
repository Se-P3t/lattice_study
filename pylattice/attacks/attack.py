"""

"""
import os
import time
from fractions import Fraction

from fpylll import IntegerMatrix

from ..util import matrix_overview, str_mat
from pylattice.algorithms.lll import run_LLL
from pylattice.algorithms.bkz import run_BKZ, run_BKZ2, run_DBKZ


class BasicLatticeObject:
    """
    """
    _epsilon = 0.1
    _delta = 0.99
    _eta = 0.501

    def __init__(self, verbose=0, debug=False):
        self.verbose = verbose
        self.debug = debug
        self.B = None

    def set_basis(self, basis):
        """
        """
        if isinstance(basis, list):
            self.B = IntegerMatrix.from_matrix(basis, int_type="mpz")
        elif isinstance(basis, IntegerMatrix):
            self.B = basis
        else:
            raise TypeError("basis type '%s' unknown." % type(basis))

    def basis_overview(self):
        """
        overview of basis matrix
        """
        assert self.B is not None
        matrix_overview(self.B)

    @property
    def dim(self):
        """
        """
        assert self.B is not None
        return self.B.nrows

    def randomize(self):
        """
        TODO
        """
        raise NotImplementedError

    def LLL(self):
        """
        """
        assert self.B

        if self.verbose == 2:
            t0 = time.perf_counter()
            print("[-] LLL started")

        run_LLL(
            self.B,
            delta=self._delta,
            eta=self._eta,
            verbose=(self.verbose >= 3),
        )

        if self.verbose == 2:
            t1 = time.perf_counter()
            print(f"[+] LLL done :: {t1-t0:.2f} s")

    def BKZ(self, block_size=2, method='BKZ'):
        """
        :param block_size: (default: 2)
        :param method: (default: 'BKZ') one of 'BKZ', 'BKZ2' or 'DBKZ'
        """
        assert self.B
        if method == 'BKZ2':
            bkz_func = run_BKZ2
        elif method == 'DBKZ':
            bkz_func = run_DBKZ
        else:
            bkz_func = run_BKZ

        if self.verbose == 2:
            t0 = time.perf_counter()
            print(f"[-] BKZ started")

        bkz_func(
            self.B,
            block_size=block_size,
            verbose=(self.verbose >= 3),
        )

        if self.verbose == 2:
            t1 = time.perf_counter()
            print(f"[+] BKZ done :: {t1-t0:.2f} s")
