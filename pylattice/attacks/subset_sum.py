"""
solve subset sum problem

References:
    Improved low-density subset sum algorithms
    https://link.springer.com/content/pdf/10.1007/BF01201999.pdf

    A Note on the Density of the Multiple Subset Sum Problems
    https://eprint.iacr.org/2011/525.pdf

    Solving low-density multiple subset sum problems with SVP oracle
    http://www.sysmath.com/jweb_xtkxyfzx/CN/article/downloadArticleFile.do?attachType=PDF&id=12706
"""
__all__ = ["FPLLL", "SSP"]

import os
import math
import random

from pylattice import FPLLL

from ..util import matrix_overview, str_mat
from pylattice.algorithms.lll import run_LLL
from pylattice.algorithms.bkz import run_BKZ, run_BKZ2, run_DBKZ
from pylattice.algorithms.sieve_asvp import solve_asvp
from pylattice.algorithms.sieve_svp import solve_svp


SEED = int.from_bytes(os.urandom(8), 'big')
random.seed(SEED)
FPLLL.set_random_seed(SEED)
FPLLL.set_precision(0)


class SSP:
    """
    subset sum problem
    """
    _epsilon = 0.1
    _delta = 0.99
    _eta = 0.501

    def __init__(self, numss=None, sums=None, key=None, verbose=0):
        """
        :param numss:
        :param sums:
        :param key:
        :param verbose: (default: 0)
        """
        if numss and sums:
            if isinstance(sums, int):
                numss = [numss]
                sums = [sums]

        self.numss = numss
        self.sums = sums
        self.key = key
        self.verbose = verbose
        self.zeros = tuple()
        self.ones = tuple()
        self.method = None
        self.B = None

    @property
    def n(self):
        """
        number of the given positive integers (weights)
        """
        assert self.numss and self.sums
        return len(self.numss[0])

    @property
    def k(self):
        """
        number of instances of SSP with same key
        """
        assert self.numss and self.sums
        return len(self.sums)

    @property
    def d(self):
        """
        density of SSP
        """
        assert self.numss and self.sums
        bound = float(sum(map(lambda nums: max(nums).bit_length(), self.numss)))
        return float(self.n) / bound

    def __repr__(self):
        r = 'SSP(None)'
        if self.numss and self.sums:
            r = f"SSP(n={self.n}, k={self.k}) of density {self.d:.4f}"
        if self.key:
            r += ' (known key)'
        return r

    @staticmethod
    def gen(n, d, k=1, e = None, key = None):
        """
        generate (`k`-multiple) subset sum problem

        :param n:
        :param d:
        :param k:
        :param e: hamming weight of secret key (0 < `e` <= `n`)
        :param key:
        """
        if not (0 < d < 1):
            raise ValueError(f"Parameter d({d}) does not satisfy constraint  0 < d < 1")
        if e is not None and not (0 < e <= n):
            raise ValueError(f"Parameter e({e}) does not satisfy constraint  0 < e <= n = {n}")
        if key is not None and len(key) != n:
            raise ValueError(f"Parameter key(len {len(key)}) must have length of n({n})")

        bound = 2**math.ceil(n/d/k)

        if e is None:
            key = [random.randint(0, 1) for _ in range(n)]
        else:
            key = [0] * n
            for i in random.sample(range(n), e):
                key[i] = 1

        nums = [random.randint(1, bound) for _ in range(n)]
        s = sum(ai*xi for ai,xi in zip(nums, key))

        if k == 1:
            return (nums, s), tuple(key)
        else:
            numss = [nums]
            sums = [s]

            while len(sums) < k:
                nums = [random.randint(1, bound) for _ in range(n)]
                s = sum(ai*xi for ai,xi in zip(nums, key))
                if s not in sums:
                    numss.append(nums)
                    sums.append(s)

            return (numss, sums), tuple(key)

    @classmethod
    def new_random(cls, n, d, k=1, e = None, key = None):
        """
        """
        (numss, sums), key = cls.gen(n=n, d=d, k=k, e=e, key=key)
        return cls(numss=numss, sums=sums, key=key)

    def remove_partial_key(self, zeros=tuple(), ones=tuple()):
        """
        remove partial guess/known keys

        !!! TODO !!! can only be used once

        :param zeros: index of zeros, count from 0
        :param ones: index of ones, count from 0
        """
        assert self.numss and self.sums
        assert (not self.zeros) and (not self.ones)
        assert zeros or ones

        new_numss = [
            [num for idx, num in enumerate(nums) if idx not in zeros and idx not in ones]
            for nums in self.numss
        ]
        new_sums = [
            s - sum(num for idx, num in enumerate(nums) if idx in ones)
            for nums, s in zip(self.numss, self.sums)
        ]

        self.numss = new_numss
        self.sums = new_sums
        self.zeros = tuple(zeros)
        self.ones = tuple(ones)

        if self.verbose >= 1:
            print('partial key: ' + self.partial_key)

    @property
    def partial_key(self):
        """
        """
        key_s = ''
        for idx in range(self.n):
            if idx in self.zeros:
                key_s += '0'
            elif idx in self.ones:
                key_s += '1'
            else:
                key_s += '.'
        return key_s

    def construct_matrix(self, method=1):
        """
        :param method: (default: 1)
        """
        assert self.numss and self.sums
        n = self.n
        k = self.k

        if method == 0:
            """
            Brickell and Lagarias-Odlyzko
            """
            if k != 1:
                raise ValueError("not considered for solving multiple subset sum problems")

            if self.verbose >= 1 and self.d > 0.6463:
                print("density {d:.4f} > 0.6463, solutions might not be found.")

            # N^2 > sum(key)
            N = math.ceil((1+self._epsilon) * math.sqrt(n))

            B = [[0]*(n+1) for _ in range(n + 1)]
            B[0][n] = N * self.sums[0]
            for i in range(n):
                B[i+1][i] = 1
                B[i+1][n] = N * self.numss[0][i]

        elif method == 1:
            """
            Coster et al. *frequently used*
            """
            if self.verbose >= 1 and self.d > 0.9408:
                print("density {d:.4f} > 0.9408, solutions might not be found.")

            if k == 1:
                # N^2 > sum(key)/4
                N = math.ceil((1+self._epsilon) * math.sqrt(n)/2.)
            else:
                # ?
                N = math.ceil((1+self._epsilon) * math.sqrt((n+1)/4.))

            B = [[0]*(n+k) for _ in range(n + 1)]
            for j in range(k):
                B[0][n+j] = 2 * N * self.sums[j]
            for i in range(n):
                B[0][i] = 1
                B[i+1][i] = 2
                for j in range(k):
                    B[i+1][n+j] = 2 * N * self.numss[j][i]

        elif method == 2:
            """
            Coster et al.
            """
            if self.verbose >= 1 and self.d > 0.9408:
                print("density {d:.4f} > 0.9408, solutions might not be found.")

            # N big enough
            N = math.ceil((1+self._epsilon) * n**2)

            B = [[-1]*(n+1+k) for _ in range(n + 1)]
            B[0][n] = n + 1
            for j in range(k):
                B[0][n+1+j] = -N * self.sums[j]
            for i in range(n):
                B[i+1][i] = n + 1
                for j in range(k):
                    B[i+1][n+1+j] = N * self.numss[j][i]

        else:
            raise RuntimeError(f"unknown method id {method}")

        self.method = method
        self.B = B

        if self.verbose >= 2:
            matrix_overview(self.B)

    @property
    def basis(self):
        """
        overview of basis matrix
        """
        matrix_overview(self.B)

    def extract_solution(self):
        """
        """
        assert self.method is not None
        n = self.n

        if self.method == 0:
            """
            (e_1, ..., e_n, 0)
            """
            for i, b in enumerate(self.B):
                if self.verbose >= 1:
                    print(f"\rchecking... {i+1}/{n+1}", end='')

                if any(x != 0 for x in b[n:]):
                    continue
                if any(abs(x) > 1 for x in b[:n]):
                    continue

                sol = tuple(map(abs, b[:n]))
                #if self.verbose >= 5:
                #    print(i, sol)

                if self.check_sol(sol):
                    if self.verbose >= 1:
                        print('\nfind solution:', sol)
                    return sol

            else:
                if self.verbose >= 1:
                    print("\nnot found")
                return None

        elif self.method == 1:
            """
            (2e_1-1, ..., 2e_n-1, 0, ...
            """
            for i, b in enumerate(self.B):
                if self.verbose >= 1:
                    print(f"\rchecking... {i+1}/{n+1}", end='')

                if any(x != 0 for x in b[n:]):
                    continue
                if any(abs(x) != 1 for x in b[:n]):
                    continue

                sol_0 = tuple([0 if x == 1 else 1 for x in b[:n]])
                sol_1 = tuple([1 if x == 1 else 0 for x in b[:n]])

                if self.check_sol(sol_0):
                    if self.verbose >= 1:
                        print('\nfind solution:', sol_0)
                    return sol_0
                if self.check_sol(sol_1):
                    if self.verbose >= 1:
                        print('\nfind solution:', sol_1)
                    return sol_1

            else:
                if self.verbose >= 1:
                    print("\nnot found")
                return None

        elif self.method == 2:
            """
            """
            for i, b in enumerate(self.B):
                if self.verbose >= 1:
                    print(f"\rchecking... {i+1}/{n+1}", end='')

                if any(x != 0 for x in b[n+1:]):
                    continue
                if len(set(b[:n+1])) != 2:
                    continue

                sol = tuple(1 if x == b[n] else 0 for x in b[:n])

                if self.check_sol(sol):
                    if self.verbose >= 1:
                        print('\nfind solution:', sol)
                    return sol

            else:
                if self.verbose >= 1:
                    print("\nnot found")
                return None

        else:
            raise RuntimeError(f"unknown method id {self.method}")

    def check_sol(self, sol):
        """
        """
        return all(sum(num_i*sol_i for num_i, sol_i in zip(nums, sol)) == s
                   for nums, s in zip(self.numss, self.sums))

    def shuffle(self, start=0, end = None):
        """
        :param start: (default: 0)
        :param end: (default: `n-1`)
        """
        B = self.B
        if end is None:
            end = self.n - 1

        self.B = (B[:start]
                  +random.shuffle(B[start:end])
                  +B[end+1:])

        if self.verbose >= 2:
            matrix_overview(self.B)

    def LLL(self):
        """
        """
        assert self.B

        self.B = run_LLL(
            self.B,
            delta=self._delta,
            eta=self._eta,
            verbose=(self.verbose >= 3),
        )

        if self.verbose >= 2:
            matrix_overview(self.B)

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

        self.B = bkz_func(
            self.B,
            block_size=block_size,
            verbose=(self.verbose >= 3),
        )

        if self.verbose >= 2:
            matrix_overview(self.B)

    def _sieve(self, method='svp', **kwds):
        """
        !!! TODO !!!
        """
        kwds['verbose'] = kwds.get('verbose', self.verbose >= 3)

        if method == 'svp':
            if self.verbose >= 1:
                print("may need lots of time to find norm")
            self.B = solve_svp(self.B, **kwds)

        elif method == 'asvp':
            if self.verbose >= 1:
                print("need to guess the goal_r0/gh")
            self.B = solve_asvp(self.B, **kwds)

        else:
            raise ValueError(f"unknown method {method}")
