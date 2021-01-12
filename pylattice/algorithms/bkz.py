"""
References:
    fpylll.BKZ.reduction
    fpylll.algorithms
"""
__all__ = ['run_BKZ', 'run_BKZ2', 'run_DBKZ']

from fpylll import BKZ, IntegerMatrix
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
from fpylll.algorithms.simple_dbkz import BKZReduction as DBKZ

#from .pbkz import BKZReduction as pBKZ # deleted after 0.5.2dev version


def run_BKZ(A, block_size=2, **kwds):
    u"""run BKZ reduction for a given matrix

    :param A: Integer matrix, represent in `list`
    :param int_type: (default: 'mpz') an element of `fpylll.config.int_types`
    :param block_size: (default: 2) an integer from 1 to ``nrows``
    :param delta: (default: 0.99) LLL parameter `0.25 < δ ≤ 1`
    :param verbose: (default: `False`) print verbose outputs
    :param float_type: an element of `fpylll.config.float_types` or `None`
    :param precision: bit precision to use if `float_type` is 'mpfr'

    For more params, see docs of `fpylll.BKZ.Param`

    :returns: reduced matrix ``B``, represent in `list`
    """
    int_type = kwds.get('int_type', 'mpz')
    float_type = kwds.get('float_type', None)
    precision = kwds.get('precision', 0)
    kwds['delta'] = kwds.get('delta', 0.99)
    kwds["strategies"] = kwds.get('strategies', BKZ.DEFAULT_STRATEGY)
    kwds['flags'] = BKZ.DEFAULT
    if kwds.get('verbose', False):
        kwds['flags'] |= BKZ.VERBOSE
    if kwds.get("auto_abort", False):
        kwds["flags"] |= BKZ.AUTO_ABORT

    A = IntegerMatrix.from_matrix(A, int_type=int_type)
    BKZ.reduction(A, BKZ.Param(block_size=block_size, **kwds),
        float_type=float_type, precision=precision)
    B = [[A[i,j] for j in range(A.ncols)] for i in range(A.nrows)]

    return B


def run_BKZ2(A, block_size=2, verbose=False):
    """
    """
    flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
    if verbose:
        flags |= BKZ.VERBOSE

    A = IntegerMatrix.from_matrix(A)
    _ = BKZ2(A)(BKZ.EasyParam(block_size=block_size, flags=flags))
    B = [[A[i,j] for j in range(A.ncols)] for i in range(A.nrows)]

    return B


def run_DBKZ(A, block_size=2, verbose=False):
    """
    """
    A = IntegerMatrix.from_matrix(A)
    _ = DBKZ(A)(block_size=block_size) # minimal implementation
    B = [[A[i,j] for j in range(A.ncols)] for i in range(A.nrows)]

    return B
