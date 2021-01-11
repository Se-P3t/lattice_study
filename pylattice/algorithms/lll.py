"""
References:
    fpylll.LLL.reduction
"""
__all__ = ['run_LLL']

from fpylll import GSO, IntegerMatrix, LLL

def run_LLL(A, delta=0.99, eta=0.501, **kwds):
    u"""run LLL reduction for a given matrix

    :param A: Integer matrix, represent in `list`
    :param delta: (default: 0.99) LLL parameter `0.25 < δ ≤ 1`
    :param eta: (default: 0.501) LLL parameter `0 ≤ η < √δ`
    :param int_type: (default: 'mpz') an element of `fpylll.config.int_types`
    :param method: one of 'wrapper', 'proved', 'heuristic', 'fast' or `None`
    :param float_type: an element of `fpylll.config.float_types` or `None`
    :param precision: bit precision to use if `float_type` is 'mpfr'
    :param verbose: (default: `False`) print verbose outputs
    :param use_siegel: (default: `False`) use Siegel's condition
        instead of Lovász's condition
    :param early_red: (default: `False`) perform early reduction

    :returns: reduced matrix ``B``, represent in `list`
    """
    kwds['delta'] = delta
    kwds['eta'] = eta
    int_type = kwds.pop('int_type', 'mpz')
    kwds['flags'] = LLL.DEFAULT
    if kwds.pop('verbose', False):
        kwds['flags'] |= LLL.VERBOSE
    if kwds.pop('use_siegel', False):
        kwds['flags'] |= LLL.SIEGEL
    if kwds.pop('early_red', False):
        kwds['flags'] |= LLL.EARLY_RED

    A = IntegerMatrix.from_matrix(A, int_type=int_type)
    LLL.reduction(A, **kwds)
    B = [[A[i,j] for j in range(A.ncols)] for i in range(A.nrows)]

    return B
