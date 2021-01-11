"""
"""
__all__ = ['babai_nearest_plane', 'embedding_technique', 'solve_cvp']

from fpylll import IntegerMatrix, GSO, CVP


def babai_nearest_plane(A, w, **kwds):
    """Return lattice vector close to `v` using Babai's nearest plane algorithm.

    :param A: **reduced** basis, an instance of `list` or `IntegerMatrix`
    :param w: a tuple-like object
    :param int_type: (default: 'mpz') an element of `fpylll.config.int_types`
    :param float_type: (default: 'mpfr') an element of `fpylll.config.float_types`
    :param int_gram: See docs of `GSO.MatGSO`
    :param row_expo: See docs of `GSO.MatGSO`
    """
    int_type = kwds.get('int_type', 'mpz')
    if isinstance(A, list):
        A = IntegerMatrix.from_matrix(A, int_type=int_type)
    elif not isinstance(A, IntegerMatrix):
        raise TypeError("Matrix `A` type '%s' unknown." % type(A))

    float_type = kwds.get('float_type', 'mpfr')
    flags = GSO.DEFAULT
    if kwds.get('int_gram', False):
        flags |= GSO.INT_GRAM
    elif kwds.get('row_expo', False):
        flags |= GSO.ROW_EXPO

    M = GSO.Mat(A, flags=flags, float_type=float_type, update=True)
    v = M.babai(w, gso=True)

    return v


def embedding_technique(*args):
    raise NotImplementedError("do it by yourself")


def solve_cvp(A, w, int_type='mpz', method='proved', verbose=False):
    """fplll closest_vector

    :param A: **reduced** basis, an instance of `list` or `IntegerMatrix`
    :param w: a tuple-like object
    :param int_type: (default: 'mpz') an element of `fpylll.config.int_types`
    :param method: (default: 'proved') One of "fast" or "proved"
    :param verbose: (default: `False`) print verbose outputs
    """
    if isinstance(A, list):
        A = IntegerMatrix.from_matrix(A, int_type=int_type)
    elif not isinstance(A, IntegerMatrix):
        raise TypeError("Matrix `A` type '%s' unknown." % type(A))

    flags = CVP.VERBOSE if verbose else CVP.DEFAULT

    v = CVP.closest_vector(A, w, method=method, flags=flags)

    return v
