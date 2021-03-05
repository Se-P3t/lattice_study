"""
"""
__all__ = ['solve_svp']

from fpylll import IntegerMatrix, SVP


def solve_svp(A, **kwds):
    """fplll shortest_vector

    :param A: an instance of `list` or `IntegerMatrix`
    :param int_type: (default: 'mpz') an element of `fpylll.config.int_types`
    :param method: (default: 'fast') One of "fast" or "proved"
    :param pruning: (default: `True`) If ``True`` pruning parameters are computed by this function
    :param preprocess: (default: `True`) Blocksize used for preprocessing;
        if  ``True`` a block size is picked
    :param max_aux_solutions: maximum number of additional short-ish solutions to return
    :param verbose: (default: `False`) print verbose outputs
    """
    int_type = kwds.get('int_type', 'mpz')
    if isinstance(A, list):
        A = IntegerMatrix.from_matrix(A, int_type=int_type)
    elif not isinstance(A, IntegerMatrix):
        raise TypeError("Matrix `A` type '%s' unknown." % type(A))

    verbose = kwds.get('verbose', False)
    kwds['flags'] = SVP.VERBOSE if verbose else SVP.DEFAULT

    v = SVP.shortest_vector(A, **kwds)

    return v
