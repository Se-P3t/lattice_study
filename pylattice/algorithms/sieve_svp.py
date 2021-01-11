"""
G6K Exact-SVP Solver

References:
    g6k/svp_exact_find_norm.py
    g6k/svp_exact.py
"""
__all__ = ['find_norm', 'solve_svp']

import os
import sys
import copy
from collections import OrderedDict

from fpylll.util import gaussian_heuristic
from fpylll import BKZ as fplll_bkz
from fpylll.tools.bkz_stats import dummy_tracer
from fpylll import Enumeration, EnumerationError

from g6k.algorithms.ducas18 import ducas18
from g6k.algorithms.workout import workout
from g6k.siever import Siever
from g6k import SieverParams
from g6k.utils.cli import run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import load_matrix_file, db_stats
from g6k.utils.util import sanitize_params_names

from ..util import str_mat
from ..mod_func import print_stats



def find_norm_kernel_trial(arg0, params=None, seed=None, goal_r0=None):
    # Pool.map only supports a single parameter
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)
    dim4free_dec = params.pop("workout/dim4free_dec")
    pump_params = pop_prefixed_params("pump", params)
    load_matrix = params.pop("load_matrix")
    verbose = params.pop("verbose")

    A, _ = load_matrix_file(load_matrix, randomize=True, seed=None, float_type="double")

    if A.nrows != n:
        raise ValueError(f"wrong dimension:: Expected dim(A) = {A.nrows}, got n = {n}")

    g6k = Siever(A, params, seed=seed)
    tracer = SieveTreeTracer(g6k, root_label=("svp-challenge", n), start_clocks=True)

    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])
    ds = list(range(0, n - 40, dim4free_dec))[::-1] + 10*[0]

    if goal_r0 is None:
        goal_r0 = 1.1 * gh

    if verbose and n < 90:
        verbose = False

    for d in ds:
        workout(g6k, tracer, 0, n, dim4free_dec=dim4free_dec, goal_r0=goal_r0*1.001,
                pump_params=pump_params, verbose=verbose)

    tracer.exit()
    return int(g6k.M.get_r(0, 0)), gh


def find_norm_kernel(arg0, params=None, seed=None):
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    n_matches = params.pop('n_matches')

    goal_r0 = None
    matches = 1
    trials = 0
    while matches < n_matches:
        trials += 1
        found_r0, gh = find_norm_kernel_trial(arg0, goal_r0=goal_r0)
        if found_r0 == goal_r0:
            matches += 1
        else:
            matches = 1
            goal_r0 = found_r0
        print("\t", n, "Trial %3d, found norm %10d = %.4f*gh, consec matches %d/%d" % (
              trials, goal_r0, goal_r0/gh, matches, n_matches))

    return goal_r0


def find_norm(load_matrix, n, threads=4, **kwds):
    """
    find `goal_r0` for a given matrix

    :param load_matrix: file that keeps the matrix
    :param n: dimension of the matrix
    :param threads: ... (default: 4)
    :param workout__dim4free_dec: By how much do we decreaseee dim4free at each iteration (default: 2)
    :param verbose: ... (default: False)
    :param n_matches: max match trials for kernel (default: 5)
    """
    lower_bound = n # lowest lattice dimension to consider (inclusive)
    upper_bound = 0 # upper bound on lattice dimension to consider (exclusive)
    step_size = 2 # increment lattice dimension in these steps
    trials = 1 # number of experiments to run per dimension
    workers = 1 # number of parallel experiments to run
    seed = int.from_bytes(os.urandom(8), 'big') # randomness seed

    workout__dim4free_dec = kwds.get('workout__dim4free_dec', 2)
    verbose = kwds.get('verbose', False)
    n_matches = kwds.get('n_matches', 5)

    params = SieverParams(load_matrix=load_matrix,
                          n_matches=n_matches,
                          threads=threads,
                          verbose=verbose)
    params['workout/dim4free_dec'] = workout__dim4free_dec

    all_params = OrderedDict({f"'threads': {threads}, ": params})

    res = run_all(find_norm_kernel, list(all_params.values()),
                  lower_bound=lower_bound,
                  upper_bound=upper_bound,
                  step_size=step_size,
                  trials=trials,
                  workers=workers,
                  seed=seed)

    # __import__('IPython').embed()
    goal_r0 = list(res.values())[0][0]

    return goal_r0


GRADIENT_BLOCKSIZE = 31
NPS = 60*[2.**29] + 5 * [2.**27] + 5 * [2.**26] + 1000 * [2.**25]


# Re-implement bkz2.svp_reduction, with a precise radius goal rather than success proba
def svp_enum(bkz, params, goal):
    n = bkz.M.d
    r = [bkz.M.get_r(i, i) for i in range(0, n)]
    gh = gaussian_heuristic(r)

    rerandomize = False
    while bkz.M.get_r(0, 0) > goal:
        if rerandomize:
            bkz.randomize_block(0, n)
        bkz.svp_preprocessing(0, n, params)

        strategy = params.strategies[n]
        radius = goal
        pruning = strategy.get_pruning(goal, gh)

        try:
            enum_obj = Enumeration(bkz.M)
            max_dist, solution = enum_obj.enumerate(0, n, radius, 0, pruning=pruning.coefficients)[0]
            bkz.svp_postprocessing(0, n, solution, tracer=dummy_tracer)
            rerandomize = False
        except EnumerationError:
            rerandomize = True

        bkz.lll_obj()

    return


def svp_kernel(arg0, params=None, seed=None):
    # Pool.map only supports a single parameter
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)
    load_matrix = params.pop("load_matrix")
    alg = params.pop("svp/alg")
    goal_r0 = 1.001 * params.pop('goal_r0')
    workout_params = pop_prefixed_params("workout/", params)
    pump_params = pop_prefixed_params("pump/", params)
    verbose = params.pop("verbose")
    if verbose and alg == "workout":
        workout_params["verbose"] = True

    A, bkz = load_matrix_file(load_matrix, randomize=True, seed=seed, float_type="double")
    if verbose:
        print(("Loaded file '%s'" % load_matrix))
    g6k = Siever(A, params, seed=seed)
    tracer = SieveTreeTracer(g6k, root_label=("svp-exact", n), start_clocks=True)

    if alg == "enum":
        assert len(workout_params) + len(pump_params) == 0
        bkz_params = fplll_bkz.Param(block_size=n, max_loops=1, strategies=fplll_bkz.DEFAULT_STRATEGY,
                                     flags=fplll_bkz.GH_BND)
        svp_enum(bkz, bkz_params, goal_r0)
        flast = -1
    elif alg == "duc18":
        assert len(workout_params) + len(pump_params) == 0
        flast = ducas18(g6k, tracer, goal=goal_r0)
    elif alg == "workout":
        flast = workout(g6k, tracer, 0, n, goal_r0=goal_r0, pump_params=pump_params, **workout_params)
    else:
        raise ValueError("Unrecognized algorithm for SVP")

    r0 = bkz.M.get_r(0, 0) if alg == "enum" else g6k.M.get_r(0, 0)
    if r0 > goal_r0:
        raise ValueError('Did not reach the goal')
    if 1.002 * r0 < goal_r0:
        raise ValueError('Found a vector(%d) shorter than the goal(%d) for n=%d.'%(r0, goal_r0, n))

    if verbose:
        print("sol %d, %s" % (n, A[0]))

    tracer.exit()
    stat = tracer.trace
    stat.data["flast"] = flast
    tracer.trace.data['res'] = A

    return stat


def solve_svp(A, goal_r0 = None, threads=4, **kwds):
    """
    A G6K Exact-SVP Solver

    :param A: basis matrix (repr. in list)
    :param goal_r0: ... Quit when this is reached. If it's None,
        use the result of `find_norm` with `n_matches` = 3
        !!! may take long time running `find_norm` !!!
    :param threads: ... (default: 4)
    :param keep_tmpfile: keep the reduced matrix (default: False)
    :param load_matrix: filename for temp matrix file
    :param verbose: ... (default: True)
    :param alg: algorithm used to solve svp, choosen in 'enum', 'duc18', 'workout'(default)
    :param debug: ... (default: False)
    """
    n = len(A)
    keep_tmpfile = kwds.get('keep_tmpfile', False)
    load_matrix = kwds.get('load_matrix', f'svpchallenge-{n}.txt')
    verbose = kwds.get('verbose', True)
    alg = kwds.get('alg', 'workout')
    debug = kwds.get('debug', False)

    lower_bound = n # lowest lattice dimension to consider (inclusive)
    upper_bound = 0 # upper bound on lattice dimension to consider (exclusive)
    step_size = 2 # increment lattice dimension in these steps
    trials = 1 # number of experiments to run per dimension
    workers = 1 # number of parallel experiments to run
    seed = int.from_bytes(os.urandom(8), 'big') # randomness seed

    with open(load_matrix, 'w') as f:
        f.write(str_mat(A))

    if goal_r0 is None:
        goal_r0 = find_norm(load_matrix, n, verbose=verbose, n_matches=3)

    if verbose:
        print(f"goal_r0 = {goal_r0}")

    params = SieverParams(load_matrix=load_matrix,
                          threads=threads,
                          goal_r0=goal_r0,
                          verbose=verbose)
    params['svp/alg'] = alg

    all_params = OrderedDict({f"'threads': {threads}, ": params})

    stats = run_all(svp_kernel, list(all_params.values()),
                    lower_bound=lower_bound,
                    upper_bound=upper_bound,
                    step_size=step_size,
                    trials=trials,
                    workers=workers,
                    seed=seed)

    inverse_all_params = OrderedDict([(v, k) for (k, v) in all_params.items()])
    stats = sanitize_params_names(stats, inverse_all_params)

    fmt = "{name:20s} :: n: {n:2d}, cputime {cputime:7.4f}s, walltime: {walltime:7.4f}s, "\
          "flast: {flast:3.2f}, |db|: 2^{avg_max:.2f}"
    print_stats(fmt, stats, ("cputime", "walltime", "flast", "avg_max"),
                extractf={"avg_max": lambda n, params, stat: db_stats(stat)[0]})

    res = list(stats.values())[0][0].data['res']

    if debug:
        __import__('IPython').embed()

    if keep_tmpfile:
        with open(load_matrix, 'w') as f:
            f.write(str(res))
    else:
        os.system(f'rm -f {load_matrix}')

    return tuple(res[0])
