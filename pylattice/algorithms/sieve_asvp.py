"""
G6K Approx-SVP Solver

References:
    g6k/svp_challenge.py
"""
__all__ = ['solve_asvp']

import os
import sys
import copy
from collections import OrderedDict

from fpylll.util import gaussian_heuristic

from g6k.algorithms.workout import workout
from g6k.siever import Siever
from g6k import SieverParams
from g6k.utils.cli import run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import load_matrix_file, db_stats
from g6k.utils.util import sanitize_params_names

from ..util import str_mat
from ..mod_func import print_stats



def asvp_kernel(arg0, params=None, seed=None):
    # Pool.map only supports a single parameter
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)

    load_matrix = params.pop("load_matrix")
    goal_r0__gh = params.pop('goal_r0__gh')
    pump_params = pop_prefixed_params("pump", params)
    workout_params = pop_prefixed_params("workout", params)
    verbose = params.pop("verbose")
    if verbose:
        workout_params["verbose"] = True

    A, _ = load_matrix_file(load_matrix, randomize=False, seed=None, float_type="double")
    if verbose:
        print(("Loaded file '%s'" % load_matrix))

    g6k = Siever(A, params, seed=seed)
    tracer = SieveTreeTracer(g6k, root_label=("svp-challenge", n), start_clocks=True)

    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])
    goal_r0 = (goal_r0__gh**2) * gh
    if verbose:
        print(("gh = %f, goal_r0/gh = %f, r0/gh = %f" % (gh, goal_r0/gh, sum([x*x for x in A[0]])/gh)))
 
    flast = workout(g6k, tracer, 0, n, goal_r0=goal_r0, **workout_params)

    tracer.exit()
    stat = tracer.trace
    sol = tuple(A[0])
    stat.data["flast"] = flast
    tracer.trace.data['res'] = A

    if verbose:
        print(f"svp: sol {sol}")

    norm = sum([x*x for x in sol])
    if verbose:
        print("svp: norm %.1f ,hf %.5f" % (norm**.5, (norm/gh)**.5))

    return tracer.trace


def asvp(n, params, threads):
    lower_bound = n # lowest lattice dimension to consider (inclusive)
    upper_bound = 0 # upper bound on lattice dimension to consider (exclusive)
    step_size = 2 # increment lattice dimension in these steps
    trials = 1 # number of experiments to run per dimension
    workers = 1 # number of parallel experiments to run
    seed = int.from_bytes(os.urandom(8), 'big') # randomness seed

    all_params = OrderedDict({f"'threads': {threads}, ": params})

    stats = run_all(asvp_kernel, list(all_params.values()),
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

    return res


def solve_asvp(A, **kwds):
    """
    A G6K Approx-SVP Solver

    :param A: basis matrix (repr. in list)
    :param keep_tmpfile: keep the reduced matrix (default: False)
    :param load_matrix: filename for temp matrix file
    :param verbose: ... (default: True)
    :param workout__dim4free_dec: By how much do we decreaseee dim4free at each iteration (default: 3)
    :param goal_r0__gh: ... Quit when this is reached (default: 1.05)
    """
    n = len(A)
    keep_tmpfile = kwds.get('keep_tmpfile', False)
    threads = kwds.get('threads', 4)
    load_matrix = kwds.get('load_matrix', f'svpchallenge-{n}.txt')
    verbose = kwds.get('verbose', True)
    workout__dim4free_dec = kwds.get('workout__dim4free_dec', 3)
    goal_r0__gh = kwds.get('goal_r0__gh', 1.05)

    params = SieverParams(threads=threads,
                          load_matrix=load_matrix,
                          verbose=verbose)
    params['workout/dim4free_dec'] = workout__dim4free_dec
    params['goal_r0__gh'] = goal_r0__gh

    with open(load_matrix, 'w') as f:
        f.write(str_mat(A))

    res = asvp(n, params, threads)

    if keep_tmpfile:
        with open(load_matrix, 'w') as f:
            f.write(str(res))
    else:
        os.system(f'rm -f {load_matrix}')

    res = [[res[i,j] for j in range(res.ncols)] for i in range(res.nrows)]
    return res
