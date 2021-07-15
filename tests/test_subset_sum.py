# -*- coding: utf-8 -*-
"""
"""
from pylattice.attacks.subset_sum import SSP


def test_reduction():
    n = 80
    print()
    for d in range(4, 7):
        d *= 0.1
        ssp = SSP.new_random(n, d)
        for method in range(3):
            ssp.construct_matrix(method)
            ssp.B.randomize_block()

            ssp.LLL()
            ssp.BKZ(block_size=20)

            sol = ssp.extract_solution()
            print(f"n={n} d={d:.1f} method {method} :: {sol==ssp.key}")


def test_multi():
    n = 80
    d = 0.5
    print()
    for k in range(2, 5):
        ssp = SSP.new_random(n, d, k)

        ssp.construct_matrix()
        ssp.B.randomize_block()

        ssp.BKZ(block_size=10)

        sol = ssp.extract_solution()
        print(f"n={n} d={d:.1f} k={k} :: {sol==ssp.key}")


def test_sieve():
    n = 80
    print()
    for d in range(7, 10):
        d *= 0.1
        ssp = SSP.new_random(n, d)

        ssp.construct_matrix(method=1)
        ssp.B.randomize_block()

        ssp.LLL()
        ssp.BKZ(block_size=20)

        ssp.verbose = 2
        ssp.sieve(1.0, threads=4)

        sol = ssp.extract_solution()
        print(f"n={n} d={d:.1f} :: {sol==ssp.key}")
