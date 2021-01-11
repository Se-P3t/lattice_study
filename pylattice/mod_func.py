# -*- coding: utf-8 -*-
"""
"""
from collections import OrderedDict


def print_stats(fmt, stats, keys, extractf=None):
    if extractf is None:
        extractf = {}
    for (n, params), stat in stats.items():
        kv = OrderedDict()
        for key in keys:
            if key in extractf:
                value = extractf[key](n, params, stat)
            else:
                value = sum([float(node[key]) for node in stat]) / len(stat)
            kv[key] = value

        print(fmt.format(name=params, n=n, **kv))

    return
