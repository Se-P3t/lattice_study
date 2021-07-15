# -*- coding: utf-8 -*-
"""
"""
from fpylll import IntegerMatrix

def matrix_overview(BB):
    for ii in range(BB.nrows):
        a = ('%03d ' % ii)
        for jj in range(BB.ncols):
            if BB[ii, jj] == 0:
                a += ' '
            elif abs(BB[ii, jj]) == 1:
                a += '1'
            elif abs(BB[ii, jj]) < 1:
                a += '.'
            else:
                a += 'X'
            # if len(BB) < 60:
            #     a += ' '
        print(a, flush=True)
    print('', flush=True)


def str_mat(m):
    if isinstance(m, IntegerMatrix):
        return str(m)
    elif isinstance(m, list):
        s = ''
        for v in m:
            s += str(v).replace(',', '')
            s += '\n'
        return s
    raise TypeError(f"unknown type ({type(m)}) of input")
