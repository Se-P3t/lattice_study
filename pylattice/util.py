# -*- coding: utf-8 -*-
"""
"""

def matrix_overview(BB):
    for ii in range(len(BB)):
        a = ('%03d ' % ii)
        for jj in range(len(BB[0])):
            if BB[ii][jj] == 0:
                a += ' '
            elif abs(BB[ii][jj]) == 1:
                a += '1'
            elif abs(BB[ii][jj]) < 1:
                a += '.'
            else:
                a += 'X'
            # if len(BB) < 60:
            #     a += ' '
        print(a, flush=True)
    print('', flush=True)


def str_mat(m):
    s = ''
    for v in m:
        s += str(v).replace(',', '')
        s += '\n'
    return s

