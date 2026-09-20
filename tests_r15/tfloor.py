#!/usr/bin/env python3
"""tests_r15: floor-atmosphere temperature history (i >= 100) and max T over i >= 96.

usage: tfloor.py <arm> ...
"""
import os
import sys

import numpy as np

import cliff

TURN = 4705.0


def main():
    for arm in sys.argv[1:]:
        P = cliff.rd(os.path.join(arm, 'rt_profile.bin'))
        print('== %s' % arm)
        print('   t/turn   <T>_i>=100    Tmin_i>=100   Tmax_i>=100   maxT_i>=96'
              '    T(i=100)    T(i=120)    rho(i=120)')
        for k in range(0, len(P), max(1, len(P)//8)):
            t, r, q = P[k]
            T, rho = q[5], q[0]
            print('  %7.3f  %11.4e  %11.4e  %11.4e  %11.4e  %10.3e  %10.3e  %10.3e'
                  % (t/TURN, T[100:].mean(), T[100:].min(), T[100:].max(),
                     T[96:].max(), T[100], T[120], rho[120]))
        t, r, q = P[-1]
        print('  %7.3f  %11.4e  %11.4e  %11.4e  %11.4e  %10.3e  %10.3e  %10.3e  (last)'
              % (t/TURN, q[5][100:].mean(), q[5][100:].min(), q[5][100:].max(),
                 q[5][96:].max(), q[5][100], q[5][120], q[0][120]))


if __name__ == '__main__':
    main()
