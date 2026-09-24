#!/usr/bin/env python3
"""Transverse spread of E per radial index over the bin dumps of a run (the thin-cell
seed test): max_i (max - min)/mean of E over (x2, x3).
usage: spread.py rundir [rundir ...]"""
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'runs_5b_sp_s2'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'vis', 'python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import s2lib  # noqa: E402

for rd in sys.argv[1:]:
    out = []
    for f in sorted(glob.glob(rd + '/bin/*.m1.*.bin')):
        d = s2lib.load(f)
        e = d['m1_e']
        sp = (e.max(axis=(0, 1)) - e.min(axis=(0, 1)))/e.mean(axis=(0, 1))
        out.append(f"{d['cycle']}:{sp.max():.1e}")
    print(os.path.basename(rd), ' '.join(out))
