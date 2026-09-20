#!/usr/bin/env python3
"""Reduce the mltfaces.txt dumps of the optically thick spherical two-stream test."""
import sys
import numpy as np

COLS = dict(i=0, r=1, T=2, p=3, Freq=10, Fdiff=11, w=12, F2s=17)


def load(path):
    d = np.loadtxt(path)
    return d


def table(path, rin, every=8):
    d = load(path)
    r, T = d[:, 1], d[:, 2]
    freq, fdiff, f2s, w = d[:, 10], d[:, 11], d[:, 17], d[:, 12]
    # H_T = T/|dT/dr| by centred differences on the DUMPED T
    dTdr = np.gradient(T, r)
    HT = T/np.abs(dTdr)
    pred = 1.0 - HT/(2.0*r)
    meas = f2s/freq
    rd = fdiff/freq
    rows = []
    for k in range(1, len(r)-1, every):
        rows.append((r[k]/rin, T[k], HT[k]/r[k], pred[k], meas[k], rd[k]))
    return np.array(rows), (r, T, HT, pred, meas, rd, w)


def mid(path, rin):
    _, (r, T, HT, pred, meas, rd, w) = table(path, rin)
    k = len(r)//2
    return r[k]/rin, pred[k], meas[k], meas[k]-pred[k], rd[k], w.max()


if __name__ == "__main__":
    rin = 1.0e12
    for p in sys.argv[1:]:
        rows, _ = table(p, rin)
        print("### " + p)
        print("  r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req")
        for a in rows:
            print("  %7.4f  %9.1f  %8.4f  %8.4f  %8.4f  %8.4f" % tuple(a))
        print("  MID: r/rin=%.3f pred=%.4f meas=%.4f resid=%+.4f Frd/Freq=%.4f wmax=%.2e"
              % mid(p, rin))
