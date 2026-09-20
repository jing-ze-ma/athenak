"""tests_r26: shared loaders for the two settled r22 wedge runs.  READ-ONLY on data."""
import os, re, struct, sys
import numpy as np

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
R22 = B + '/tests_r22'
RUNS = {'lr1t': R22 + '/lr1t', 'lr_f100t': R22 + '/lr_f100t'}
LOGS = {'lr1t': R22 + '/lr1t.log', 'lr_f100t': R22 + '/lr_f100t.log'}
RS = 2.3717e11
TURN = 4705.0
LW = 4.0776e37
LSTAR = 2.3066e38          # L_rad,out/L uses the FULL-star L
OMFRAC = 0.17678
GM = None                  # filled from the input if needed


def hst(run):
    return np.loadtxt(RUNS[run] + '/he4.hydro.hst')


_PAT = re.compile(r'erg in=(\S+) out=(\S+) \| g in=(\S+) out=(\S+)'
                  r'(?: \| L_rad,out/L = (\S+) L_rad,cut/L = (\S+))? \(t = (\S+) s\)')


def facebudget(run):
    """returns t, Ein_cum, Eout_cum, Lout_over_Lw, gin_cum, gout_cum"""
    rows = []
    for l in open(LOGS[run], errors='ignore'):
        m = _PAT.search(l)
        if m:
            rows.append([float(m.group(7)), float(m.group(1)), float(m.group(2)),
                         float(m.group(5) or 'nan'), float(m.group(3)),
                         float(m.group(4))])
    a = np.array(rows).T
    return a[0], a[1], a[2], a[3] * LSTAR / LW, a[4], a[5]


def events(run):
    ev = np.loadtxt(RUNS[run] + '/he4.log')
    cl = np.array([(int(m.group(1)), float(m.group(2)), float(m.group(3))) for m in
                   (re.search(r'cycle=(\d+) time=(\S+) dt=(\S+)', l)
                    for l in open(LOGS[run], errors='ignore')) if m])
    tev = np.interp(ev[:, 0], cl[:, 0], cl[:, 1])
    return ev, tev, cl


def rtprof(run):
    """rt_profile.bin -> t (nt,), r (n1,), a (nt, nv, n1).
    slots: 0 rho, 1 v1, 5 T, 6 eint"""
    P = []
    with open(RUNS[run] + '/rt_profile.bin', 'rb') as f:
        while True:
            hh = f.read(16)
            if len(hh) < 16:
                break
            t, n1, nv = struct.unpack('<dii', hh)
            b = f.read(8 * (n1 + nv * n1))
            if len(b) < 8 * (n1 + nv * n1):
                break
            a = np.frombuffer(b, '<f8')
            P.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1).copy()))
    t = np.array([x[0] for x in P])
    r = P[0][1]
    return t, r, np.array([x[2] for x in P])


def shellvol(r):
    rf = np.empty(len(r) + 1)
    rf[1:-1] = .5 * (r[1:] + r[:-1])
    rf[0] = 2 * r[0] - rf[1]
    rf[-1] = 2 * r[-1] - rf[-2]
    return OMFRAC * 4 * np.pi / 3 * (rf[1:]**3 - rf[:-1]**3), rf


def stats(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return dict(n=len(x), mean=x.mean(), rms=x.std(), min=x.min(), max=x.max(),
                med=np.median(x))
