"""idort-1019 replacement of the gone wt-littests tests_m1/common.py: load_dump(path) -> object with .time, .x1, .data{name: (NK, NJ, NI)}
(bin via bin_convert, assembled from the MeshBlock geometry; tab via the column header) and load_tab (alias)."""
import sys, re
import numpy as np
sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python')
class Dump:
    pass
def _tab(path):
    hdr = []
    with open(path) as f:
        for ln in f:
            if not ln.startswith('#'): break
            hdr.append(ln)
    t = float(re.search(r'time=\s*([0-9eE.+-]+)', hdr[0]).group(1))
    names = hdr[-1].replace('#', '').split()
    vs = [n for n in names if n not in ('gid', 'i', 'j', 'k', 'x1v', 'x2v', 'x3v')]
    a = np.loadtxt(path, comments='#', ndmin=2); off = a.shape[1] - len(vs)
    d = Dump(); d.time = t; d.x1 = a[:, 2]
    d.data = {v: a[:, off + q][None, None, :] for q, v in enumerate(vs)}
    d.var = lambda n: d.data[n].ravel()
    return d
def _bin(path):
    from bin_convert import read_binary
    g = read_binary(path)
    geo = np.array(g['mb_geometry'], dtype=float)   # x1min x1max x2min x2max x3min x3max
    names = list(g['mb_data'].keys())
    a0 = np.array(g['mb_data'][names[0]]); nb, nk, nj, ni = a0.shape
    u1, u2, u3 = (np.unique(np.round(geo[:, c], 12)) for c in (0, 2, 4))
    NI, NJ, NK = ni * len(u1), nj * len(u2), nk * len(u3)
    d = Dump(); d.time = float(g['time']); d.data = {}
    for n in names:
        A = np.array(g['mb_data'][n], dtype=float); out = np.zeros((NK, NJ, NI))
        for m in range(nb):
            i0 = int(np.searchsorted(u1, round(geo[m, 0], 12))) * ni; j0 = int(np.searchsorted(u2, round(geo[m, 2], 12))) * nj; k0 = int(np.searchsorted(u3, round(geo[m, 4], 12))) * nk
            out[k0:k0 + nk, j0:j0 + nj, i0:i0 + ni] = A[m]
        d.data[n] = out
    x1f = np.linspace(geo[:, 0].min(), geo[:, 1].max(), NI + 1); d.x1 = 0.5 * (x1f[1:] + x1f[:-1])
    d.var = lambda n: d.data[n]
    return d
def load_dump(path):
    return _tab(path) if path.endswith('.tab') else _bin(path)
load_tab = load_dump
