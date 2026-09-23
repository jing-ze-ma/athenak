"""Reader for the rt_test_out record written by DhjRtFreezeSplit (problem/rt_test_out).

Lines: X (cell radii of the reference column), D (density), I (cut index, implicit arms
only), and per call S (scalars, see the '# S' header), T [K], E [erg/cm^3], C [erg/g/K].
"""
import numpy as np

SCOLS = ("ncall trt dt ndiffT ndiffE maxrelT passes nonconv res gap capped fallback "
         "bud0 dE0 dtS0 budmax resfinal sumVe0").split()


def read(path):
    r = {"S": [], "T": {}, "E": {}, "C": {}, "icut": 0}
    for ln in open(path):
        if ln.startswith("#"):
            continue
        t = ln.split()
        if not t:
            continue
        if t[0] == "X":
            r["x"] = np.array(t[1:], float)
        elif t[0] == "D":
            r["rho"] = np.array(t[1:], float)
        elif t[0] == "I":
            r["icut"] = int(t[1])
        elif t[0] == "S":
            r["S"].append([float(v) for v in t[1:]])
        elif t[0] in "TEC":
            r[t[0]][int(t[1])] = np.array(t[2:], float)
    s = np.array(r["S"])
    r["s"] = {k: s[:, n] for n, k in enumerate(SCOLS[:s.shape[1]])} if len(s) else {}
    return r


def last(r, q="T"):
    n = max(r[q])
    return r[q][n]
