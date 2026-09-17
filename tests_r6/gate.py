#!/usr/bin/env python3
"""Gate table for the he4_presn 3-D smoke arms.  Usage: gate.py <armdir> [label]"""
import sys
import os
import glob
import numpy as np

TURN = 4.705e3          # s, one turnover on the sphere
LSTAR = 2.3066e38
RSTAR = 2.3717e11
RTOP = 2.4057e11
VMLT = 1.45e7
LAYERS = [("FeCZ base", 0.635), ("FeCZ peak", 0.76), ("FeCZ top", 0.97)]


def read_prof(fn):
    """-> (t[nrec], x1v[nx1], q[nrec,nvar,nx1])"""
    ts, qs, x1v = [], [], None
    with open(fn, "rb") as f:
        while True:
            b = f.read(8)
            if len(b) < 8:
                break
            t = np.frombuffer(b, "<f8")[0]
            nx1, nvar = np.frombuffer(f.read(8), "<i4")
            x = np.frombuffer(f.read(8 * nx1), "<f8")
            q = np.frombuffer(f.read(8 * nvar * nx1), "<f8").reshape(nvar, nx1)
            if x1v is None:
                x1v = x.copy()
            ts.append(t)
            qs.append(q.copy())
    return np.array(ts), x1v, np.array(qs)


def read_surf(fn):
    ts, fs = [], []
    with open(fn, "rb") as f:
        while True:
            b = f.read(8)
            if len(b) < 8:
                break
            n = np.frombuffer(b, "<i8")[0]
            t = np.frombuffer(f.read(8), "<f8")[0]
            row = np.frombuffer(f.read(8 * 4 * n), "<f8").reshape(n, 4)
            # gnomonic cells are NOT equal area: dA ~ (1+x2^2+x3^2)^(-3/2)
            w = (1.0 + row[:, 1]**2 + row[:, 2]**2)**(-1.5)
            ts.append(t)
            fs.append(float((w * row[:, 3]).sum() / w.sum()))
    return np.array(ts), np.array(fs)


def main():
    d = sys.argv[1]
    lab = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(d.rstrip("/"))
    hst = sorted(glob.glob(os.path.join(d, "*.hydro.hst")))
    H = np.concatenate([np.loadtxt(h, comments="#", ndmin=2) for h in hst])
    H = H[np.argsort(H[:, 0])]
    t, dt, mass = H[:, 0], H[:, 1], H[:, 2]
    ke1, keh = H[:, 7], H[:, 8] + H[:, 9]

    ev = os.path.join(d, "he4.log")
    E = np.loadtxt(ev, comments="#", ndmin=2) if os.path.exists(ev) else None
    # the event log is indexed by CYCLE and is written only when something fired, so
    # its rows must be put on the time axis through the driver's own cycle/time prints
    logs = sorted(glob.glob(os.path.join(os.path.dirname(d.rstrip("/")) or ".",
                                         os.path.basename(d.rstrip("/")) + ".*.log")))
    import re as _re
    cyc, tim = [], []
    for lg in logs:
        for ln in open(lg, errors="ignore"):
            mm = _re.search(r"cycle=(\d+) time=(\S+)", ln)
            if mm:
                cyc.append(int(mm.group(1)))
                tim.append(float(mm.group(2)))
    Et = None
    if E is not None and len(E) and len(cyc) > 1:
        o = np.argsort(cyc)
        cyc = np.array(cyc)[o]
        tim = np.array(tim)[o]
        Et = np.interp(E[:, 0], cyc, tim)

    pt, x1v, Q = read_prof(os.path.join(d, "rt_profile.bin"))
    # the emergent luminosity: the pgen's own area-weighted face budget print,
    # L_rad,out/L = sum(area * F_top)/L, is exact where a weighted rt_surface mean is not
    import re as _re
    st, sf = [], []
    for lg in sorted(glob.glob(os.path.join(os.path.dirname(d.rstrip("/")) or ".",
                                            os.path.basename(d.rstrip("/")) + ".*.log"))):
        for ln in open(lg, errors="ignore"):
            mm = _re.search(r"L_rad,out/L = (\S+).*\(t = (\S+) s\)", ln)
            if mm:
                st.append(float(mm.group(2)))
                sf.append(float(mm.group(1)))
    o = np.argsort(st) if st else []
    st = np.array(st)[o] if len(st) else np.array([])
    sf = np.array(sf)[o] if len(sf) else np.array([])

    idx = [int(np.argmin(np.abs(x1v - f * RSTAR))) for _, f in LAYERS]
    rho0 = Q[0, 0, :]
    T0 = Q[0, 5, :]

    print("\n===== ARM %s =====" % lab)
    print("t_end = %.4g s = %.3f turnovers ; cycles(hst rows) = %d" %
          (t[-1], t[-1] / TURN, len(t)))
    print("radial nodes used: " + ", ".join(
        "%s r=%.4e (i=%d)" % (n, x1v[i], i) for (n, _), i in zip(LAYERS, idx)))

    hdr = ("  t/turn      t[s]     dt[s]    lnKE1   lnKEh   dmass    "
           "eos_fail floors  fofc  vceil  tclamp    tset   L/Lstar   " +
           "  ".join("d_rho_%d d_T_%d" % (k, k) for k in range(3)) + "   vr_rms/vMLT")
    print(hdr)
    step = 0.25 if t[-1] > TURN else 0.02
    nstep = int(np.floor(t[-1] / (step * TURN))) + 1
    for k in list(range(nstep)) + [None]:
        tt = t[-1] if k is None else k * step * TURN
        i = int(np.argmin(np.abs(t - tt)))
        j = int(np.argmin(np.abs(pt - tt)))
        row = "%7.2f %9.4g %9.3g %7.2f %7.2f %8.1e" % (
            t[i] / TURN, t[i], dt[i], np.log(max(ke1[i], 1e-300)),
            np.log(max(keh[i], 1e-300)), mass[i] / mass[0] - 1.0)
        if E is not None and len(E):
            # the event log is written on the SAME dt = 47 s cadence as the hst, so
            # row i of one is row i of the other; counts are per-row, so cumulate.
            sub = E[Et <= t[i] + 1e-9] if Et is not None else E[:i + 1]
            fails = sub[:, 5].sum()
            floors = sub[:, 1:5].sum()
            fofc = sub[:, 7].sum()
            # column 9 is eos_tclamp, column 10 eos_tset (APPENDED by this commit, so
            # the indices above are the ones tests_3d/arms/gate.py already used)
            tcl = sub[:, 9].sum() if sub.shape[1] > 9 else 0
            tst = sub[:, 10].sum() if sub.shape[1] > 10 else 0
            row += " %8d %7d %6d %6d %7d %7d" % (
                fails, floors, fofc, sub[:, 4].sum(), tcl, tst)
        else:
            row += " %8s %7s %6s %6s %7s %7s" % ("-", "-", "-", "-", "-", "-")
        if len(st):
            s = int(np.argmin(np.abs(st - tt)))
            row += " %9.3f" % sf[s]
        else:
            row += " %9s" % "-"
        for i3 in idx:
            row += " %7.3f %7.3f" % (Q[j, 0, i3] / rho0[i3] - 1.0,
                                     Q[j, 5, i3] / T0[i3] - 1.0)
        vr = np.sqrt(max(Q[j, 3, idx[1]], 0.0))
        row += " %9.3f" % (vr / VMLT)
        print(row)
    if E is not None and len(E):
        print("event-log totals: eos_fail=%d dfloor=%d efloor=%d tfloor=%d vceil=%d "
              "fofc=%d max c2p it=%d tclamp=%d tset=%d"
              % (E[:, 5].sum(), E[:, 1].sum(), E[:, 2].sum(), E[:, 3].sum(),
                 E[:, 4].sum(), E[:, 7].sum(), E[:, 6].max(),
                 E[:, 9].sum() if E.shape[1] > 9 else 0,
                 E[:, 10].sum() if E.shape[1] > 10 else 0))


main()
