#!/usr/bin/env python3
"""sp_excess: read every arm's history and event log and print the energy excess.

The quantity under test is the TOTAL ENERGY of a spherical-polar MHD run minus that of
the matched hydro run at the same time.  The memory note `sp-mhd-energy-excess` records
+2.96e35 erg at 0.1 rotations, +3.3e36 at 1, +1.4e37 at 5, with the horizontal kinetic
energy 3-4x the hydro value from 0.5 rotations, while the total magnetic energy is only
~1e33 -- so the excess is not field energy and has to be made by the scheme.

Usage:  python3 analysis/excess.py [arm ...]      (default: every staged arm)
Run it from bench/sp_excess.
"""

import glob
import os
import sys

ROT = 3.05e5           # rotation period [s]
SAMPLES = [0.1, 0.25, 0.5]   # rotations at which every arm is reported

ARMS = ["mhd_ctl", "hyd_ctl", "mhd_b0", "mhd_thuni", "hyd_thuni", "mhd_lowbeta"]
# each arm is compared against this hydro control (None = it IS a control)
PAIR = {
    "mhd_ctl": "hyd_ctl",
    "hyd_ctl": None,
    "mhd_b0": "hyd_ctl",
    "mhd_thuni": "hyd_thuni",
    "hyd_thuni": None,
    "mhd_lowbeta": "hyd_ctl",
}


def read_hst(arm):
    """Return (cols, rows) of the arm's .hst, or None if it has not run."""
    cand = sorted(glob.glob(os.path.join(arm, "out", "*.hst")) +
                  glob.glob(os.path.join(arm, "*.hst")))
    if not cand:
        return None
    cols, rows = [], []
    for line in open(cand[0]):
        if line.startswith("#"):
            if "[1]=time" in line:
                cols = [f.split("=")[1] for f in line[1:].split() if "=" in f]
            continue
        f = line.split()
        if f:
            rows.append([float(x) for x in f])
    if not rows:
        return None
    return cols, rows


def at(hst, t):
    """The history row nearest time t (history is dumped every 0.01 rotations)."""
    cols, rows = hst
    row = min(rows, key=lambda r: abs(r[0] - t))
    return dict(zip(cols, row))


def read_log(arm):
    """Cumulative floor / C2P counters from the event log, as {name: [(cycle, n)]}."""
    cand = sorted(glob.glob(os.path.join(arm, "out", "*.log")) +
                  glob.glob(os.path.join(arm, "*.log")))
    if not cand:
        return None, []
    names, rows = [], []
    for line in open(cand[0]):
        if line.startswith("#"):
            if "cycle" in line:
                names = line[1:].split()
            continue
        f = line.split()
        if f:
            rows.append([int(x) for x in f])
    return names, rows


def horiz_ke(d):
    return d.get("2-KE", 0.0) + d.get("3-KE", 0.0)


def tot_me(d):
    return d.get("1-ME", 0.0) + d.get("2-ME", 0.0) + d.get("3-ME", 0.0)


def main():
    want = sys.argv[1:] or ARMS
    hst = {}
    for a in ARMS:
        h = read_hst(a)
        if h is not None:
            hst[a] = h

    missing = [a for a in want if a not in hst]
    if missing:
        print("NOT RUN (no .hst):", " ".join(missing))
    print()

    print("=== absolute history, at the sampled times " + "=" * 34)
    hdr = "{:<12s} {:>5s} {:>8s} {:>12s} {:>12s} {:>12s} {:>10s}"
    print(hdr.format("arm", "rot", "dt[s]", "tot-E", "KE2+KE3", "tot-ME", "mass"))
    for a in want:
        if a not in hst:
            continue
        for s in SAMPLES:
            d = at(hst[a], s * ROT)
            print("{:<12s} {:>5.2f} {:>8.2f} {:>12.5e} {:>12.5e} {:>12.5e} {:>10.4e}"
                  .format(a, d["time"] / ROT, d["dt"], d["tot-E"], horiz_ke(d),
                          tot_me(d), d["mass"]))
    print()

    print("=== THE MEASUREMENT: arm minus its hydro control " + "=" * 26)
    print("(E excess in erg; KE ratio is (KE2+KE3)/(KE2+KE3)_hydro; the reference number")
    print(" from sp_mhd_prod vs sp_dhj_hyd is +2.96e35 at 0.1 rot and KE ratio 1.12)")
    hdr2 = "{:<12s} {:<10s} {:>5s} {:>13s} {:>10s} {:>12s}"
    print(hdr2.format("arm", "control", "rot", "E - E_ctl", "KE ratio", "tot-ME"))
    for a in want:
        c = PAIR.get(a)
        if c is None or a not in hst or c not in hst:
            continue
        for s in SAMPLES:
            d, dc = at(hst[a], s * ROT), at(hst[c], s * ROT)
            kr = horiz_ke(d) / horiz_ke(dc) if horiz_ke(dc) else float("nan")
            print("{:<12s} {:<10s} {:>5.2f} {:>+13.4e} {:>10.3f} {:>12.4e}"
                  .format(a, c, s, d["tot-E"] - dc["tot-E"], kr, tot_me(d)))
    print()

    print("=== arm minus mhd_ctl (what the change bought) " + "=" * 28)
    if "mhd_ctl" in hst:
        for a in want:
            if a in ("mhd_ctl", "hyd_ctl", "hyd_thuni") or a not in hst:
                continue
            for s in SAMPLES:
                d, dc = at(hst[a], s * ROT), at(hst["mhd_ctl"], s * ROT)
                print("{:<12s} rot {:>4.2f}  dE = {:>+12.4e}  dKE_h = {:>+12.4e}"
                      .format(a, s, d["tot-E"] - dc["tot-E"],
                              horiz_ke(d) - horiz_ke(dc)))
    print()

    print("=== floor / C2P counters (cumulative over the run) " + "=" * 24)
    for a in want:
        names, rows = read_log(a)
        if not rows:
            continue
        tot = [0] * (len(names) - 1)
        for r in rows:
            for i, v in enumerate(r[1:]):
                tot[i] += v
        print(a, " ".join("{}={}".format(n, v)
                          for n, v in zip(names[1:], tot) if v))
    print()
    print("READ IT LIKE THIS:")
    print("  mhd_b0 == hyd_ctl to round-off        -> the MHD hydro path is sound and")
    print("                                           the excess needs a field")
    print("  mhd_b0 carries the excess             -> the MHD path's HYDRO is the bug;")
    print("                                           no field required")
    print("  mhd_thuni vs hyd_thuni loses it       -> the theta stretch")
    print("  mhd_lowbeta loses it                  -> the low-beta MHD flux path")
    print("  none of the above                     -> field-dependent, grid-independent,")
    print("                                           not the low-beta branch")


if __name__ == "__main__":
    main()
