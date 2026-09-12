#!/usr/bin/env python3
"""sp_excess2: one-flip-at-a-time arms for the spherical-polar MHD-module energy excess.

sp_excess established that the excess needs NO magnetic field: `mhd_b0` (bbot = 0, total
ME identically zero) reproduced the whole thing against `hyd_ctl` -- +4.53e35 erg at 0.1
rotations and a horizontal-KE ratio of 2.064.  Those two numbers are the GATE here.

Each arm flips ONE thing relative to its module's control, on BOTH modules wherever the
key exists in both.  The mechanism is named by the arm PAIR whose gap collapses: if
`X_nowb` minus `hyd_nowb` is ~0 while `mhdb0_ctl` minus `hyd_ctl` is +4.5e35, the
well-balanced source is the asymmetry.  A flip that changes BOTH arms by the same amount
and leaves the gap alone is innocent, however much it moves the absolute energy.

Usage:  python3 analysis/excess2.py [arm ...]     (run it from bench/sp_excess2)
"""

import glob
import os
import sys

ROT = 3.05e5
SAMPLE = 0.1                       # rotations -- the whole run
REF_DE, REF_KR = 4.53e35, 2.064    # sp_excess mhd_b0 - hyd_ctl at 0.1 rot

# (mhd arm, hydro arm, what the pair tests)
PAIRS = [
    ("mhdb0_ctl",       "hyd_ctl",    "control: nothing flipped (the gate)"),
    ("mhdb0_hlle",      "hyd_hlle",   "rsolver hlle in BOTH (removes hllc-vs-hlld)"),
    ("mhdb0_nowb",      "hyd_nowb",   "wellbalance_dynamic=false, wb_x1=false"),
    ("mhdb0_nocond",    "hyd_nocond", "isotropic_conduction key removed"),
    ("mhdb0_norot",     "hyd_norot",  "omega = 0 (no rotation source)"),
    ("mhdb0_noemfdiss", "hyd_ctl",    "polar_emf_diss=false (MHD-only; vs hyd_ctl)"),
]
ARMS = sorted({a for p in PAIRS for a in p[:2]})


def read_hst(arm):
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
    return (cols, rows) if rows else None


def at(hst, t):
    cols, rows = hst
    return dict(zip(cols, min(rows, key=lambda r: abs(r[0] - t))))


def keh(d):
    return d.get("2-KE", 0.0) + d.get("3-KE", 0.0)


def read_log(arm):
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


def main():
    want = set(sys.argv[1:]) or set(ARMS)
    hst = {}
    for a in ARMS:
        h = read_hst(a)
        if h is not None:
            hst[a] = h
    missing = [a for a in sorted(want) if a not in hst]
    if missing:
        print("NOT RUN (no .hst):", " ".join(missing))
    print()

    t = SAMPLE * ROT
    print("=== every arm at 0.1 rotations " + "=" * 50)
    print("{:<16s} {:>5s} {:>8s} {:>13s} {:>12s} {:>12s} {:>11s}".format(
        "arm", "rot", "dt[s]", "tot-E", "KE2+KE3", "tot-ME", "mass"))
    for a in ARMS:
        if a not in hst or a not in want:
            continue
        d = at(hst[a], t)
        print("{:<16s} {:>5.3f} {:>8.2f} {:>13.6e} {:>12.5e} {:>12.4e} {:>11.5e}".format(
            a, d["time"] / ROT, d["dt"], d["tot-E"], keh(d),
            d.get("1-ME", 0.0) + d.get("2-ME", 0.0) + d.get("3-ME", 0.0), d["mass"]))
    print()

    print("=== THE GATE: MHD(B=0) minus HYDRO within each pair " + "=" * 29)
    print("reference (sp_excess, nothing flipped): dE = {:+.3e}, KE_h ratio = {:.3f}"
          .format(REF_DE, REF_KR))
    print("a pair whose dE and KE ratio COLLAPSE toward 0 / 1.0 names the mechanism.")
    print()
    print("{:<16s} {:<12s} {:>13s} {:>9s} {:>9s} {:>11s}  {}".format(
        "mhd arm", "hydro arm", "dE = E_m-E_h", "KE ratio", "%of ref", "d_mass", "flip"))
    base = None
    for m, h, what in PAIRS:
        if m not in hst or h not in hst:
            print("{:<16s} {:<12s} {:>13s}   -- not run --   {}".format(m, h, "", what))
            continue
        dm, dh = at(hst[m], t), at(hst[h], t)
        de = dm["tot-E"] - dh["tot-E"]
        kr = keh(dm) / keh(dh) if keh(dh) else float("nan")
        if base is None:
            base = de
        print("{:<16s} {:<12s} {:>+13.4e} {:>9.3f} {:>8.1f}% {:>+11.3e}  {}".format(
            m, h, de, kr, 100.0 * de / (base if base else REF_DE),
            dm["mass"] - dh["mass"], what))
    print()

    print("=== what each flip did to the ABSOLUTE energy (arm minus its own control) ===")
    print("(a big move here with NO move in the gate above = the flip matters but is not")
    print(" the asymmetry; a small move here with a collapsed gate = the mechanism)")
    for ctl, arms in (("hyd_ctl", [a for a in ARMS if a.startswith("hyd_")]),
                      ("mhdb0_ctl", [a for a in ARMS if a.startswith("mhdb0_")])):
        if ctl not in hst:
            continue
        dc = at(hst[ctl], t)
        for a in arms:
            if a == ctl or a not in hst:
                continue
            d = at(hst[a], t)
            print("  {:<16s} - {:<11s} dE = {:>+11.3e}  dKE_h = {:>+11.3e}  "
                  "dt {:>6.2f} vs {:>6.2f}".format(
                      a, ctl, d["tot-E"] - dc["tot-E"], keh(d) - keh(dc),
                      d["dt"], dc["dt"]))
    print()

    print("=== floor / C2P counters at the end of the run " + "=" * 34)
    for a in ARMS:
        if a not in want:
            continue
        names, rows = read_log(a)
        if not rows:
            continue
        last = rows[-1]
        cum = {n: 0 for n in names[1:]}
        for r in rows:
            for n, v in zip(names[1:], r[1:]):
                cum[n] += v
        sel = [n for n in cum if "floor" in n or "c2p" in n.lower()]
        print("  {:<16s} ".format(a) +
              " ".join("{}={}".format(n, cum[n]) for n in sel))
    print()

    print("READ IT LIKE THIS:")
    print("  the control pair must reproduce {:+.2e} / {:.2f}; if it does not, the".format(
        REF_DE, REF_KR))
    print("  binary or the input drifted and nothing below means anything.")
    print("  gap collapses on exactly one pair  -> that flip IS the asymmetry")
    print("  gap survives every pair            -> it is none of these; next suspects")
    print("     are the reconstruction actually applied on sp and the c2p floor path")


if __name__ == "__main__":
    main()
