"""cmp.py A B: how far two runs of the M1 He box / slab (tests_m1/gates) are apart.

Two classes of quantity, each measured as the max relative difference over the
history (every row; and the first row after t = 0 alone) and over the active cells of
the state in the final restart file:
  cons  the radiation and the conserved totals: history mass, tot-E, Etot, F1top,
        F1mid, F1bot, Fres; state dens, gas energy, E_rad, F_rad (vector);
  dyn   the flow: history 1/2/3-mom, 1/2/3-KE, V1max, V1mid, KEcol; state momentum
        (vector).  The velocity is the small residual of large balanced forces (gravity,
        pressure, the radiation force), so it amplifies any difference in E_rad or F_rad
        by orders of magnitude (tests_m1/gates/README.md, 2).
Scales: each column or field against its own max |a|, except the net momenta of the
history (sums that cancel to ~0: against sqrt(2 mass (1-KE + 2-KE + 3-KE))), the signed
V1mid (against V1max), and vectors (against max |vector|).
Also prints whether the final restart files are identical after <par_end> (the header
echoes the input, which differs between arms).  The bin outputs are single precision
and cannot resolve these levels.

Restart layout used: per MeshBlock the hydro conserved variables (dens, 3 momenta,
energy), then the <rad_m1> ones (E, 3 fluxes), each over the block with its ghost zones,
then face and auxiliary arrays (not compared); sizes from the parameter dump in the file.
"""
import glob
import os
import sys

import numpy as np

CONS = ("mass", "tot-E", "Etot", "F1top", "F1mid", "F1bot", "Fres")
DYN = ("1-mom", "2-mom", "3-mom", "1-KE", "2-KE", "3-KE", "V1max", "V1mid", "KEcol")


def _names(path):
    with open(path) as f:
        for ln in f:
            if "[1]=" in ln:
                return [c.split("=")[1] for c in ln[1:].split()]
    return []


def _hst(path):
    rows = [ln.split() for ln in open(path) if ln.strip() and not ln.startswith("#")]
    return np.array(rows, dtype=float)


def _rel(d, s):
    d = float(np.max(np.abs(d)))
    if s == 0.0:
        return 0.0 if d == 0.0 else np.inf
    return d / s


def _rst_state(path):
    """(nmb, 9, nx3, nx2, nx1) active cells: dens, 3 mom, E, E_rad, 3 F_rad."""
    b = open(path, "rb").read()
    i = b.find(b"<par_end>")
    par, blk = {}, ""
    for ln in b[:i].decode(errors="ignore").splitlines():
        ln = ln.split("#")[0].strip()
        if ln.startswith("<"):
            blk = ln
        elif "=" in ln:
            k, v = [x.strip() for x in ln.split("=", 1)]
            par[(blk, k)] = v
    ng = int(par[("<mesh>", "nghost")])
    nx = [int(par[("<meshblock>", "nx%d" % d)]) for d in (1, 2, 3)]
    no = [(n + 2*ng) if n > 1 else 1 for n in nx]
    nc = no[0]*no[1]*no[2]
    p = i + len(b"<par_end>\n")
    nmb = int(np.frombuffer(b[p:p + 4], dtype=np.int32)[0])
    for q in range(p, len(b) - 8):   # the data size precedes nmb blocks of that size
        ds = int(np.frombuffer(b[q:q + 8], dtype=np.uint64)[0])
        if ds > 0 and q + 8 + nmb*ds == len(b):
            break
    d = np.frombuffer(b[q + 8:], dtype=np.float64).reshape(nmb, ds//8)
    d = d[:, :9*nc].reshape(nmb, 9, no[2], no[1], no[0])
    s = [slice(ng, ng + n) if n > 1 else slice(0, 1) for n in (nx[2], nx[1], nx[0])]
    return d[:, :, s[0], s[1], s[2]]


def _last(d, pat):
    f = sorted(glob.glob(os.path.join(d, pat)))
    return f[-1] if f else None


def diff(a, b):
    """hst_cons, hst_dyn (all rows), row1_cons (the first row after t = 0: one cycle,
    before the flow amplifies anything), rst_cons, rst_dyn, rst_bitwise; each *_at
    names the column or field of the maximum."""
    out = {}
    for k in ("hst_cons", "hst_dyn", "row1_cons", "rst_cons", "rst_dyn"):
        out[k], out[k + "_at"] = 0.0, ""
    out["rst_bitwise"] = None

    def put(kind, r, where):
        if r > out[kind]:
            out[kind], out[kind + "_at"] = r, where

    for fa in sorted(glob.glob(os.path.join(a, "*.hst"))):
        fb = os.path.join(b, os.path.basename(fa))
        ha, hb = _hst(fa), _hst(fb)
        nm = _names(fa)
        if ha.shape != hb.shape:
            put("hst_cons", np.inf, os.path.basename(fa) + " rows")
            continue
        col = {n: c for c, n in enumerate(nm)}
        scale = {}
        if all(k in col for k in ("mass", "1-KE", "2-KE", "3-KE")):
            ke = sum(ha[:, col[k]] for k in ("1-KE", "2-KE", "3-KE"))
            p = np.max(np.sqrt(2.0*np.abs(ha[:, col["mass"]]*ke)))
            scale.update({"1-mom": p, "2-mom": p, "3-mom": p})
        if "V1max" in col:
            scale["V1mid"] = np.max(np.abs(ha[:, col["V1max"]]))
        for n, c in col.items():
            kind = "cons" if n in CONS else ("dyn" if n in DYN else None)
            if kind is None:
                continue
            s = scale.get(n, np.max(np.abs(ha[:, c])))
            put("hst_" + kind, _rel(ha[:, c] - hb[:, c], s), n)
            if kind == "cons" and ha.shape[0] > 1:
                put("row1_cons", _rel(ha[1, c] - hb[1, c], s), n)
    fa = _last(os.path.join(a, "rst"), "*.rst")
    if fa is not None:
        fb = os.path.join(b, "rst", os.path.basename(fa))
        pa = open(fa, "rb").read().split(b"<par_end>", 1)[-1]
        pb = open(fb, "rb").read().split(b"<par_end>", 1)[-1]
        out["rst_bitwise"] = (pa == pb)
        sa, sb = _rst_state(fa), _rst_state(fb)
        for n, grp, kind in (("dens", [0], "cons"), ("mom", [1, 2, 3], "dyn"),
                             ("E", [4], "cons"), ("Erad", [5], "cons"),
                             ("Frad", [6, 7, 8], "cons")):
            s = np.max(np.sqrt(sum(sa[:, g]**2 for g in grp)))
            put("rst_" + kind, _rel(sa[:, grp] - sb[:, grp], s), n)
    return out


if __name__ == "__main__":
    d = diff(sys.argv[1], sys.argv[2])
    print("  ".join("%s %.2e (%s)" % (k, d[k], d[k + "_at"]) for k in (
        "hst_cons", "hst_dyn", "row1_cons", "rst_cons", "rst_dyn"))
          + "  rst bitwise %s" % d["rst_bitwise"])
