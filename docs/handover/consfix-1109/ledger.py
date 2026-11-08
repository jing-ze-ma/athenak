"""Energy ledger residual from he_star_m1 user.hst:
dEtot/dt = L_in - L_top + P_esrc - P_sponge  (pgen he_star_m1.cpp HeStarHist doc).
Residual R = dEtot/dt - RHS, trapezoid between outputs, reported as fraction of L_in."""
import sys, re
import numpy as np


def load(fn):
    hdr = None
    rows = []
    for ln in open(fn):
        if ln.startswith('#'):
            if '[1]=' in ln:
                hdr = [s.split('=')[1] for s in re.findall(r'\[\d+\]=\S+', ln)]
            continue
        p = ln.split()
        if hdr and len(p) == len(hdr):
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    a = np.array(rows)
    t = a[:, 0]
    # keep the LAST occurrence of each time after restarts (drop overwritten tails)
    keep = np.ones(len(t), bool)
    tmin = np.inf
    for i in range(len(t) - 1, -1, -1):
        if t[i] >= tmin:
            keep[i] = False
        else:
            tmin = t[i]
    return hdr, a[keep]


def main(fn, nwin=6):
    hdr, a = load(fn)
    c = {h: i for i, h in enumerate(hdr)}
    t = a[:, 0]
    E = a[:, c['Etot']]
    rhs = a[:, c['L_in']] - a[:, c['L_top']]
    if 'P_esrc' in c:
        rhs = rhs + a[:, c['P_esrc']]
    if 'P_sponge' in c:
        rhs = rhs - a[:, c['P_sponge']]
    dt = np.diff(t)
    dE = np.diff(E)
    irhs = 0.5 * (rhs[1:] + rhs[:-1]) * dt
    Lin = a[:, c['L_in']]
    Lbar = 0.5 * (Lin[1:] + Lin[:-1])
    print(f"{fn}: {len(t)} rows, t {t[0]:.4e}..{t[-1]:.4e} s")
    print(" window         t0          t1        <L_in>     dE/dt-RHS (fraction of L_in)  "
          "<L_top>/L_in  <P_esrc>/L_in  <Mdot_top>")
    edges = np.linspace(0, len(dt), nwin + 1).astype(int)
    for w in range(nwin):
        s = slice(edges[w], edges[w + 1])
        T = dt[s].sum()
        if T <= 0:
            continue
        Lm = (Lbar[s] * dt[s]).sum() / T
        res = (dE[s].sum() - irhs[s].sum()) / T
        ltop = (0.5 * (a[1:, c['L_top']] + a[:-1, c['L_top']])[s] * dt[s]).sum() / T
        pes = 0.0
        if 'P_esrc' in c:
            pes = (0.5 * (a[1:, c['P_esrc']] + a[:-1, c['P_esrc']])[s] * dt[s]).sum() / T
        md = (0.5 * (a[1:, c['Mdot_top']] + a[:-1, c['Mdot_top']])[s] * dt[s]).sum() / T
        print(f" {w:3d} {t[edges[w]]:11.4e} {t[edges[w + 1]]:11.4e} {Lm:11.4e} "
              f"{res / Lm:+.4e}   {ltop / Lm:.4f}  {pes / Lm:+.3e}  {md:+.3e}")
    T = dt.sum()
    Lm = (Lbar * dt).sum() / T
    print(f" ALL residual/L_in = {(dE.sum() - irhs.sum()) / T / Lm:+.4e}; "
          f"dEtot/(<L_in> T) = {dE.sum() / T / Lm:+.4e}; Etot/(L_in) = {E[-1] / Lm:.3e} s")


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 6)
