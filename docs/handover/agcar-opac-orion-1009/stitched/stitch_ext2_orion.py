#!/usr/bin/env python3
"""Stitch the orion 10-09 LTE low-density Planck SHAPE onto ext2 below ext2's data floor.

usage: stitch_ext2_orion.py <ext2_dir> <orion_tables_dir> <outdir>
  ext2_dir  = docs/handover/agcar-opac-1009 (branch agcar-opac-1009): planck_ext2_*.txt,
              sources/ (TOPS, Ferguson 2005), scripts/planck_tools/compare.py (read_ferg)
  orion_dir = docs/handover/agcar-opac-orion-1009/tables (planck_lowrho_orion_*.txt)

Floor rho_f(T) per ext2 log T node = lowest ext2 density node holding real source data,
reproduced from scripts/build_ext.py (fergR):
  low-T source (Ferguson 2005 Planck, used for log T < 4.2): the grid node is inside the
    Ferguson box, (RRg = lD - 3 lT + 18) >= rs[0] = -8 and ts[0] <= lT <= ts[-1]; inside
    the box build_ext.py overwrites the T-interpolated rows with bilin() of the data.
  TOPS (used for log T > 4.0): per TOPS row the clamp floor jf (first density at which kR
    or kP changes, scanning up from 1e-21), je = min(jf, j(-14)); the node is real when
    lD >= lds[je] in BOTH TOPS rows bracketing lT (tinterp is linear in log T).
  blend 4.0 < lT < 4.2: real only where both sources are real (max of the two floors);
    w = 0 at lT = 4.0 (Ferguson only), w = 1 at lT = 4.2 (TOPS only).
Stitch (orion nodes are a subset of ext2 nodes): for j < j_f(T)
  log kP_new = log kP_ext2(T, rho_f) + [log kP_orion(T, rho) - log kP_orion(T, rho_f)];
  cells with rho >= rho_f unchanged (original text). The correction delta = new - ext2 is
  multiplied by a T taper: 1/3 at the orion edge nodes log T 3.450 and 4.300, 2/3 at 3.475
  and 4.275, 1 in between, 0 outside (linear ramp from 0 at the first node outside).
  Orion starts at log rho -21 = ext2's lowest node, so no hold below -21 is needed.
Rosseland: not touched; use rosseland_ext2_gs98_x0.36_z0.02.txt (b28e97c5) as is.
"""
import sys
import numpy as np

ext2_dir, orion_dir, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, ext2_dir + '/scripts/planck_tools')
from compare import read_ferg     # noqa: E402

KEV_K = 1.1604518e7
LT0, DLT, NT = 2.6, 0.025, 217
LD0, DLD, ND = -21.0, 0.05, 421
LB0, LB1 = 4.0, 4.2
lT = np.round(LT0 + DLT*np.arange(NT), 6)
lD = np.round(LD0 + DLD*np.arange(ND), 6)
TAG = 'gs98_x0.36_z0.02'


def read_table(fn):
    """mirror of HsReadOpacityTable (he_star_m1.cpp l. 218): skip empty lines, '#' lines;
    the first comment parsing as 'int int real real real real' is the grid line."""
    hdr, vals, raw, grid = [], [], [], None
    for line in open(fn):
        s = line.rstrip('\n')
        if not s:
            continue
        if s[0] == '#':
            hdr.append(line)
            if grid is None:
                t = s[1:].split()
                try:
                    grid = (int(t[0]), int(t[1])) + tuple(float(x) for x in t[2:6])
                except (ValueError, IndexError):
                    pass
            continue
        raw.append(line)
        vals.append(float(s))
    assert grid is not None and len(vals) == grid[0]*grid[1], fn
    return hdr, raw, np.array(vals).reshape(grid[0], grid[1]), grid


def tops_floors():
    """build_ext.tops_field(): per TOPS row (log T, log rho of the first real node)."""
    lo = np.loadtxt(ext2_dir + '/sources/tops_lowrho/tops_lo.dat')
    hi = np.loadtxt(ext2_dir + '/sources/tops_only/tops_gs98_x0.36_z0.02.dat')
    d = np.vstack([lo, hi[hi[:, 1] > 1.0001e-14]])
    T, Rh = np.unique(d[:, 0]), np.unique(d[:, 1])
    kR = np.full((len(T), len(Rh)), np.nan)
    kP = np.full((len(T), len(Rh)), np.nan)
    ti = {v: i for i, v in enumerate(T)}
    ri = {v: i for i, v in enumerate(Rh)}
    for t, r, a, b in d:
        kR[ti[t], ri[r]] = a
        kP[ti[t], ri[r]] = b
    lts, lds = np.log10(T*KEV_K), np.log10(Rh)
    j14 = int(np.argmin(abs(lds + 14.0)))
    fl = []
    for i in range(len(T)):
        j = 0
        while j + 1 < len(lds) and kR[i, j+1] == kR[i, 0] and kP[i, j+1] == kP[i, 0]:
            j += 1
        fl.append(lds[min(j, j14)])
    return lts, np.array(fl)


def real_mask():
    """(NT, ND) True where ext2 holds real source data (Ferguson, TOPS or their blend)."""
    TT, DD = np.meshgrid(lT, lD, indexing='ij')
    ts, rs, _ = read_ferg(ext2_dir + '/sources/ferguson05/g98.pl.35.02.tpon')
    RRg = DD - 3*TT + 18.0                         # same expression as build_ext.py
    ferg = (TT >= ts[0]) & (TT <= ts[-1]) & (RRg >= rs[0])   # upper edge irrelevant here
    lts, fl = tops_floors()
    tfl = np.full(NT, np.inf)
    for i, t in enumerate(lT):
        if t < lts[0] - 1e-9 or t > lts[-1] + 1e-9:
            continue                               # edge-filled TOPS: not data
        b = int(np.searchsorted(lts, t - 1e-9))
        a = b if abs(lts[b] - t) < 1e-9 else b - 1
        tfl[i] = max(fl[a], fl[b])
    tops = DD >= tfl[:, None] - 1e-9
    m = np.where(TT[:, :1] < LB0 + 1e-9, ferg,
                 np.where(TT[:, :1] > LB1 - 1e-9, tops, ferg & tops))
    return m, lts, fl


def main():
    hdr, raw, E, g = read_table(ext2_dir + '/planck_ext2_%s.txt' % TAG)
    _, _, O, go = read_table(orion_dir + '/planck_lowrho_orion_%s.txt' % TAG)
    assert g == (NT, ND, LT0, DLT, LD0, DLD)
    i0 = int(round((go[2] - LT0)/DLT))
    j0 = int(round((go[4] - LD0)/DLD))
    nTo, nDo = go[0], go[1]
    assert j0 == 0
    m, lts, fl = real_mask()
    jf = np.full(NT, -1)
    for i in range(NT):
        # real data must form one contiguous block in rho reaching up to the top
        js = np.nonzero(m[i])[0]
        if len(js):
            jf[i] = js[0]
            assert m[i, js[0]:].all() or not (3.4 < lT[i] < 4.35), i
    delta = np.zeros_like(E)
    Dfull = np.zeros_like(E)
    wT = np.zeros(NT)
    for io in range(nTo):
        i = i0 + io
        wT[i] = min(1.0, (io + 1)/3.0, (nTo - io)/3.0)
        f = jf[i]
        assert 0 < f < j0 + nDo, (lT[i], f)
        fo = f - j0
        new = E[i, f] + O[io, :fo] - O[io, fo]
        Dfull[i, :f] = new - E[i, :f]
        delta[i, :f] = wT[i]*Dfull[i, :f]
    N = E + delta
    S = np.array(raw, dtype=object).reshape(NT, ND)
    chg = delta != 0.0
    for i, j in zip(*np.nonzero(chg)):
        S[i, j] = '%.5f\n' % N[i, j]
    # header: same lines, same count; line 3 and the VALID BOX line amended
    h = list(hdr)
    k3 = [k for k, s in enumerate(h) if 'EXTENDED table' in s][0]
    h[k3] = ('# STITCHED ext2+orion table (10-09): ext2 = build_ext.py fergR (dc482452); below '
             'ext2\'s data floor rho_f(T) at log T 3.45..4.30 log kP = ext2(T,rho_f) + '
             'orion_LTE(T,rho) - orion_LTE(T,rho_f) (orion 10-09 planck_lowrho_orion, '
             'agcar-opac-orion-1009/stitched/stitch_ext2_orion.py), correction tapered 1/3, '
             '2/3 at the two edge T nodes; elsewhere bitwise ext2\n')
    kv = [k for k, s in enumerate(h) if s.startswith('# VALID BOX')][0]
    h[kv] = h[kv].rstrip('\n') + ('; below the floors at log T 3.45..4.30: orion LTE '
                                  'density shape (calculation, not the power law)\n')
    fn = outdir + '/planck_ext2orion_%s.txt' % TAG
    with open(fn, 'w') as f:
        f.writelines(h)
        for i in range(NT):
            f.write(''.join(S[i]))
    np.savez(outdir + '/stitch_fields.npz', lT=lT, lD=lD, E=E, N=N, delta=delta,
             Dfull=Dfull, jf=jf, wT=wT, O=O, i0=i0)
    # floor listing
    with open(outdir + '/floors.txt', 'w') as f:
        f.write('# log T, rho_f node log rho (lowest ext2 node with real data), source, '
                'taper weight; TOPS rows: log T and first real log rho\n')
        for io in range(nTo):
            i = i0 + io
            src = ('Ferguson' if lT[i] < LB0 + 1e-9 else 'TOPS' if lT[i] > LB1 - 1e-9
                   else 'blend')
            f.write('%.3f %.2f %s %.4f\n' % (lT[i], lD[jf[i]], src, wT[i]))
        for a, b in zip(lts, fl):
            if a < 4.6:
                f.write('# TOPS row %.4f floor %.2f\n' % (a, b))
    print('wrote', fn)


if __name__ == '__main__':
    main()
