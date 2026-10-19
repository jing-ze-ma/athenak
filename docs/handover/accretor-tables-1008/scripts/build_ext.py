#!/usr/bin/env python3
"""Build the EXTENDED accretor opacity tables (X 0.70, Z 0.02, GS98 metals) with the AG Car ext2
procedure (lbv_1008/agcar/tables_ext/scripts/build_ext.py fergR, copied and only re-pointed):
TOPS/ATOMIC (high T) + AESOPUS 2.1 gas-only Rosseland / Ferguson+2005 Planck (low T),
blended in log T over [LB0, LB1], with a physical low-density extension below each
source's lowest VALID density, on the grid log T 2.6..8.0 (0.025), log rho -21..0 (0.05).

usage: build_ext.py <aesopus_src> <outdir>
  aesopus_src = 'web'  -> raw/aesopus_form/x036.dat (on-demand X=0.36, exact)
              = 'xint' -> precomputed GS98 gas tables X=0.35 and X=0.5, linear in X to 0.36

Low-density extension below the edge rho_e(T) of a source's valid data (per source T row):
  Rosseland: kR = kes(T,rho) + kabs_e (rho/rho_e)^s,  kabs_e = kR_e - kes(T,rho_e) (> 0),
             kes = LTE Saha electron scattering (saha_es.py), s = log-log slope of kR - kes
             over the 0.6 dex above the edge, clamped to [0, 1].  If kabs_e <= 0 the edge
             value is carried by the es scaling: kR = kR_e kes(rho)/kes(rho_e).
  Planck:    kP = kP_e (rho/rho_e)^s, s = log-log slope of kP over the 0.6 dex above the
             edge, clamped to [0, 1] (0 = bound-bound of the dominant ion stage, density-
             independent per gram; 1 = two-body continuum, ff / bf / H-).  No scattering.
TOPS rho_e(T) = its own density floor: below a T-dependent density TOPS returns the SAME
numbers (measured by querying to 1e-21: flat below log rho -15.6 at 5800 K ... -12 at
> 2.6 MK).  Where that floor lies above log rho -14 (T > 1e5 K) the old table's clamped
values for -14 <= log rho < floor are kept (bitwise), and the extension starts at -14.
AESOPUS / Ferguson rho_e(T) = log R -8 (their minimum).
"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/accretor_tables_1009/scripts')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/planck_tools')
from saha_es import Saha                       # noqa: E402
from compare import read_ferg, read_repo       # noqa: E402

KEV_K = 1.1604518e7
BASE = '/viper/ptmp2/jinma/accretor_tables_1009'
RAW = BASE + '/raw'
LT0, DLT, NT = 2.6, 0.025, 217
LD0, DLD, ND = -21.0, 0.05, 421
LB0, LB1 = 4.0, 4.2            # blend window in log T
SLOPE_DEX = 0.6
lT = np.round(LT0 + DLT*np.arange(NT), 6)
lD = np.round(LD0 + DLD*np.arange(ND), 6)
saha = Saha(BASE + '/mixture_x0.70_z0.02.txt')


def kes_row(t, lds):
    """log10 kes on the log-rho points lds at temperature t [K]."""
    lds = np.asarray(lds, float)
    return np.log10(saha.kes(t, 10**lds))


def slope(ld_e, xs, ys):
    """log-log slope of y over [ld_e, ld_e + SLOPE_DEX] (linear interpolation in xs)."""
    y0 = np.interp(ld_e, xs, ys)
    y1 = np.interp(ld_e + SLOPE_DEX, xs, ys)
    return (y1 - y0)/SLOPE_DEX


def smooth(lts, s, half=0.1):
    """boxcar mean of the finite slopes over source rows within +-half dex in log T: the
    per-row slopes are noisy and the extension reaches up to ~7 dex below the edge."""
    s = np.asarray(s, float)
    out = s.copy()
    for i, t in enumerate(lts):
        m = (abs(lts - t) <= half + 1e-9) & np.isfinite(s)
        if m.any() and np.isfinite(s[i]):
            out[i] = s[m].mean()
    return out


def extend_row(t, xs, kR, kP, ld_e, sPo=None, sRo=None):
    """one source row (log rho points xs ascending, log10 kR / kP) valid for log rho >= ld_e:
    returns log10 kR, kP on the target lD (edge-filled above xs[-1]) and slopes."""
    R = np.interp(lD, xs, kR)
    P = np.interp(lD, xs, kP)
    lo = lD < ld_e - 1e-9
    kRe = np.interp(ld_e, xs, kR)
    kPe = np.interp(ld_e, xs, kP)
    sP = float(np.clip(slope(ld_e, xs, kP), 0.0, 1.0)) if sPo is None else sPo
    P[lo] = kPe + sP*(lD[lo] - ld_e)
    # Rosseland: es + absorption
    xe = np.array([ld_e, ld_e + SLOPE_DEX])
    ke = kes_row(t, np.concatenate([xe, lD[lo]]))
    kRxe = np.array([kRe, np.interp(ld_e + SLOPE_DEX, xs, kR)])
    abs_x = 10**kRxe - 10**ke[:2]
    if abs_x[0] > 0.0:
        sR = (np.log10(abs_x[1]) - np.log10(abs_x[0]))/SLOPE_DEX if abs_x[1] > 0 else 1.0
        sR = float(np.clip(sR, 0.0, 1.0)) if sRo is None else sRo
        R[lo] = np.log10(10**ke[2:] + abs_x[0]*10**(sR*(lD[lo] - ld_e)))
    else:
        sR = np.nan
        R[lo] = kRe + ke[2:] - ke[0]
    return R, P, sP, sR


# ---------------------------------------------------------------- TOPS
def tops_field():
    lo = np.loadtxt(RAW + '/tops_lo_x0.70_10keV.dat')
    hi = np.loadtxt(RAW + '/tops_gs98_x0.70_z0.02_10keV.dat')
    d = np.vstack([lo, hi[hi[:, 1] > 1.0001e-14]])
    T = np.unique(d[:, 0])
    Rh = np.unique(d[:, 1])
    kR = np.full((len(T), len(Rh)), np.nan)
    kP = np.full((len(T), len(Rh)), np.nan)
    ti = {v: i for i, v in enumerate(T)}
    ri = {v: i for i, v in enumerate(Rh)}
    for t, r, a, b in d:
        kR[ti[t], ri[r]] = a
        kP[ti[t], ri[r]] = b
    assert np.isfinite(kR).all()
    lts, lds = np.log10(T*KEV_K), np.log10(Rh)
    j14 = int(np.argmin(abs(lds + 14.0)))
    rowsR, rowsP, info = [], [], []
    # pass 1: per-row floors and raw slopes; pass 2: apply the T-smoothed slopes
    rows = []
    for i in range(len(T)):
        j = 0
        while j + 1 < len(lds) and kR[i, j+1] == kR[i, 0] and kP[i, j+1] == kP[i, 0]:
            j += 1
        jf = j
        je = min(jf, j14)
        _, _, sP, sR = extend_row(T[i]*KEV_K, lds[je:], np.log10(kR[i, je:]),
                                  np.log10(kP[i, je:]), lds[je])
        if jf > j14:   # slopes measured on the unclamped data when the floor is above -14
            sP = float(np.clip(slope(lds[jf], lds[jf:], np.log10(kP[i, jf:])), 0, 1))
        rows.append((jf, je, sP, sR))
    sPs = smooth(lts, [r[2] for r in rows])
    sRs = smooth(lts, [r[3] for r in rows])
    for i in range(len(T)):
        jf, je, sP0, sR0 = rows[i]
        R, P, sP, sR = extend_row(T[i]*KEV_K, lds[je:], np.log10(kR[i, je:]),
                                  np.log10(kP[i, je:]), lds[je], sPs[i], sRs[i])
        rowsR.append(R)
        rowsP.append(P)
        info.append((lts[i], lds[jf], sR, sP))
    return lts, np.array(rowsR), np.array(rowsP), info


# ---------------------------------------------------------------- AESOPUS
def read_aes(fn):
    ln = open(fn).read().splitlines()
    import re
    for s in ln:
        m = re.search(r'log10\(R\) range: nre=\s*(\d+)\s+values from\s+([-0-9.]+)\s+to\s+'
                      r'([-0-9.]+)', s)
        if m:
            n, a, b = int(m.group(1)), float(m.group(2)), float(m.group(3))
    lR = np.linspace(a, b, n)
    rows = {}
    for s in ln:
        if s.startswith('#') or not s.strip():
            continue
        v = s.split()
        if len(v) == n + 1:
            rows[round(float(v[0]), 4)] = np.array(v[1:], float)
    t = np.array(sorted(rows))
    return t, lR, np.array([rows[x] for x in t])


def aes_merged(x):
    dd = RAW + '/aes21_gas/03-GS98/GS98_a0.0_OPALZ_%s/aesopus2.0_gasbroad_GS98_Z0.020000_X%s.tab'
    tA, rA, KA = read_aes(dd % ('highT', x))
    tL, rL, KL = read_aes(dd % ('lowT', x))
    tH, rH, KH = read_aes(dd % ('highR', x))
    t = np.concatenate([tL, tA])
    K1 = np.vstack([KL, KA])
    assert np.allclose(t, tH)
    return t, np.concatenate([rA, rH]), np.hstack([K1, KH])


def aes_field(src):
    if src == 'fergR':
        fd = '/viper/u2/jinma/ATHENAK/bench/m1_opac/ferguson05/ross/'
        t1, r1, F1 = read_ferg(fd + 'g98.7.02.tron')
        return t1, r1, F1, ('Ferguson et al. 2005 Rosseland g98.7.02.tron (GS98, X 0.7 Z 0.02 '
                            'directly, no X interpolation; grains included as in the Planck)')
    if src == 'web':
        t, lR, K = read_aes(RAW + '/aesopus_form/x036.dat')
        note = 'AESOPUS 2.1 web on-demand, GS98, Zref 0.02, X 0.36 (exact), gas only'
    else:
        t, lR, K1 = aes_merged('0.35')
        t2, lR2, K2 = aes_merged('0.5')
        assert np.allclose(t, t2) and np.allclose(lR, lR2)
        w = (0.36 - 0.35)/(0.5 - 0.35)
        K = (1 - w)*K1 + w*K2
        note = ('AESOPUS 2.1 precomputed GS98 gas-only (gasbroad) Z 0.02, X 0.35 and 0.5, '
                'linear in X to 0.36 (w=%.4f)' % w)
    return t, lR, K, note


def ferg_planck():
    fd = '/viper/u2/jinma/ATHENAK/bench/m1_opac/ferguson05/'
    t1, r1, F1 = read_ferg(fd + 'g98.pl.7.02.tpon')
    return t1, r1, F1


def lowT_field(src):
    tA, rA, KA, note = aes_field(src)
    tF, rF, KF = ferg_planck()
    # one row per AESOPUS T for R, per Ferguson T for P; build on target lD each
    def rows(ts, rs, K, which):
        sl = []
        for i, t in enumerate(ts):
            xs = rs + 3*t - 18.0
            _, _, sP, sR = extend_row(10**t, xs, K[i], K[i], xs[0])
            sl.append(sR if which == 'R' else sP)
        sl = smooth(ts, sl)
        out = []
        for i, t in enumerate(ts):
            xs = rs + 3*t - 18.0
            if which == 'R':
                R, _, _, _ = extend_row(10**t, xs, K[i], K[i], xs[0], None, sl[i])
                out.append(R)
            else:
                _, P, _, _ = extend_row(10**t, xs, K[i], K[i], xs[0], sl[i], None)
                out.append(P)
        return np.array(out)
    return (tA, rows(tA, rA, KA, 'R'), tF, rows(tF, rF, KF, 'P'), note,
            (rA[0], rA[-1], rF[0], rF[-1], tA[0], tA[-1], tF[0], tF[-1]))


def tinterp(ts, F, lt):
    """linear in log T between source rows, edge-filled; F (nsrc, ND) -> (len(lt), ND)."""
    out = np.empty((len(lt), F.shape[1]))
    for j in range(F.shape[1]):
        out[:, j] = np.interp(lt, ts, F[:, j])
    return out


def main():
    src, outdir = sys.argv[1], sys.argv[2]
    tag = 'ext2' if src == 'fergR' else 'ext'
    lts, TR, TP, info = tops_field()
    KR_t, KP_t = tinterp(lts, TR, lT), tinterp(lts, TP, lT)
    tA, AR, tF, FP, note, rng = lowT_field(src)
    KR_l, KP_l = tinterp(tA, AR, lT), tinterp(tF, FP, lT)
    if src == 'fergR':
        # ext2: inside the Ferguson data box use convert_ferguson.py's own interpolation
        # (bilinear in (log T, log R), its bilin) for BOTH means; the extension stays outside
        from convert_tops import bilin
        TT, DD = np.meshgrid(lT, lD, indexing='ij')
        RRg = DD - 3*TT + 18.0
        for K_, (ts, rs, F) in ((KR_l, aes_field('fergR')[:3]), (KP_l, ferg_planck())):
            m = (TT >= ts[0]) & (TT <= ts[-1]) & (RRg >= rs[0]) & (RRg <= rs[-1])
            K_[m] = bilin(ts, rs, F, TT, RRg)[m]
    w = np.clip((lT - LB0)/(LB1 - LB0), 0.0, 1.0)[:, None]
    KR = (1 - w)*KR_l + w*KR_t
    KP = (1 - w)*KP_l + w*KP_t
    # bitwise patch: log T >= LB1 and log rho >= -14 take the old table's strings
    oldR = [s for s in open(BASE + '/tables/rosseland_tops_gs98_x0.70_z0.02.txt')
            if not s.startswith('#')]
    oldP = [s for s in open(BASE + '/tables/planck_tops_gs98_x0.70_z0.02.txt')
            if not s.startswith('#')]
    sR = np.array([['%.5f\n' % v for v in row] for row in KR], dtype=object)
    sP = np.array([['%.5f\n' % v for v in row] for row in KP], dtype=object)
    j0 = int(round((-14.0 - LD0)/DLD))
    maxdev = [0.0, 0.0]
    for i in range(NT):
        if lT[i] < LB1 - 1e-9:
            continue
        for jj in range(281):
            a, b = oldR[i*281 + jj], oldP[i*281 + jj]
            maxdev[0] = max(maxdev[0], abs(float(a) - KR[i, j0 + jj]))
            maxdev[1] = max(maxdev[1], abs(float(b) - KP[i, j0 + jj]))
            sR[i, j0 + jj] = a
            sP[i, j0 + jj] = b
    print('patched region: max |computed - old| R %.2e P %.2e dex' % tuple(maxdev))
    hdr_common = [
        'GS98 metals, X=0.700 Y=0.280 Z=0.020 (accretor / Plaskett gainer); mixture '
        'mixture_x0.70_z0.02.txt (planck_tools/mixspec.py 0.70 0.02, TOPS number fractions; '
        'Ferguson GS98 scaled solar, same mixture)',
        'EXTENDED table (accretor-tables-1008, 10-09) with the AG Car ext2 procedure. Built by '
        'scripts/build_ext.py %s' % src,
        'HIGH T: LANL TOPS/ATOMIC (fetch_tops.py, T 0.0005-10 keV, rho 1e-21..1); its data '
        'are flat (clamped) below a T-dependent floor (log rho -15.8 at 5800 K, -14 at 1e5 K, '
        '-12.6 at 1 MK, -11.2 at 8 MK, -9.4 at 116 MK; there the TOPS-only table values (clamped, kept '
        'bitwise) fill -14..floor)',
        'BLEND: linear in log T of log10 kappa over log T %.2f..%.2f (low-T source below, '
        'TOPS above); log T >= %.2f and log rho >= -14 are BITWISE the TOPS-only table '
        '(tables/*_tops_gs98_x0.70_z0.02.txt, convert_tops.py)' % (LB0, LB1, LB1),
        'LOW-RHO EXTENSION below each source floor: Rosseland = Saha e-scattering + '
        'absorption x (rho/rho_e)^s; Planck = kP_e (rho/rho_e)^s; s = local log slope over '
        '0.6 dex, clamped [0,1], boxcar-averaged over +-0.1 dex in log T (see build_ext.py)']
    if src == 'fergR':
        hdr_common.append('ext2 = FERGUSON PAIR: low-T Rosseland AND Planck from Ferguson et al. 2005 '
                          'GS98 (f05.gs98.tar.gz md5 96a1b1db, f05.g98.pl.tar.gz md5 69d8b6d3, '
                          'Wichita State), interpolated as data/stellar_opac/planck_tools/'
                          'convert_ferguson.py (bilinear in log T, log R); grains included '
                          '(only below ~2000 K, under both tfloors)')
    write(outdir + '/rosseland_%s_gs98_x0.70_z0.02.txt' % tag, sR,
          ['AthenaK stellar Rosseland opacity table, log10 kappa_R [cm^2/g] on (log10 T[K], '
           'log10 rho[g/cm3]); INCLUDES electron scattering'] + hdr_common +
          ['LOW T (Rosseland): %s; valid log T %.2f..%.2f, log R %.1f..%.1f '
           '(R = rho/T6^3)' % (note, rng[4], rng[5], rng[0], rng[1]),
           'VALID BOX (data): log T %.2f..8.0 (TOPS to 8.065) [%s below log T 4.2, TOPS above 4.0, blended '
           'between]; log rho >= TOPS floor (log T >= 4.0) or log R >= %g (log T < 4.0); below: '
           'the physical extension (model, not data); log T < %.2f and log R > %g '
           '(log T < 4.0) EDGE-FILLED' % (max(rng[4], 2.6), 'Ferguson' if src == 'fergR' else
                                         'AESOPUS', rng[0], max(rng[4], 2.6), rng[1]),
           'Read with he_opac_table (he_star_m1), same format as before.'])
    write(outdir + '/planck_%s_gs98_x0.70_z0.02.txt' % tag, sP,
          ['AthenaK stellar Planck opacity table, log10 kappa_P [cm^2/g] on (log10 T[K], '
           'log10 rho[g/cm3]); ABSORPTION ONLY (no scattering)'] + hdr_common +
          ['LOW T (Planck): Ferguson et al. 2005 g98.pl.7.02.tpon (X 0.7 Z 0.02 directly); valid log T %.2f..%.2f, log R %.1f..%.1f; Ferguson INCLUDES GRAINS (matters '
           'only below ~1500-2000 K, log T < 3.3, below the run tfloor 3000 K)'
           % (rng[6], rng[7], rng[2], rng[3]),
           'VALID BOX (data): log T 2.70..8.0, log rho >= TOPS floor (log T >= 4.0) or '
           'log R >= -8 (log T < 4.0); log R > 1 at log T < 4.0 EDGE-FILLED; '
           'below the floors the physical extension (model, not data)',
           'Read with he_planck_table (he_star_m1), same format as before.'])
    np.savez(outdir + '/build_fields.npz', lT=lT, lD=lD, KR=KR, KP=KP, KR_t=KR_t,
             KP_t=KP_t, KR_l=KR_l, KP_l=KP_l,
             tops_info=np.array([(a, b, c if c == c else -9, d) for a, b, c, d in info]))
    for a, b, c, d in info[:20]:
        print('TOPS row log T %.3f floor %.2f sR %.2f sP %.2f' % (a, b, c, d))


def write(fn, S, comments):
    with open(fn, 'w') as f:
        for c in comments:
            f.write('# %s\n' % c)
        f.write('# grid: nT nD lTmin dlT lDmin dlD\n')
        f.write('# %d %d %g %g %g %g\n' % (NT, ND, LT0, DLT, LD0, DLD))
        f.write('# then nT*nD rows, T slowest: log10 kappa [cm^2/g]\n')
        for i in range(NT):
            f.write(''.join(S[i]))
    print('wrote', fn)


if __name__ == '__main__':
    main()
