"""table_coverage.py (dhj_planets_1101/wasp189b): EOS / ck / CE / CIA table coverage for
given T(p) columns.  usage: table_coverage.py label1=prof1.txt [label2=prof2.txt ...]
(prof: rows 'log10 p [barye]  T [K]', comments '#').  Read-only on every table."""
import sys
import numpy as np
import os
MET = os.environ.get('TC_MET', '3x')    # 3x (default) | 1x | 10x
CK = {'3x': '/viper/ptmp2/jinma/dhj_planets_1101/wasp189b/ck/CK3X',
      '1x': '/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck',
      '10x': '/viper/ptmp2/jinma/kelt20_0925/ckdata10'}[MET]
EOSF = {'3x': '/viper/ptmp2/jinma/w121_remap3x_1003/eosdump/eos_table_3x.txt',
        '1x': '/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt',
        '10x': '/viper/ptmp2/jinma/wasp121_0925/met10/dump/eos_table.txt'}[MET]
print('metallicity %s: ck/CE %s, EOS dump %s' % (MET, CK, EOSF))
KB, MH = 1.380649e-16, 1.6605390666e-24


def read_ck(fn):
    with open(fn) as fh:
        fh.readline()
        v = np.array(fh.read().split(), float)
    nT, nP, nb, ng = v[:4].astype(int)
    p = 4
    T = v[p:p+nT]; p += nT
    P = v[p:p+nP]; p += nP
    p += 2*(nb + 1) + 2*ng
    return T, P, v[p:p+nT*nP*nb*ng].reshape(nT, nP, nb, ng)


def read_ce(fn):
    with open(fn) as fh:
        hdr = fh.readline().split()
        fh.readline()
        w = np.array(fh.read().split(), float)
    cT, cP = int(hdr[0]), int(hdr[1])
    return w[:cT], w[cT:cT+cP], w[cT+cP:cT+cP+cT*cP*6].reshape(cT, cP, 6)


# ---- 1. ck tables: ranges and identity of hiT2 with the native table below 6100 K
T0, P0, k0 = read_ck(CK + '/ck/Premixed_%s_g8_11.txt' % MET)
T2, P2, k2 = read_ck(CK + '/ck/Premixed_%s_g8_11_hiT2.txt' % MET)
print('ck native  : T %g..%g K (%d), p %g..%g bar (%d)' % (T0[0], T0[-1], len(T0), P0[0], P0[-1], len(P0)))
print('ck hiT2    : T %g..%g K (%d), p %g..%g bar (%d)' % (T2[0], T2[-1], len(T2), P2[0], P2[-1], len(P2)))
idx = [int(np.where(T2 == t)[0][0]) for t in T0]
d = np.abs(np.log10(np.maximum(k2[idx], 1e-99)) - np.log10(np.maximum(k0, 1e-99)))
print('hiT2 vs native at the native T nodes (<= 6100 K): max |dlog10 k| = %.3e' % d.max())
c0T, c0P, ce0 = read_ce(CK + '/CE_tables/FastChem_ck_%s_int.txt' % MET)
c2T, c2P, ce2 = read_ce(CK + '/CE_tables/FastChem_ck_%s_int_hiT2.txt' % MET)
print('CE native  : T %g..%g, p %g..%g bar; hiT2: T %g..%g, p %g..%g bar'
      % (c0T[0], c0T[-1], c0P[0], c0P[-1], c2T[0], c2T[-1], c2P[0], c2P[-1]))
j = [int(np.where(c2T == t)[0][0]) for t in c0T]
print('CE hiT2 vs native (<= 6100 K): max |rel diff| = %.3e'
      % np.max(np.abs(ce2[j] - ce0)/np.maximum(np.abs(ce0), 1e-300)))
# CE composition along T at a few p (mu, VMR H2 He H e- H-)
print('\nFastChem hiT2 equilibrium (cols: mu, x_H2, x_H, x_e-, x_H-):')
for pb in (1e-6, 1e-3, 1e-1, 1.0, 100.0):
    jp = int(np.argmin(np.abs(np.log10(c2P/pb))))
    for t in (2500, 3000, 3500, 4000, 4500, 5000, 6000):
        it = int(np.argmin(np.abs(c2T - t)))
        r = ce2[it, jp]
        print('  p %8.1e bar T %5.0f : mu %.3f H2 %.2e H %.2e e- %.2e H- %.2e'
              % (c2P[jp], c2T[it], r[0], r[1], r[3], r[4], r[5]))

# ---- 2. EOS table (dump of the w121 3x deck build; deck range log rho -19..0, log T 1.5..6)
with open(EOSF) as fh:
    fh.readline(); fh.readline()
    nx, ny, xmin, dx, ymin, dy = [float(s) for s in fh.readline()[1:].split()]
nx, ny = int(nx), int(ny)
tab = np.genfromtxt(EOSF, comments='#').reshape(ny, nx, 2)
lrho = xmin + dx*np.arange(nx); lT = ymin + dy*np.arange(ny)
bad = ~np.isfinite(tab[..., 0]) | ~np.isfinite(tab[..., 1])
rowbad = bad.any(axis=1)
print('\nEOS dump: log rho %.2f..%.2f, log T %.2f..%.2f; rows with a non-finite entry: log T <= %.2f'
      % (lrho[0], lrho[-1], lT[0], lT[-1], lT[rowbad].max() if rowbad.any() else -1))
mu = 10**lT[:, None]*KB/(MH*10**tab[..., 1])        # mu = kT / (m_u p/rho)


def eos_at(pbar, T):
    """mu and e/rho at (p, T) from the table: root in log rho along the T row."""
    iy = int(round((np.log10(T) - ymin)/dy))
    lp = np.log10(pbar*1e6)
    lpr = lrho + tab[iy, :, 1]
    if not (lpr[0] <= lp <= lpr[-1]):
        return np.nan, np.nan, 'OFF (log p range %.2f..%.2f)' % (lpr[0]/1, lpr[-1])
    x = np.interp(lp, lpr, lrho)
    return np.interp(x, lrho, mu[iy]), x, 'in'


print('EOS mu(p, T) [the H2 -> H dissociation and H ionization, from the table itself]:')
for pb in (1e-8, 1e-6, 1e-3, 1.0, 100.0):
    s = '  p %7.0e bar:' % pb
    for t in (2000, 3000, 3500, 4000, 5000, 6000, 8000):
        m, x, st = eos_at(pb, t)
        s += '  %5d K mu %.3f (lrho %.2f)' % (t, m, x) if st == 'in' else '  %5d K %s' % (t, st)
    print(s)

# ---- 3. profiles
CIA_TMAX = {'H2-H2': 3000.0, 'H2-He': 9900.0, 'H2-H': 2500.0, 'He-H': 10000.0}
for arg in sys.argv[1:]:
    lab, fn = arg.split('=')
    r = np.loadtxt(fn)
    pbar = 10**r[:, 0]/1e6; T = r[:, 1]
    print('\n== %s (%s): p %.2e..%.2e bar, T %.0f..%.0f K' % (lab, fn, pbar.min(), pbar.max(), T.min(), T.max()))
    for lo, hi in ((1e-10, 1e-8), (1e-8, 1e-6), (1e-6, 1e-3), (1e-3, 1e-1), (1e-1, 10.0), (10.0, 1e4)):
        m = (pbar >= lo) & (pbar < hi)
        if m.any():
            print('   %7.0e..%7.0e bar: T %5.0f..%5.0f K' % (lo, hi, T[m].min(), T[m].max()))
    print('   rows T > 6100 K (beyond native ck/CE data, hiT2 extension): %d of %d'
          % ((T > 6100).sum(), len(T)))
    print('   rows T > 10100 K (ck/CE clamped in T): %d' % (T > 10100).sum())
    print('   rows p < 1e-8 bar (ck/CE clamped in p): %d (T there %s)'
          % ((pbar < 1e-8).sum(), ('%.0f..%.0f K' % (T[pbar < 1e-8].min(), T[pbar < 1e-8].max()))
             if (pbar < 1e-8).any() else '-'))
    print('   rows p > 1000 bar (clamped): %d' % (pbar > 1000).sum())
    print('   rows T > 3000 K (H2-H2 CIA held flat; H2 dissociated there): %d' % (T > 3000).sum())
    off, offd, lr = 0, 0, []
    for pb, t in zip(pbar, T):
        m, x, st = eos_at(pb, t)
        if st != 'in':            # below the dump's first row (log rho -14): ideal gas with
            offd += 1             # the mu of that row (the deck table itself goes to -19)
            iy = int(round((np.log10(t) - ymin)/dy))
            x = np.log10(pb*1e6*mu[iy, 0]*MH/(KB*t))
        lr.append(x)
        off += (x < -19.0 or x > 0.0 or t < 10**1.5 or t > 1e6)
    lr = np.array(lr)
    print('   log rho range %.2f..%.2f; rows below the dump (log rho < -14, inside the deck '
          'table -19..0): %d; rows off the deck EOS table: %d' % (lr.min(), lr.max(), offd, off))
