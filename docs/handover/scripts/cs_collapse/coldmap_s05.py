"""si/s05 post-mortem: cold / drained column maps at rot 1..5, coldmap.py's method,
plus the meshblock-gid mapping needed to point problem/diag_gid at the abort column."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/analysis')
import bin_convert
import dhjcs
import csgeom

B = '/viper/u2/jinma/ATHENAK/bench/cs_ens/'
OUT = B + 'analysis/si_s05/'
EOS = dhjcs.EOS(B + 'analysis/eos_table.txt')
ROT = 3.05e5
NP, NC = 6, 32
RV, _ = csgeom.x1v(128)
TCOLD, FCOLD = 1500.0, 0.30
IREF, RHOREF = 52, 1e-9

JJ, KK = np.meshgrid(np.arange(NC), np.arange(NC), indexing='ij')
DV = np.maximum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK))
DE = np.minimum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK))


def classify(dv, de):
    c = np.full(dv.shape, 2, int)
    c[de <= 2] = 1
    c[dv <= 2] = 0
    return c


CLASS = classify(DV, DE)

GIDMAP = None   # (p, J//16, K//16) -> m


def load(run, n):
    global GIDMAP
    r = bin_convert.read_binary('%s%s/bin/dhj.hydro_w.%05d.bin' % (B, run, n))
    mb = r['mb_data']
    d = np.asarray(mb['dens']).astype(np.float64)
    e = np.asarray(mb['eint']).astype(np.float64)
    v1 = np.asarray(mb['velx']).astype(np.float64)
    g = np.asarray(r['mb_geometry'])
    dens = np.zeros((NP, NC, NC, 128))
    eint = np.zeros_like(dens)
    vr = np.zeros_like(dens)
    gm = np.zeros((NP, 2, 2), int)
    for m in range(24):
        p = m // 4
        j0 = 0 if g[m, 2] < -0.5 else 16
        k0 = 0 if g[m, 4] < -0.5 else 16
        gm[p, j0//16, k0//16] = m
        dens[p, j0:j0+16, k0:k0+16] = np.transpose(d[m], (1, 0, 2))
        eint[p, j0:j0+16, k0:k0+16] = np.transpose(e[m], (1, 0, 2))
        vr[p, j0:j0+16, k0:k0+16] = np.transpose(v1[m], (1, 0, 2))
    GIDMAP = gm
    T, pr = EOS.invert(dens, eint)
    return r['time'], dens, T, vr


def hist(mask):
    c = np.broadcast_to(CLASS[None, :, :], (NP, NC, NC))
    return [int((mask & (c == q)).sum()) for q in range(3)]


ALLH = hist(np.ones((NP, NC, NC), bool))
TOT = sum(ALLH)


def gid_of(p, J, K):
    return int(GIDMAP[p, J//16, K//16])


def local_jk(J, K):
    # ghost-inclusive block indices: active j = 2..17 for a 16-cell block
    return (J % 16) + 2, (K % 16) + 2


log = open(OUT + 'coldmap_s05.txt', 'w')


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    log.write(s+'\n')


P('si/s05 (semi-implicit RT, seed 5) -- cold / drained column maps, coldmap.py method')
P('cold column = >30%% of its 128 radial cells with T < 1500 K')
P('drained column = dens[i=%d] < %.0e   r = %.4e cm' % (IREF, RHOREF, RV[IREF]))
P('ALL columns (6*32*32=6144):  V=%d  E=%d  I=%d' % tuple(ALLH))
P('')
store = {}
P('=== TASK 1: cold / drained column counts per dump (compare coldmap.txt) ===')
P('%-8s %4s %9s %9s   %7s %7s   %9s %9s' %
  ('arm', 'n', 'time', 'rot', 'ncold', 'ndrain', 'minrho', 'meanfrac'))
for n in range(1, 6):
    t, dens, T, vr = load('si/s05', n)
    frac = (T < TCOLD).mean(axis=3)
    dr = dens[:, :, :, IREF] < RHOREF
    store[n] = (t, dens, T, vr, frac, dr)
    P('%-8s %4d %9.4e %9.4f   %7d %7d   %9.3e %9.3f' %
      ('si/s05', n, t, t/ROT, int((frac > FCOLD).sum()), int(dr.sum()),
       dens.min(), frac.mean()))
P('')
P('=== TASK 2: geography (V/E/I) of cold and drained columns, si/s05 ===')
P('%-24s %6s %6s %6s   %6s %6s %6s' % ('set', 'V', 'E', 'I', 'V enr', 'E enr', 'I enr'))
P('%-24s %6d %6d %6d   %6.2f %6.2f %6.2f' %
  ('ALL columns', ALLH[0], ALLH[1], ALLH[2], 1, 1, 1))
for n in range(1, 6):
    t, dens, T, vr, frac, dr = store[n]
    for nm, msk in (('cold', frac > FCOLD), ('drained', dr)):
        h = hist(msk)
        s = sum(h)
        if s == 0:
            P('%-24s %6d %6d %6d   %6s %6s %6s' %
              ('si/s05 rot%d %s' % (n, nm), 0, 0, 0, '-', '-', '-'))
        else:
            P('%-24s %6d %6d %6d   %6.2f %6.2f %6.2f' % ('si/s05 rot%d %s' % (n, nm),
              h[0], h[1], h[2], (h[0]/s)/(ALLH[0]/TOT), (h[1]/s)/(ALLH[1]/TOT),
                (h[2]/s)/(ALLH[2]/TOT)))
P('')
P('=== TASK 5: layer temperature statistics (compare coldmap.txt si/s01, ctl2/s01) ===')
P('%-8s %4s | %9s %9s %9s | %9s %9s' %
  ('arm', 'rot', 'A med', 'A p10', 'A f<1500', 'B med', 'B p10'))
for n in range(1, 6):
    t, dens, T, vr, frac, dr = store[n]
    A = T[:, :, :, 67:103].ravel()
    Bl = T[:, :, :, 29:40].ravel()
    P('%-8s %4d | %9.1f %9.1f %9.5f | %9.1f %9.1f' %
      ('si/s05', n, np.median(A), np.percentile(A, 10), (A < TCOLD).mean(),
       np.median(Bl), np.percentile(Bl, 10)))
P('')
P('=== per-panel mean cold fraction (panel 0 = antistellar) ===')
for n in range(1, 6):
    t, dens, T, vr, frac, dr = store[n]
    P('si/s05 rot%d  %s' % (n, ' '.join('%.4f' % x for x in frac.mean(axis=(1, 2)))))
P('')
P('=== TASK: the MOST DRAINED columns at rot 5 (candidates for the abort column) ===')
t, dens, T, vr, frac, dr = store[5]
r52 = dens[:, :, :, IREF]
idx = np.dstack(np.unravel_index(np.argsort(r52.ravel())[:25], r52.shape))[0]
P('%4s %3s %3s %3s %5s %4s %4s  %10s %10s %10s %8s %10s %6s %5s' %
  ('#', 'p', 'J', 'K', 'class', 'dV', 'dE', 'rho[52]', 'minrho', 'v_r[52]', 'Tmin',
      'T[52]', 'gid', 'j,k'))
for r_, (p, J, K) in enumerate(idx):
    jj, kk = local_jk(J, K)
    P('%4d %3d %3d %3d %5s %4d %4d  %10.3e %10.3e %10.3e %8.1f %10.1f %6d %5s' %
      (r_+1, p, J, K, 'VEI'[CLASS[J, K]], DV[J, K], DE[J, K], r52[p, J, K],
       dens[p, J, K].min(), vr[p, J, K, IREF], T[p, J, K].min(), T[p, J, K, IREF],
       gid_of(p, J, K), '%d,%d' % (jj, kk)))
P('')
P('=== drained-column count per gid at rot 5 (dr mask) ===')
cnt = {}
for p in range(NP):
    for J in range(NC):
        for K in range(NC):
            if dr[p, J, K]:
                cnt[gid_of(p, J, K)] = cnt.get(gid_of(p, J, K), 0) + 1
for gd in sorted(cnt, key=lambda x: -cnt[x]):
    P('  gid %2d : %d drained columns' % (gd, cnt[gd]))
P('')
P('=== gid map (panel, J-half, K-half) -> gid ===')
for p in range(NP):
    for a in range(2):
        for b in range(2):
            P('  panel %d  J %2d-%2d  K %2d-%2d  -> gid %2d' %
              (p, a*16, a*16+15, b*16, b*16+15, GIDMAP[p, a, b]))
log.close()
