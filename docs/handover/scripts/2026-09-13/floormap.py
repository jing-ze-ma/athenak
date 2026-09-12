"""RG_fofc: WHERE are the density-floor hits?  Geography (radius / panel / cube vertex /
panel seam / top) of the low-density cells in the newest hydro_w dumps, plus the top
radial profile at the first and newest dump.  Analysis only: writes only into this dir."""
import sys, os
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert

B   = '/viper/u2/jinma/ATHENAK/bench/RG_fofc/'
OUT = B + 'analysis/'
NP, NC, NR = 6, 32, 480
R0, R1 = 9.6e10, 4.16e12
CPOLY = np.array([-0.196600, +8.811622, -17.210081, +8.883710])
DFLOOR = 1.0e-18
import glob
NEWEST = sorted(int(f[-9:-4]) for f in glob.glob('/viper/u2/jinma/ATHENAK/bench/RG_fofc/bin/rg.hydro_w.*.bin'))[-3:]
NEW = NEWEST[-1]
RGATE  = 6.0e-17          # <hydro>/rad_gate_rho population split quoted in the input
RJOIN  = 3.652820e12      # star / corona join face quoted in the input
RSTAR  = 3.2e12

# ---------------- radial grid (identical map to coordinates.cpp cs_coord1d_1) ----------
def stretch(r, c=CPOLY, r0=R0, r1=R1):
    xi = (np.asarray(r, float) - r0)/(r1 - r0)
    u = xi.copy(); xik = xi.copy()
    for k in range(1, len(c)+1):
        u = u + c[k-1]*xik*(1.0-xi); xik = xik*xi
    return r0 + (r1-r0)*u
XE = stretch(np.linspace(R0, R1, NR+1))
Q  = XE[:-1]/XE[1:]
RV = 0.25*(Q*Q+1.0)/((1.0/3.0)*(Q*Q+Q+1.0))*(XE[1:]+XE[:-1])   # volume centroids
DR = XE[1:]-XE[:-1]
IJOIN = int(np.searchsorted(RV, RJOIN))        # first corona cell

# ---------------- angular geography on the 32x32 panel grid ---------------------------
JJ, KK = np.meshgrid(np.arange(NC), np.arange(NC), indexing='ij')
DE = np.minimum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK))          # dist to a seam
DV = np.maximum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK))          # dist to a vertex
def classify(dv, de, w=2):
    c = np.full(dv.shape, 2, int); c[de <= w] = 1; c[dv <= w] = 0
    return c
CLASS = classify(DV, DE)                       # 0 vertex, 1 seam, 2 interior
CNAME = ('vertex', 'seam', 'interior')
ALLH = np.array([(CLASS == q).sum()*NP for q in range(3)], float)

def hist3(mask3d):                             # mask (6,32,32) -> counts per class
    c = np.broadcast_to(CLASS[None], (NP, NC, NC))
    return np.array([int((mask3d & (c == q)).sum()) for q in range(3)])

# ---------------- loader ---------------------------------------------------------------
def load(n):
    r = bin_convert.read_binary(B + 'bin/rg.hydro_w.%05d.bin' % n)
    mb, g = r['mb_data'], np.asarray(r['mb_geometry'])
    out = {}
    for v in ('dens', 'velx', 'vely', 'velz', 'eint'):
        a = np.asarray(mb[v]).astype(np.float64)          # (24, k, j, i)
        f = np.zeros((NP, NC, NC, NR))
        for m in range(24):
            p  = m // 4
            j0 = 0 if g[m, 2] < -0.5 else 16
            k0 = 0 if g[m, 4] < -0.5 else 16
            f[p, j0:j0+16, k0:k0+16] = np.transpose(a[m], (1, 0, 2))
        out[v] = f
    return r['time'], r['cycle'], out

def main():
    log = open(OUT + 'floormap.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); log.write(s + '\n')

    P('RG_fofc floor geography.  dfloor = %.1e   grid 6 x %d x %d x %d' % (DFLOOR,NC,NC,NR))
    P('r[0] = %.4e  r[%d] = %.4e  join r = %.4e -> first corona cell i = %d (r = %.4e)'
      % (RV[0], NR-1, RV[-1], RJOIN, IJOIN, RV[IJOIN]))
    P('geography: vertex = within 2 cells of a panel corner (cube vertex), seam = within 2'
      ' of a panel edge but not a vertex, interior = rest')
    P('all columns 6*32*32 = %d:  vertex %d  seam %d  interior %d'
      % (ALLH.sum(), ALLH[0], ALLH[1], ALLH[2]))
    P('')

    global S, DUMPS, NEW
    DUMPS = [0] + NEWEST
    S = {}
    for n in DUMPS:
        S[n] = load(n)

    # ---------- 1. counts per dump ----------
    P('=== 1. low-density cell counts per dump (ACTIVE cells only; the event log counts'
      ' ghosts too) ===')
    P('%5s %10s %7s %12s %9s %9s %9s %9s %9s %9s' %
      ('dump', 'time', 'cycle', 'min dens', 'n<1.01f', 'n<2f', 'n<10f', 'n<1e-16',
       'n<6e-17', 'ncell'))
    for n in DUMPS:
        t, cyc, f = S[n]; d = f['dens']
        P('%5d %10.4e %7d %12.5e %9d %9d %9d %9d %9d %9d' %
          (n, t, cyc, d.min(),
           (d < 1.01*DFLOOR).sum(), (d < 2*DFLOOR).sum(), (d < 10*DFLOOR).sum(),
           (d < 1e-16).sum(), (d < RGATE).sum(), d.size))
    P('')

    # ---------- 2. radial distribution of the lowest-density cells ----------
    P('=== 2. radial location of the lowest-density cells, newest dump ===')
    t, cyc, f = S[NEW]; dN = f['dens']
    P('min density per radial index, last 25 cells, and the count below 1e-16 there')
    P('%4s %12s %11s %12s %12s %12s %8s' %
      ('i', 'r[cm]', 'r/Rstar', 'dr[cm]', 'min dens', 'mean dens', 'n<1e-16'))
    for i in range(NR-25, NR):
        col = dN[:, :, :, i]
        P('%4d %12.5e %11.4f %12.4e %12.5e %12.5e %8d' %
          (i, RV[i], RV[i]/RSTAR, DR[i], col.min(), col.mean(), (col < 1e-16).sum()))
    P('')
    P('radial histogram of the n<6e-17 cells (%d total), newest dump:' %
      (dN < RGATE).sum())
    idx = np.where((dN < RGATE).any(axis=(0, 1, 2)))[0]
    if idx.size:
        P('  first such radial index %d (r=%.4e), last %d (r=%.4e)'
          % (idx[0], RV[idx[0]], idx[-1], RV[idx[-1]]))
        P('%6s %12s %10s' % ('i', 'r[cm]', 'count'))
        for i in idx:
            c = (dN[:, :, :, i] < RGATE).sum()
            if c: P('%6d %12.5e %10d' % (i, RV[i], c))
    P('')

    # ---------- 3. panel + angular geography ----------
    P('=== 3. panel and angular geography of the low-density cells ===')
    for n in NEWEST:
        t, cyc, f = S[n]; d = f['dens']
        for thr, nm in ((1e-16, '<1e-16'), (RGATE, '<6e-17')):
            m4 = d < thr
            if not m4.any():
                P('dump %d %s : none' % (n, nm)); continue
            col = m4.any(axis=3)                      # (6,32,32) columns with any such cell
            h = hist3(col); s = h.sum()
            pan = m4.sum(axis=(1, 2, 3))
            P('dump %2d %s  cells %7d   columns %5d   per panel %s' %
              (n, nm, m4.sum(), s, ' '.join('%d:%d' % (p, pan[p]) for p in range(NP))))
            P('          columns by class  vertex %5d  seam %5d  interior %5d   '
              'enrichment %.2f / %.2f / %.2f' %
              (h[0], h[1], h[2], *(h/s/(ALLH/ALLH.sum()))))
    P('')

    # ---------- 4. the minimum-density cell itself ----------
    P('=== 4. the 10 lowest-density cells, newest dump ===')
    t, cyc, f = S[NEW]
    d, vx, vy, vz, ei = (f['dens'], f['velx'], f['vely'], f['velz'], f['eint'])
    flat = d.ravel().argsort()[:10]
    P('%8s %4s %4s %4s %4s %12s %12s %12s %12s %12s %12s' %
      ('rank', 'p', 'J', 'K', 'i', 'r[cm]', 'dens', 'dens/floor', 'v_r[cm/s]',
       'vtan[cm/s]', 'e/rho'))
    for q, ix in enumerate(flat):
        p, J, K, i = np.unravel_index(ix, d.shape)
        P('%8d %4d %4d %4d %4d %12.5e %12.5e %12.4f %12.4e %12.4e %12.5e' %
          (q, p, J, K, i, RV[i], d[p,J,K,i], d[p,J,K,i]/DFLOOR, vx[p,J,K,i],
           np.hypot(vy[p,J,K,i], vz[p,J,K,i]), ei[p,J,K,i]/d[p,J,K,i]))
    P('')

    # ---------- 5. cells that DROPPED most since t=0 ----------
    P('=== 5. where has the density dropped most since dump 0? ===')
    d0 = S[0][2]['dens']
    rat = dN/d0
    P('%4s %12s %12s %12s %12s %10s' %
      ('i', 'r[cm]', 'min ratio', 'mean ratio', 'max ratio', 'n<0.5'))
    for i in list(range(NR-30, NR)):
        c = rat[:, :, :, i]
        P('%4d %12.5e %12.4e %12.4e %12.4e %10d' %
          (i, RV[i], c.min(), c.mean(), c.max(), (c < 0.5).sum()))
    P('')
    drop = rat.min(axis=3) < 0.2
    h = hist3(drop); s = max(h.sum(), 1)
    P('columns whose density fell below 0.2x its t=0 value somewhere: %d  '
      '(vertex %d seam %d interior %d, enrichment %.2f / %.2f / %.2f)' %
      (h.sum(), h[0], h[1], h[2], *(h/s/(ALLH/ALLH.sum()))))
    P('')

    # ---------- 6. top radial profiles ----------
    P('=== 6. horizontally averaged profiles over the last 30 radial cells ===')
    P('(v_r = velx, the cubed-sphere x1 = radial component; e/rho = eint/dens)')
    P('%4s %12s %6s | %12s %12s %11s | %12s %12s %11s' %
      ('i', 'r[cm]', 'r/R*', 'dens t0', 'dens tN', 'tN/t0',
       'v_r t0', 'v_r tN', 'e/rho tN'))
    f0, fN = S[0][2], S[NEW][2]
    for i in range(NR-30, NR):
        a0, aN = f0['dens'][:,:,:,i], fN['dens'][:,:,:,i]
        v0, vN = f0['velx'][:,:,:,i], fN['velx'][:,:,:,i]
        eN = fN['eint'][:,:,:,i]
        P('%4d %12.5e %6.3f | %12.5e %12.5e %11.4f | %12.4e %12.4e %11.4e' %
          (i, RV[i], RV[i]/RSTAR, a0.mean(), aN.mean(), aN.mean()/a0.mean(),
           v0.mean(), vN.mean(), (eN/aN).mean()))
    P('')
    P('mass-flux check (4 pi r^2 <rho v_r>, g/s) at the newest dump:')
    P('%4s %12s %14s %14s' % ('i', 'r[cm]', 'rho*vr mean', '4pi r^2 flux'))
    for i in (IJOIN-40, IJOIN-20, IJOIN-5, IJOIN, IJOIN+10, NR-20, NR-5, NR-1):
        fl = (fN['dens'][:,:,:,i]*fN['velx'][:,:,:,i]).mean()
        P('%4d %12.5e %14.5e %14.5e' % (i, RV[i], fl, 4*np.pi*RV[i]**2*fl))
    log.close()

main()


# ======================================================================================
# PART 7 (appended): can the CUBE-VERTEX corner ghost fill undershoot the density floor?
# bvals_cc.cpp:783 cs_fill_corners_cc extrapolates the ng x ng corner ghost block of a
# TRUE cube vertex QUADRATICALLY (weights 10,-15,6 at the third ghost), which is not
# positivity preserving.  ConsToPrim runs over the full range INCLUDING ghosts
# (hydro_tasks.cpp:688), so a negative extrapolant is counted as a dfloor hit and never
# appears in a dump.  Proxy test: apply the same weights to the three ACTIVE cells running
# inward from each of the 24 panel corners and count the sub-floor results.
# ======================================================================================
def part7():
    log = open(OUT + 'vertexfill.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); log.write(s + '\n')
    P('cube-vertex corner-ghost extrapolation proxy (bvals_cc.cpp cs_fill_corners_cc)')
    P('weights for ghost offset d: ((d+1)(d+2)/2, -d(d+2), d(d+1)/2) on the 3 cells inward')
    P('24 panel corners = 8 cube vertices x 3 panels; ng = 3 so 9 ghost cells each')
    P('')
    for n in [0]+NEWEST:
        t, cyc, f = S[n]; d = f['dens']
        P('--- dump %d  t = %.4e ---' % (n, t))
        tot_neg = tot_sub = 0
        worst = (1e99, None)
        radial = {}
        for p in range(NP):
            for (J0, sj) in ((0, +1), (NC-1, -1)):
                for (K0, sk) in ((0, +1), (NC-1, -1)):
                    aj = d[p, J0:J0+3*sj:sj, K0, :] if sj > 0 else d[p, J0:J0-3:-1, K0, :]
                    ak = d[p, J0, K0:K0+3*sk:sk, :] if sk > 0 else d[p, J0, K0:K0-3:-1, :]
                    for dj in (1, 2, 3):
                        w = np.array([0.5*(dj+1)*(dj+2), -dj*(dj+2), 0.5*dj*(dj+1)])
                        for arr, lbl in ((aj, 'j'), (ak, 'k')):
                            e = w[0]*arr[0] + w[1]*arr[1] + w[2]*arr[2]
                            bad = e < DFLOOR
                            tot_sub += int(bad.sum()); tot_neg += int((e < 0).sum())
                            if e.min() < worst[0]:
                                worst = (e.min(), (p, J0, K0, lbl, dj, int(e.argmin())))
                            for i in np.where(bad)[0]:
                                radial[i] = radial.get(i, 0) + 1
        P('extrapolants below dfloor: %d   negative: %d   (of %d tested)'
          % (tot_sub, tot_neg, 24*3*2*NR))
        P('most negative extrapolant %.4e at (panel,J0,K0,dir,d,i) = %s' % worst)
        if radial:
            ks = sorted(radial)
            P('radial indices with sub-floor extrapolants: i = %d..%d  (r = %.4e..%.4e)'
              % (ks[0], ks[-1], RV[ks[0]], RV[ks[-1]]))
            P('  top 12 by count: ' + ' '.join('i=%d(r=%.3e):%d' % (i, RV[i], radial[i])
              for i in sorted(radial, key=radial.get, reverse=True)[:12]))
        P('')
    log.close()
part7()
# ---- PART 8: is the corona draining or launching? ------------------------------------
def part8():
    log = open(OUT + 'corona.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); log.write(s + '\n')
    # shell volume per column: (4 pi /(6*32*32)) * (r_{i+1}^3 - r_i^3)/3
    vol = (4*np.pi/(NP*NC*NC))*(XE[1:]**3 - XE[:-1]**3)/3.0
    P('corona = cells above the join r = %.4e (i >= %d)' % (RJOIN, IJOIN))
    P('%5s %10s %14s %14s %12s %12s %10s' %
      ('dump', 'time', 'M(r>join)[g]', 'M(r>4.0e12)', '<v_r> i=440',
       '<v_r> i=479', 'f(v_r>0) top'))
    for n in [0]+NEWEST:
        t, cyc, f = S[n]; d, v = f['dens'], f['velx']
        m = (d*vol[None,None,None,:])
        P('%5d %10.4e %14.6e %14.6e %12.4e %12.4e %10.4f' %
          (n, t, m[:,:,:,IJOIN:].sum(), m[:,:,:,455:].sum(),
           v[:,:,:,440].mean(), v[:,:,:,479].mean(),
           (v[:,:,:,450:] > 0).mean()))
    P('')
    P('mean radial velocity profile, corona (every 8th cell):')
    P('%4s %12s' % ('i', 'r[cm]') + ''.join(' %13s' % ('<v_r> d%d' % n) for n in [0]+NEWEST))
    for i in range(IJOIN, NR, 8):
        P('%4d %12.5e' % (i, RV[i]) +
          ''.join(' %13.4e' % S[n][2]['velx'][:,:,:,i].mean() for n in [0]+NEWEST))
    log.close()
part8()
