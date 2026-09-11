"""Cold-column geography on the cubed sphere: are the cold/drained columns at the
cube vertices, at the panel edges, or everywhere?  ctl2/s01 (explicit RT, dies rot 5.71)
vs si/s01 (semi-implicit, survives).  Analysis only -- writes only into this directory."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/analysis')
import bin_convert
import dhjcs
import csgeom

B = '/viper/u2/jinma/ATHENAK/bench/cs_ens/'
OUT = B + 'analysis/coldmap/'
EOS = dhjcs.EOS(B + 'analysis/eos_table.txt')
ROT = 3.05e5
NP, NC = 6, 32          # 6 panels, 32x32 columns each
RV, _ = csgeom.x1v(128)
TCOLD, FCOLD = 1500.0, 0.30
IREF, RHOREF = 52, 1e-9

# ---- column grid geometry (J = xi index, K = eta index, both 0..31 per panel) ----
JJ, KK = np.meshgrid(np.arange(NC), np.arange(NC), indexing='ij')   # (J,K)
# Chebyshev distance in cells to the nearest of the 4 panel corners = the 4 cube vertices
DV = np.minimum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK)) * 0 + \
    np.maximum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK))
DE = np.minimum(np.minimum(JJ, NC-1-JJ), np.minimum(KK, NC-1-KK)
                )   # dist to nearest panel edge


def classify(dv, de):
    """0 = within 2 cells of a cube vertex; 1 = within 2 of an edge but not a vertex; 2 = interior."""
    c = np.full(dv.shape, 2, int)
    c[de <= 2] = 1
    c[dv <= 2] = 0
    return c


CLASS = classify(DV, DE)     # (32,32) per panel, same for all panels


def load(run, n):
    r = bin_convert.read_binary('%s%s/bin/dhj.hydro_w.%05d.bin' % (B, run, n))
    mb = r['mb_data']
    d = np.asarray(mb['dens']).astype(np.float64)     # (24, k, j, i)
    e = np.asarray(mb['eint']).astype(np.float64)
    g = np.asarray(r['mb_geometry'])
    # regrid meshblocks -> (panel, J, K, i)
    dens = np.zeros((NP, NC, NC, 128))
    eint = np.zeros_like(dens)
    for m in range(24):
        p = m // 4
        j0 = 0 if g[m, 2] < -0.5 else 16
        k0 = 0 if g[m, 4] < -0.5 else 16
        # file index [m,k,j,i] -> panel [J=j, K=k]
        dens[p, j0:j0+16, k0:k0+16] = np.transpose(d[m], (1, 0, 2))
        eint[p, j0:j0+16, k0:k0+16] = np.transpose(e[m], (1, 0, 2))
    T, pr = EOS.invert(dens, eint)
    return r['time']/ROT, dens, T, pr


def stats(dens, T):
    frac = (T < TCOLD).mean(axis=3)          # (6,32,32) cold fraction of the column
    rmin = dens.min(axis=3)
    drained = dens[:, :, :, IREF] < RHOREF
    return frac, rmin, drained


def hist(mask):
    """counts in the 3 geography classes for a (6,32,32) boolean mask."""
    c = np.broadcast_to(CLASS[None, :, :], (NP, NC, NC))
    return [int((mask & (c == q)).sum()) for q in range(3)]


ALLH = hist(np.ones((NP, NC, NC), bool))


def main():
    log = open(OUT + 'coldmap.txt', 'w')

    def P(*a):
        s = ' '.join(str(x) for x in a)
        print(s)
        log.write(s + '\n')

    P('cold column = >%.0f%% of its 128 radial cells with T < %.0f K' % (100*FCOLD,
                                                                         TCOLD))
    P('drained column = dens[i=%d] < %.0e  (healthy ~3e-8)   r(i=%d) = %.4e cm' %
      (IREF, RHOREF, IREF, RV[IREF]))
    P('geography classes (Chebyshev cell distance on the 32x32 panel grid):')
    P('   V = within 2 cells of a cube vertex (panel corner)')
    P('   E = within 2 cells of a panel edge but not a vertex')
    P('   I = panel interior')
    P('ALL columns (6*32*32=6144):  V=%d  E=%d  I=%d' % tuple(ALLH))
    P('')

    store = {}
    P('=== TASK 1: cold / drained column counts per dump ===')
    P('%-6s %4s %8s   %7s %7s   %9s %9s' %
      ('arm', 'rot', 'time/rot', 'ncold', 'ndrain', 'minfrac', 'meanfrac'))
    for run in ('ctl2/s01', 'si/s01'):
        for n in (1, 2, 3, 4, 5):
            rot, dens, T, pr = load(run, n)
            frac, rmin, dr = stats(dens, T)
            store[(run, n)] = (rot, frac, rmin, dr, dens, T)
            cold = frac > FCOLD
            P('%-6s %4d %8.4f   %7d %7d   %9.3e %9.3f' %
              (run.split('/')[0], n, rot, cold.sum(), dr.sum(), rmin.min(), frac.mean()))
    P('')

    P('=== TASK 2: geography of cold and drained columns, ctl2 rot 3,4,5 ===')
    P('%-22s %6s %6s %6s   %8s %8s %8s' %
      ('set', 'V', 'E', 'I', 'V frac', 'E frac', 'I frac'))
    tot = sum(ALLH)
    P('%-22s %6d %6d %6d   %8.4f %8.4f %8.4f' %
      ('ALL columns', ALLH[0], ALLH[1], ALLH[2], ALLH[0]/tot, ALLH[1]/tot, ALLH[2]/tot))
    for run in ('ctl2/s01', 'si/s01'):
        for n in (1, 2, 3, 4, 5):
            rot, frac, rmin, dr, dens, T = store[(run, n)]
            for nm, msk in (('cold', frac > FCOLD), ('drained', dr)):
                h = hist(msk)
                s = sum(h)
                if s == 0:
                    P('%-22s %6d %6d %6d   %8s %8s %8s' %
                      ('%s rot%d %s' % (run.split('/')[0], n, nm), 0, 0, 0, '-', '-',
                          '-'))
                else:
                    P('%-22s %6d %6d %6d   %8.4f %8.4f %8.4f' %
                      ('%s rot%d %s' % (run.split('/')[0], n, nm), h[0], h[1], h[2],
                       h[0]/s, h[1]/s, h[2]/s))
    P('')
    P('enrichment = (fraction of the set in class) / (fraction of all columns in class)')
    P('%-22s %8s %8s %8s' % ('set', 'V', 'E', 'I'))
    for run in ('ctl2/s01', 'si/s01'):
        for n in (1, 2, 3, 4, 5):
            rot, frac, rmin, dr, dens, T = store[(run, n)]
            for nm, msk in (('cold', frac > FCOLD), ('drained', dr)):
                h = hist(msk)
                s = sum(h)
                if s == 0:
                    continue
                P('%-22s %8.2f %8.2f %8.2f' % ('%s rot%d %s' % (run.split('/')[0], n, nm),
                  (h[0]/s)/(ALLH[0]/tot), (h[1]/s)/(ALLH[1]/tot), (h[2]/s)/(ALLH[2]/tot)))
    P('')

    P('=== TASK 3: ctl2 rot 2, the FIRST cold cells (T<1500 K at i>60) ===')
    rot, frac, rmin, dr, dens, T = store[('ctl2/s01', 2)]
    fup = (T[:, :, :, 61:] < TCOLD).mean(axis=3)     # cold fraction above i=60
    ncell = (T[:, :, :, 61:] < TCOLD).sum()
    P('rot %.4f: %d cells with T<1500 K at i>60, in %d distinct columns' %
      (rot, ncell, int((fup > 0).sum())))
    # tie-break the many frac=1.0 columns by the whole-column cold fraction
    key = fup.ravel() + 1e-3*frac.ravel()
    idx = np.dstack(np.unravel_index(np.argsort(-key)[:10], fup.shape))[0]
    P('%4s %3s %3s %3s  %9s %9s %6s %6s %5s  %9s %6s' %
      ('#', 'p', 'J', 'K', 'frac>i60', 'frac_all', 'd_vert', 'd_edge', 'class',
          'rho[i=52]', 'imin'))
    for r_, (p, J, K) in enumerate(idx):
        cold_i = np.where(T[p, J, K] < TCOLD)[0]
        P('%4d %3d %3d %3d  %9.4f %9.4f %6d %6d %5s  %9.3e %6d' %
          (r_+1, p, J, K, fup[p, J, K], frac[p, J, K], DV[J, K], DE[J, K],
           'VEI'[CLASS[J, K]], dens[p, J, K, IREF], cold_i.min() if len(cold_i) else -1))
    hs = hist(fup >= 1.0)
    ss = max(sum(hs), 1)
    P('columns FULLY cold above i=60 (frac>i60 == 1): n=%d  V=%d E=%d I=%d'
      '  (enrichment %.2f %.2f %.2f)' %
      (int((fup >= 1.0).sum()), hs[0], hs[1], hs[2],
       (hs[0]/ss)/(ALLH[0]/tot), (hs[1]/ss)/(ALLH[1]/tot), (hs[2]/ss)/(ALLH[2]/tot)))
    # class breakdown of every column that has ANY cold cell above i=60
    h = hist(fup > 0)
    s = max(sum(h), 1)
    P('columns with ANY cold cell at i>60: V=%d E=%d I=%d  (enrichment %.2f %.2f %.2f)' %
      (h[0], h[1], h[2], (h[0]/s)/(ALLH[0]/tot), (h[1]/s)/(ALLH[1]/tot),
          (h[2]/s)/(ALLH[2]/tot)))

    P('')
    P('=== TASK 5: layer-mean temperature bias, explicit vs semi-implicit ===')
    P('layer A = array index 67..102 (ghost-inclusive i 69..104)   r = %.4e .. %.4e cm' %
      (RV[67], RV[102]))
    P('layer B = array index 29..39  (ghost-inclusive i 31..41)    r = %.4e .. %.4e cm' %
      (RV[29], RV[39]))
    P('per-cell statistics over all 6144 columns x the layer cells')
    P('%-6s %4s | %9s %9s %9s %9s | %9s %9s' %
      ('arm', 'rot', 'A med', 'A p10', 'A f<1500', 'A f<=201', 'B med', 'B p10'))
    for run in ('ctl2/s01', 'si/s01'):
        for n in (1, 2, 3, 4, 5):
            rot, frac, rmin, dr, dens, T = store[(run, n)]
            A = T[:, :, :, 67:103].ravel()
            Bl = T[:, :, :, 29:40].ravel()
            P('%-6s %4d | %9.1f %9.1f %9.5f %9.5f | %9.1f %9.1f' %
              (run.split('/')[0], n, np.median(A), np.percentile(A, 10),
               (A < TCOLD).mean(), (A <= 201.0).mean(), np.median(Bl), np.percentile(Bl,
                                                                                     10)))
    P('')
    P('difference ctl2 - si (same rotation)')
    P('%4s | %9s %9s %9s | %9s %9s' % ('rot', 'A med', 'A p10', 'A f<1500', 'B med',
                                       'B p10'))
    for n in (1, 2, 3, 4, 5):
        Tc = store[('ctl2/s01', n)][5]
        Ts = store[('si/s01', n)][5]
        Ac, As = Tc[:, :, :, 67:103].ravel(), Ts[:, :, :, 67:103].ravel()
        Bc, Bs = Tc[:, :, :, 29:40].ravel(), Ts[:, :, :, 29:40].ravel()
        P('%4d | %9.1f %9.1f %9.5f | %9.1f %9.1f' %
          (n, np.median(Ac)-np.median(As), np.percentile(Ac, 10)-np.percentile(As, 10),
           (Ac < TCOLD).mean()-(As < TCOLD).mean(), np.median(Bc)-np.median(Bs),
           np.percentile(Bc, 10)-np.percentile(Bs, 10)))
    log.close()
    np.savez(OUT + 'coldfrac.npz',
             **{'%s_%d' % (r.split('/')[0], n): store[(r, n)][1]
                for r in ('ctl2/s01', 'si/s01') for n in range(1, 6)},
             **{'rot_%s_%d' % (r.split('/')[0], n): store[(r, n)][0]
                for r in ('ctl2/s01', 'si/s01') for n in range(1, 6)},
             cls=CLASS)


if __name__ == '__main__':
    main()
