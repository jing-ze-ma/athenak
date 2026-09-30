"""hotspot.py: masked-hotspot phase-curve test on the w1x rot-300 state (coordinator 09-30).

  python hotspot.py [nproc=12]
Mask: columns with |lat| < 20 deg and 95 < lon < 135 deg (lat = asin(z), lon = atan2(-y,-x),
0 = substellar, east +, dhjcs), cells with p < 1e-3 bar: T -> mean over the reference columns
(140 < lon < 170, |lat_ref - lat| < 3 deg) of T interpolated in log p to the cell's p.  p kept;
rho rescaled at fixed p by (mu_new/mu_old)(T_old/T_new), mu = FastChem mean molecular weight of
the 1x_eq table at (T, p); chemistry and opacity follow from the new (T, p) through the table.
Only the masked columns are re-solved: P_mask = P_full - sum_q P_q(orig) + sum_q P_q(masked).
Outputs out/hotspot.txt, out/hotspot_T1e-4.png, out/hotspot.npz.
"""
import sys
import numpy as np
from multiprocessing import Pool
import rtcore as rc
import emis_prt as ep
import products as pr
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhjcs   # noqa: E402

SD = rc.SD
BANDS = ['NRS1', 'NRS2', 'SOSS1', 'SOSS2', 'WFC3', 'IRAC1', 'IRAC2', 'CO4.5-4.8']


def masked_state():
    st = dict(np.load(SD + 'state/w1x.npz'))
    X = st['X']
    lat = np.degrees(np.arcsin(X[:, 2]))
    lon = np.degrees(np.arctan2(-X[:, 1], -X[:, 0]))
    msk = (abs(lat) < 20) & (lon > 95) & (lon < 135)
    ref = (lon > 140) & (lon < 170)
    t = np.load(SD + 'tables/prt_1x_eq.npz')
    T, p, rho = st['T'].copy(), st['p'], st['rho'].copy()
    lp = np.log10(p)

    def mu_at(Tc, pc):
        i, f = ep.interp(np.log10(Tc), t['lT'])
        j, g = ep.interp(np.log10(pc/1e6), t['lP'])
        M = t['mu']
        return (1-f)*((1-g)*M[i, j] + g*M[i, j+1]) + f*((1-g)*M[i+1, j] + g*M[i+1, j+1])
    for q in np.nonzero(msk)[0]:
        rs = np.nonzero(ref & (abs(lat - lat[q]) < 3))[0]
        cells = p[q] < 1e3                       # 1e-3 bar
        Tn = np.zeros(cells.sum())
        for r in rs:
            Tn += np.interp(lp[q, cells], lp[r, ::-1], T[r, ::-1])
        Tn /= len(rs)
        To = T[q, cells]
        rho[q, cells] *= mu_at(Tn, p[q, cells])/mu_at(To, p[q, cells])*To/Tn
        T[q, cells] = Tn
    st['T'], st['rho'] = T, rho
    np.savez(SD + 'state/w1x_mask.npz', **st)
    return msk, lat, lon, st


def contrib(q, mu):
    """band contribution functions dI/dlayer (norm.) of column q at emission cosine mu."""
    G = ep.G
    lk = G['lk']
    iT, fT = ep.interp(np.log10(G['T'][q]), G['lT'])
    iP, fP = ep.interp(np.log10(G['p'][q]), G['lP'])
    a, b = fT[:, None, None], fP[:, None, None]
    L = ((1 - a)*((1 - b)*lk[iT, iP] + b*lk[iT, iP + 1])
         + a*((1 - b)*lk[iT + 1, iP] + b*lk[iT + 1, iP + 1]))
    rf = G['rf']
    dz = np.diff(rf)
    rcen = 0.5*(rf[1:] + rf[:-1])
    dtau = (10.0**L).astype(float)*(G['rho'][q]*dz)[:, None, None]
    S = (rc.planck_lam(G['lam'][None, :]*1e-4, G['T'][q][:, None])*rcen[:, None]**2)[:, :, None]
    Sf = rc.face_source(S, dtau)
    t = dtau/mu
    e = np.exp(-t)
    aa = np.where(t > 1e-4, -np.expm1(-t)/np.maximum(t, 1e-30), 1 - 0.5*t)
    add = Sf[1:]*(1 - aa) + Sf[:-1]*(aa - e)                 # emitted in layer l
    above = np.cumsum(t[::-1], axis=0)[::-1] - t              # tau from the layer top up
    c = np.einsum('lfg,g->lf', add*np.exp(-above), G['gw'])  # (nlay, nf)
    return c


def run_cols(qs):
    return [(q, ep.column(q)) for q in qs]


if __name__ == '__main__':
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    msk, lat, lon, stm = masked_state()
    qs = list(np.nonzero(msk)[0])
    out = ['# mask: |lat| < 20, 95 < lon < 135 deg, p < 1e-3 bar; %d columns; reference '
           '140 < lon < 170, |dlat| < 3 deg' % len(qs)]
    res = {}
    for arm in ('w1x', 'w1x_mask'):
        ep.init(arm, '1x_eq')
        ch = [qs[i::nproc] for i in range(nproc)]
        with Pool(nproc) as pool:
            parts = pool.map(run_cols, ch)
        nf = len(ep.G['lam'])
        P = np.zeros((72, nf))
        for part in parts:
            for q, (_, F, vis, cc) in part:
                P[vis] += cc.T
        res[arm] = P
        if arm == 'w1x':
            # contribution functions of the knot columns at mu = 0.7 (orig state)
            C = sum(contrib(q, 0.7) for q in qs)
            pmean = np.exp(np.mean(np.log(ep.G['p'][qs]), axis=0))   # bar (ep.G p is in bar)
    full = np.load(SD + 'out/prt_w1x_1x_eq.npz')
    Pf = full['Pang']
    lam, le = full['lam'], full['lam_edges']
    Pm = Pf - res['w1x'] + res['w1x_mask']
    if le[1] < le[0]:
        lam, le, Pf, Pm, C = lam[::-1], le[::-1], Pf[:, ::-1], Pm[:, ::-1], C[:, ::-1]
    dl = np.diff(le)*1e-4
    Fst, _, _ = pr.phoenix_bins(le)
    th = pr.throughputs(lam)
    w = {nm: th[nm]*lam for nm in BANDS[:-1]}
    w['CO4.5-4.8'] = ((lam >= 4.5) & (lam <= 4.8))*lam
    a = pr.full_recipe(Pf, lam, le, w, Fst)
    b = pr.full_recipe(Pm, lam, le, w, Fst)
    phi = np.arange(72)/72
    out.append('%-10s %9s %7s %9s %9s %8s %8s %7s %7s   %s' % (
        'band', 'max|dFp|', 'at phi', 'ampl', 'd ampl', 'off2', 'd off2', 'off1', 'd off1',
        'knot-col contribution fn: p[bar] 16/50/84 %, frac from p<1e-3 bar'))
    for nm in BANDS:
        d = (b[nm] - a[nm])*1e6
        i = np.argmax(abs(d))
        sa, sb = pr.stats(phi, a[nm]), pr.stats(phi, b[nm])
        cw = (C*w[nm]*dl).sum(1)
        cw = cw/cw.sum()
        order = np.argsort(pmean)
        cc = np.cumsum(cw[order])
        pq = [pmean[order][min(np.searchsorted(cc, f), len(cc) - 1)] for f in (0.16, 0.5, 0.84)]
        out.append('%-10s %9.1f %7.3f %9.1f %+9.1f %8.2f %+8.2f %7.2f %+7.2f   %.1e %.1e %.1e  %.3f'
                   % (nm, d[i], phi[i], sa['ampl']*1e6, (sb['ampl'] - sa['ampl'])*1e6, sa['off'],
                      sb['off'] - sa['off'], sa['off1'], sb['off1'] - sa['off1'], *pq,
                      cw[pmean < 1e-3].sum()))
    # T maps at 1e-4 bar before / after
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    st = np.load(SD + 'state/w1x.npz')
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    for k, (nm, TT) in enumerate((('as is', st['T']), ('masked', stm['T']))):
        lp = np.log10(st['p']/1e6)
        Tl = dhjcs.level_interp(TT[:, None, None, :], lp[:, None, None, :], [-4])[0][:, 0, 0]
        la, lo, Mp = dhjcs.latlon_bin(Tl[:, None, None], np.radians(lat)[:, None, None],
                                      np.radians(lon)[:, None, None], 4, 4)
        im = ax[k].pcolormesh(lo, la, Mp, vmin=1000, vmax=4000, cmap='inferno')
        ax[k].set_title('w1x rot 300, T at 1e-4 bar, %s' % nm)
        ax[k].plot([95, 135, 135, 95, 95], [-20, -20, 20, 20, -20], 'c-', lw=1)
        ax[k].set_xlabel('lon (0 = substellar, east +)')
        eq = (abs(la) < 20)
        if k == 0:
            Tb = Mp[eq][:, (lo > 95) & (lo < 135)].max()
        else:
            Ta = Mp[eq][:, (lo > 95) & (lo < 135)].max()
    fig.colorbar(im, ax=ax)
    fig.savefig(SD + 'out/hotspot_T1e-4.png', dpi=110)
    out.append('# T at 1e-4 bar, max over the mask box (4-deg bins): as is %.0f K, masked %.0f K'
               % (Tb, Ta))
    np.savez(SD + 'out/hotspot.npz', phi=phi, **{'asis_' + k: v for k, v in a.items()},
             **{'mask_' + k: v for k, v in b.items()})
    open(SD + 'out/hotspot.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
