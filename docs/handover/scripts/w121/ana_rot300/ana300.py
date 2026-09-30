"""ana300.py: WASP-121b w1x (w121prod_0929) rot-300 analysis.  Pass 1 (all 151 bins, one at a
time): area-weighted isobar T at 0.1/1/10/100 bar and the bottom-cell shell (deepdrift_0927 /
ana_relax.py method: dhjcs EOS inversion, cs area weights, level_interp in log p).
Pass 2 (bins 0 and 150): horizontal-mean profile nabla vs nabla_ad (binprof.py method),
zonal-mean u(lat, p), T at 0.1 bar day/night and hot-spot longitude."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dhjcs   # noqa: E402
import dclib   # noqa: E402
S = '/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.%05d.bin'
O = '/viper/ptmp2/jinma/w121prod_0929/ana_rot300/'
PR = 1.101535e5
EOSF = '/viper/ptmp2/jinma/deepconv_0925/dump/'
eos = dhjcs.EOS(EOSF + 'eos_table.txt')
LEV = np.log10([0.1, 1.0, 10.0, 100.0]) + 6.0


def load(b):
    raw = dhjcs.bin_convert.read_binary(S % b)
    d = np.asarray(raw['mb_data']['dens'], float)
    ei = np.asarray(raw['mb_data']['eint'], float)
    T, p = eos.invert(d, ei)
    G = np.asarray(raw['mb_geometry'])
    nmb, nk, nj = d.shape[:3]
    w = np.empty((nmb, nk, nj))
    for m in range(nmb):
        x = np.tan(np.pi/4*(G[m, 2] + (G[m, 3] - G[m, 2])*(np.arange(nj) + 0.5)/nj))
        y = np.tan(np.pi/4*(G[m, 4] + (G[m, 5] - G[m, 4])*(np.arange(nk) + 0.5)/nk))
        X, Y = np.meshgrid(x, y)
        w[m] = (1 + X*X)*(1 + Y*Y)/(1 + X*X + Y*Y)**1.5
    return raw, d, T, p, w/w.sum()


def wm(q, w, msk=None):
    ww = w*np.isfinite(q)*(1.0 if msk is None else msk)
    return np.nansum(np.where(np.isfinite(q), q, 0)*ww)/ww.sum()


mode = sys.argv[1]
if mode == 'drift':
    rows = []
    for b in range(151):
        raw, d, T, p, w = load(b)
        Tl = dhjcs.level_interp(T, np.log10(p), LEV)
        rows.append([raw['time']/PR] + [wm(Tl[i], w) for i in range(len(LEV))]
                    + [wm(T[..., 0], w), np.exp(wm(np.log(p[..., 0]), w))/1e6])
        del raw, d, T, p
    A = np.array(rows)
    np.save(O + 'drift.npy', A)
    out = ['# rot  T0.1  T1  T10  T100  Tbot  pbot[bar]  (area-weighted isobar means, K)']
    for r in A:
        out.append('%6.1f ' % r[0] + ' '.join('%8.2f' % v for v in r[1:6]) + ' %7.1f' % r[6])
    for lo, hi in ((100, 200), (200, 300), (250, 300)):
        s = (A[:, 0] >= lo - 0.1) & (A[:, 0] <= hi + 0.1)
        line = '# drift rot %d-%d [K/10 rot] (rms):' % (lo, hi)
        for i, n in enumerate(['0.1', '1', '10', '100', 'bot']):
            c = np.polyfit(A[s, 0], A[s, i + 1], 1)
            line += ' %s %+.3f (%.2f)' % (n, 10*c[0], np.std(A[s, i+1] - np.polyval(c, A[s, 0])))
        out.append(line)
    open(O + 'drift.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
else:
    eg = dclib.EOS(EOSF + 'eos_grid.txt')
    PL = np.array([1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100])
    out = []
    for b in (0, 150):
        raw, d, T, p, w = load(b)
        wk = w[..., None]
        gad = eg.at('gad', d, T)
        s = eg.at('s', d, T)
        Tm = (T*wk).sum((0, 1, 2))
        lpm = (np.log(p)*wk).sum((0, 1, 2))
        rho = (d*wk).sum((0, 1, 2))
        sm = (s*d*wk).sum((0, 1, 2))/rho
        ga = (gad*wk).sum((0, 1, 2))
        nab = np.gradient(np.log(Tm), lpm)
        nl = np.gradient(np.log(T), axis=-1)/np.gradient(np.log(p), axis=-1)
        fsup = ((nl > gad)*wk).sum((0, 1, 2))
        out.append('# bin %d rot %.2f: i p[bar] T nabla nab_ad nab-nad s f_superad'
                   % (b, raw['time']/PR))
        for i in range(len(Tm) - 1, -1, -1):
            if np.exp(lpm[i])/1e6 < 0.3:
                continue
            out.append('%3d %9.3f %8.1f %7.4f %7.4f %+8.4f %9.4f %6.3f' % (
                i, np.exp(lpm[i])/1e6, Tm[i], nab[i], ga[i], nab[i] - ga[i], sm[i], fsup[i]))
        if b == 0:
            continue
        geo = dhjcs.cs_geometry(raw)
        u, _ = dhjcs.winds(raw, geo)
        ul = dhjcs.level_interp(u, np.log10(p), np.log10(PL) + 6)
        lat = np.degrees(geo['lat'])
        lon = np.degrees(geo['lon'])
        edges = np.arange(-90, 91, 6)
        out.append('# rot 300 zonal-mean u [m/s] (area-weighted, 6-deg lat bins); cols p[bar]')
        out.append('# lat  ' + ' '.join('%8.0e' % x for x in PL))
        U = np.zeros((len(edges) - 1, len(PL)))
        for il in range(len(edges) - 1):
            mk = (lat >= edges[il]) & (lat < edges[il+1])
            U[il] = [wm(ul[k], w, mk)/100 for k in range(len(PL))]
            out.append('%5.0f ' % (edges[il] + 3) + ' '.join('%8.1f' % x for x in U[il]))
        eq = np.abs(lat) < 20
        ueq = np.array([wm(ul[k], w, eq)/100 for k in range(len(PL))])
        out.append('# equatorial (|lat|<20) mean u: ' + ' '.join('%.0e:%.0f' % (a, x)
                                                                  for a, x in zip(PL, ueq)))
        # finer in p for the jet maximum
        PF = np.logspace(-5, 2.5, 61)
        uf = dhjcs.level_interp(u, np.log10(p), np.log10(PF) + 6)
        ueqf = np.array([wm(uf[k], w, np.abs(lat) < 6) for k in range(len(PF))])/100
        k = np.nanargmax(ueqf)
        out.append('# max |lat|<6 zonal-mean u = %.0f m/s at p = %.3g bar' % (ueqf[k], PF[k]))
        out.append('# max zonal-mean u any lat bin (p grid above) = %.0f m/s' % np.nanmax(U))
        T1 = dhjcs.level_interp(T, np.log10(p), [5.0])[0]
        day = np.abs(lon) < 90
        Td, Tn = wm(T1, w, day), wm(T1, w, ~day)
        lb = np.arange(-180, 181, 6)
        Te = np.array([wm(T1, w, eq & (lon >= lb[i]) & (lon < lb[i+1]))
                       for i in range(len(lb) - 1)])
        lc = lb[:-1] + 3
        a1 = np.sum(Te*np.exp(-1j*np.radians(lc)))
        out.append('# 0.1 bar: dayside mean %.0f K, nightside mean %.0f K, contrast %.0f K; '
                   'equatorial band max %.0f K at lon %+.0f, min %.0f K at lon %+.0f; '
                   'first-harmonic peak lon %+.1f deg (east +)'
                   % (Td, Tn, Td - Tn, Te.max(), lc[Te.argmax()], Te.min(), lc[Te.argmin()],
                      -np.degrees(np.angle(a1))))
        np.savez(O + 'rot300_maps.npz', U=U, PL=PL, PF=PF, ueqf=ueqf, Te=Te, lc=lc)
    open(O + 'struct.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
