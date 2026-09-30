"""T(lat,lon) maps of the WASP-121b 1x/10x rot-300 states at several pressures; secondary hot spot search.
lon = 0 substellar, east (+) = rotation/jet direction (dhjcs convention). Mean of dumps 146-150."""
import sys, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhjcs
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
eos = dhjcs.EOS('/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt')
PB = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
LEV = np.log10(PB) + 6.0
out = {}
for run in ('w1x', 'w10x'):
    acc = None
    for b in range(146, 151):
        d = dhjcs.load('/viper/ptmp2/jinma/w121prod_0929/%s/bin/dhj.hydro_w.%05d.bin' % (run, b), eos, 'cs')
        Tl = dhjcs.level_interp(d['T'], d['lp'], LEV)
        ul = dhjcs.level_interp(d['u'], d['lp'], LEV)
        la, lo, TM = dhjcs.latlon_bin(np.asarray(Tl), d['lat'], d['lon'], 4, 4)
        _, _, UM = dhjcs.latlon_bin(np.asarray(ul), d['lat'], d['lon'], 4, 4)
        acc = (TM, UM) if acc is None else (acc[0] + TM, acc[1] + UM)
    TM, UM = acc[0]/5, acc[1]/5
    out[run] = (la, lo, TM, UM)
    np.savez('maps_%s.npz' % run, lat=la, lon=lo, T=TM, u=UM, pbar=PB)
    print('==', run)
    for k, p in enumerate(PB):
        T = TM[k]; eq = T[np.abs(la) < 20].mean(0)
        # local maxima of the equatorial band profile (periodic), prominence >= 15 K
        n = len(lo); pk = []
        for i in range(n):
            w = [(i + j) % n for j in range(-6, 7)]
            if eq[i] == max(eq[w]):
                pk.append((eq[i], lo[i]))
        pk.sort(reverse=True)
        # 2-D maxima: cells hotter than all neighbours within 12 deg
        mx = []
        for a in range(1, len(la) - 1):
            for i in range(n):
                nb = T[max(a-3, 0):a+4][:, [(i + j) % n for j in range(-3, 4)]]
                if T[a, i] == nb.max():
                    mx.append((T[a, i], la[a], lo[i]))
        mx.sort(reverse=True)
        print('  %.0e bar: eq-band maxima (T K, lon) %s | 2-D maxima (T, lat, lon) %s' % (
            p, ', '.join('%.0f@%+.0f' % q for q in pk[:3]),
            ', '.join('%.0f@(%+.0f,%+.0f)' % q for q in mx[:4])))
fig, ax = plt.subplots(len(PB), 2, figsize=(11, 2.3*len(PB)), constrained_layout=True)
for c, run in enumerate(('w1x', 'w10x')):
    la, lo, TM, UM = out[run]
    for k, p in enumerate(PB):
        im = ax[k, c].pcolormesh(lo, la, TM[k], shading='auto', cmap='inferno')
        ax[k, c].contour(lo, la, UM[k]/1e5, levels=[-4, -2, 2, 4, 6, 8], colors='c', linewidths=0.5)
        ax[k, c].set_title('%s  %.0e bar' % (run, p), fontsize=9); ax[k, c].axvline(0, c='w', lw=0.4, ls=':')
        for x in (-90, 90): ax[k, c].axvline(x, c='w', lw=0.4, ls='--')
        fig.colorbar(im, ax=ax[k, c], shrink=0.8)
fig.savefig('Tmaps_rot300.png', dpi=110)
print('wrote Tmaps_rot300.png')
