"""summary.py LABEL: key-face table from deepmix_LABEL.npz (wall face, 100/10/1/0.1 bar)."""
import sys
import numpy as np
L = sys.argv[1]
z = np.load('deepmix_%s.npz' % L, allow_pickle=True)
res = z['res'].item()
led = z['ledger']
led = None if led.ndim == 0 else led
C = ['h', 'k', 'phi']
for w in sorted(res):
    r = res[w]
    pf = r['pf']
    S = lambda d: sum(d[k] for k in C)   # noqa: E731
    tot, mmc, se, tr, net = S(r['TOT']), S(r['MMC']), S(r['SE']), S(r['TR']), S(r['NET'])
    print('# %s %s rot %.0f-%.0f (%d restarts)' % (L, w, r['rot'][0], r['rot'][-1], len(r['rot'])))
    print('# face p[bar] TOTAL MMC SE TR | h k phi | cg2/T cg4/T cg8/T | ledger | oe | '
          'K_tot K/(cs dr) K/(vr dr) K/(|v| dxh) vr_rms')
    for pt in [None, 100, 10, 1, 0.1]:
        i = 0 if pt is None else int(np.argmin(np.abs(np.log(pf/pt))))
        lg = led[i] if (led is not None and w == 'w3') else np.nan
        print('%3d %8.3g %+.2e %+.2e %+.2e %+.2e | %+.2e %+.2e %+.2e | %.2f %.2f %.2f | %+.2e | '
              '%.2f | %+.1e %+.1e %+.1e %+.1e %.0f' % (
                  i + 1, pf[i], tot[i], mmc[i], se[i], tr[i], r['TOT']['h'][i],
                  r['TOT']['k'][i], r['TOT']['phi'][i], r['cg'][1, i]/tot[i],
                  r['cg'][2, i]/tot[i], r['cg'][3, i]/tot[i], lg, r['oe'][i],
                  r['Ktot'][i], r['Ktot'][i]/(r['cs'][i]*r['dr'][i]),
                  r['Ktot'][i]/(r['vr'][i]*r['dr'][i]), r['Ktot'][i]/(r['vv'][i]*r['dxh'][i]),
                  r['vr'][i]))
