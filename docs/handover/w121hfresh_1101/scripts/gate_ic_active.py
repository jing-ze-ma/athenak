import sys, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import bin_convert, dclib
eos = dclib.EOS(sys.argv[3])
raw, rho, T, p = dclib.load(sys.argv[1], eos)
tab = np.loadtxt(sys.argv[2])
Tt = np.interp(np.log10(p), tab[:, 0], tab[:, 1]); e = T/Tt - 1
print(e.shape); a = (slice(None), slice(2, -2), slice(2, -2), slice(2, 114)); ea, ra, pa, Ta, Tta = e[a], rho[a], p[a]/1e6, T[a], Tt[a]
ok = ra >= 1e-14
for thr in (0.01, 0.1):
    m = ok & (np.abs(ea) > thr)
    print('|dT/T|>%g: %d cells; p %.2e..%.2e bar; rho %.2e..%.2e' % (thr, m.sum(), pa[m].min() if m.any() else 0, pa[m].max() if m.any() else 0, ra[m].min() if m.any() else 0, ra[m].max() if m.any() else 0))
i = np.unravel_index(np.argmax(np.where(ok, np.abs(ea), 0)), ea.shape)
print('worst', i, 'p %.3e rho %.3e T %.0f table %.0f' % (pa[i], ra[i], Ta[i], Tta[i]))
col = i[:-1]
for k in range(i[-1]-3, min(i[-1]+4, 112)):
    print(k, '%.3e %.3e %.0f %.0f' % (pa[col+(k,)], ra[col+(k,)], Ta[col+(k,)], Tta[col+(k,)]))
print('ACTIVE (x1, x2, x3 ghosts excluded), rho >= 1e-14: %d cells, max |dT/T| %.2e' % (ok.sum(), np.abs(ea[ok]).max()))
for lo, hi in ((1e-9, 1e-6), (1e-6, 1e-3), (1e-3, 1), (1, 300), (300, 3000)):
    m = ok & (pa >= lo) & (pa < hi)
    if m.any(): print('  %7.0e..%-7.0e bar: n %7d max |dT/T| %.2e' % (lo, hi, m.sum(), np.abs(ea[m]).max()))
