"""hstA.py run [t1,t2,...]: user.hst summary vs time: w_mlt, L_top/L0, L_mid/L_in, L_MLT/L_in, KE_int, dt, Picard."""
import sys, re, numpy as np
d = sys.argv[1]
f = d + '/agcar3d.user.hst'; x = np.loadtxt(f); hdr = [l for l in open(f) if l.startswith('#  [')][0]
H = {n: x[:, i] for i, n in enumerate(re.findall(r'\[\d+\]=(\S+)', hdr))}
t = H['time']; L0 = H['L_top'][0]
ts = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else list(np.arange(5e5, 8.1e6, 5e5))
print('%s t_end %.4e rows %d' % (d, t[-1], len(t)))
print('      t     w_mlt  Ltop/L0  Lmid/Lin  LMLT/Lin  Lbot/Lin  KE_int    dt     Picard  Mdot_top')
for tt in ts:
    if tt > t[-1] + 1: break
    i = min(np.searchsorted(t, tt), len(t) - 1)
    print('  %.3e  %.3f  %.4f  %.4f   %.4f   %.4f  %.2e  %6.0f  %.2f  %.2e' % (t[i], H['w_mlt'][i], H['L_top'][i]/L0, H['L_mid'][i]/H['L_in'][i],
          H['L_MLT'][i]/H['L_in'][i], H['L_bot'][i]/H['L_in'][i], H['KE_int'][i], H['dt'][i], H['Picard'][i], H['Mdot_top'][i]))
