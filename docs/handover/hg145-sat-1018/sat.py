"""sat.py: deep-convection saturation check (TASK 2c703c2a) on the DeltaAI hg145 0-28 d dumps.
Per dump (every 2 d + last 3): L_conv,res = Omega r^2 (Fenth+Fke), L_rad = Omega r^2 Fr (shell means from red_*.npz, cgs),
ratio L_conv,res/(L_in - L_rad) averaged over 20-30, 30-60, 60-100 Rsun; rms v_r (npz vrrms, mass weighted) and max|v_r| (bins)
per band 8-19, 19-30, 30-60, 60-100 Rsun; KE_int e-fold time per 2-d window from the hst."""
import os, sys, glob, numpy as np
sys.path.insert(0, os.environ['VISPY']); import bin_convert as bc
RED = '/work/hdd/bivj/jma20/hg145_prod_1018/page/data'; BIN = '/work/hdd/bivj/jma20/hg145_prod_1018/bin'
HST = '/work/nvme/bivj/jma20/hg145_prod_1018/run/hegiant.user.hst'
RS = 6.957e10; OM = np.pi / 3
u = np.loadtxt(HST); _, keep = np.unique(u[::-1, 0], return_index=True); u = u[::-1][keep]   # last occurrence of each time
th, Lin, KE, Lmlt = u[:, 0], u[:, 6], u[:, 10], u[:, 23]
idx = list(range(0, 113, 8)) + [110, 111]
idx = sorted(set(idx))
LB = [(20, 30), (30, 60), (60, 100)]; VB = [(8, 19), (19, 30), (30, 60), (60, 100)]
print('# t[d] | L_conv,res/(L_in-L_rad) at ' + ' '.join('%d-%d' % b for b in LB) + ' | (L_in-L_rad)/L_in same bands | L_MLT/L_in (hst) | '
      'rms v_r / max|v_r| km/s at ' + ' '.join('%d-%d' % b for b in VB) + ' | KE_int erg')
for i in idx:
    d = np.load('%s/red_%05d.npz' % (RED, i)); t = float(d['t']); r = d['r'] / RS
    j = np.argmin(abs(th - t)); L_in = Lin[j]
    Lc = OM * d['r'] ** 2 * (d['Fenth'] + d['Fke']); Lr = OM * d['r'] ** 2 * d['Fr']
    rat, frac = [], []
    for a, b in LB:
        w = (r >= a) & (r < b)
        rat.append(np.mean(Lc[w]) / np.mean(L_in - Lr[w])); frac.append(np.mean(L_in - Lr[w]) / L_in)
    h = bc.read_binary_as_athdf('%s/hegiant.hydro_w.%05d.bin' % (BIN, i)); vr = np.asarray(h['velx'], float) / 1e5
    rr = np.asarray(h['x1v'], float) / RS
    vs = []
    for a, b in VB:
        w = (r >= a) & (r < b); wb = (rr >= a) & (rr < b)
        m = d['rho'][w]; vs.append((np.sqrt(np.sum(m * d['vrrms'][w] ** 2) / np.sum(m)) / 1e5, np.abs(vr[..., wb]).max()))
    print('%5.2f | %s | %s | %.3f | %s | %.2e' % (t / 86400, ' '.join('%+.3f' % x for x in rat), ' '.join('%+.3f' % x for x in frac),
          Lmlt[j] / L_in, '  '.join('%.2f/%.1f' % v for v in vs), KE[j]))
print('# KE_int e-folding time per 2-d window (d)')
for t0 in range(0, 28, 2):
    a, b = np.argmin(abs(th - t0 * 86400)), np.argmin(abs(th - (t0 + 2) * 86400))
    if KE[a] > 0 and KE[b] > 0: print('  %2d-%2d d: %.2f' % (t0, t0 + 2, 2.0 / np.log(KE[b] / KE[a])))
