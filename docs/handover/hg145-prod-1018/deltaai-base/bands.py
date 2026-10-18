"""bands.py RPH run: per radial band (r/R_ph) shell-mean v_r, rms |v| (mass-weighted) and KE (0.5 rho v^2 dV, spherical polar)
for the dump nearest 2.0 d and the last 3 dumps; then the share of the KE growth since 2.0 d carried by each band."""
import os, sys, glob, numpy as np
sys.path.insert(0, os.environ['VISPY']); import bin_convert as bc
R = float(sys.argv[1]); d = sys.argv[2]
B = [(0.0, 0.2), (0.2, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.01), (1.01, 1.03), (1.03, 1.1), (1.1, 1.5)]
fs = sorted(glob.glob(d + '/bin/hegiant.hydro_w.*.bin'))
def load(f):
    h = bc.read_binary_as_athdf(f)
    r, th, ph = (np.asarray(h[k], float) for k in ('x1f', 'x2f', 'x3f'))
    dV = (np.diff(ph)[:, None, None] * np.diff(-np.cos(th))[None, :, None] * np.diff(r**3 / 3)[None, None, :])
    rho = np.asarray(h['dens'], float); v = [np.asarray(h[k], float) for k in ('velx', 'vely', 'velz')]
    v2 = sum(x**2 for x in v); xr = np.asarray(h['x1v'], float) / R
    out = []
    for a, b in B:
        w = (xr >= a) & (xr < b); m = rho[..., w] * dV[..., w]
        out.append(((v[0][..., w] * dV[..., w]).sum() / dV[..., w].sum() / 1e5,
                    np.sqrt((m * v2[..., w]).sum() / m.sum()) / 1e5, 0.5 * (m * v2[..., w]).sum()))
    return float(h['Time']), out
t = np.array([bc.read_binary_as_athdf(f)['Time'] for f in fs[-25:]]) if False else None
times = {f: None for f in fs}
i2 = min(range(len(fs)), key=lambda i: abs(i * 4320.0 - 172800.0))  # dumps every 4320 s from t 0
sel = [i2] + list(range(len(fs) - 3, len(fs)))
res = [load(fs[i]) for i in sel]
print('### %s  bands r/R_ph: ' % d.split('/')[-1] + ' | '.join('%.2f-%.2f' % ab for ab in B))
for tt, o in res:
    print('  t %.4e (%.2f d)  <v_r> km/s / rms|v| km/s / KE erg:' % (tt, tt / 86400))
    print('    ' + ' | '.join('%+.2f/%.2f/%.2e' % x for x in o))
k0 = np.array([x[2] for x in res[0][1]]); k1 = np.array([x[2] for x in res[-1][1]]); dk = k1 - k0
print('  KE(band) %.2f d -> %.2f d: ' % (res[0][0] / 86400, res[-1][0] / 86400) +
      ' | '.join('%.2e->%.2e (%.0f%% of dKE)' % (a, b, 100 * c / dk.sum()) for a, b, c in zip(k0, k1, dk)))
print('  sum KE bands %.3e -> %.3e' % (k0.sum(), k1.sum()))
