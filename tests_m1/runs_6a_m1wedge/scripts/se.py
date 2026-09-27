# usage: se.py <run dir> [window s]: steady-state diagnostics of the super-Eddington wedge
# hst: Mdot_top, Mdot_bot, M_tot, L_top/L_in in consecutive windows; bins: horizontal-mean
# profiles rho, T_rad = (E/a)^1/4, v_r, Gamma = -(d(E/3)/dr)/(rho g), tau (kappa from Gamma),
# drift between the last dump and the one a window earlier.
import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/sprhd_0926/wt/vis/python')
import bin_convert as bc
d = sys.argv[1]; win = float(sys.argv[2]) if len(sys.argv) > 2 else 2000.0
GM = 4.180568e26; c = 2.99792458e10; arad = 7.5657e-15; R0 = 3.240543e10
f = glob.glob(d + '/*.user.hst')[0]
hdr = [l for l in open(f) if l.startswith('#')][-1]
names = [s.split('=')[1] for s in hdr[1:].split()]
a = np.atleast_2d(np.loadtxt(f)); h = {n: a[:, i] for i, n in enumerate(names)}
t = h['time']; M0 = h['M_tot'][0]
print('%s: t_end %.0f s, dt min/last %.4f/%.4f' % (d, t[-1], h['dt'].min(), h['dt'][-1]))
print('  window        Mdot_top    Mdot_bot  (g/s, face)   dM/M0      L_top/L_in')
edges = np.arange(0, t[-1] + 1, win)
for t0 in edges[:-1]:
    m = (t >= t0) & (t < t0 + win)
    if m.sum() < 2: continue
    print('  %6.0f-%6.0f %+.4e %+.4e  %+.4e  %.5f+-%.1e' % (t0, t0 + win, h['Mdot_top'][m].mean(),
          h['Mdot_bot'][m].mean(), h['M_tot'][m][-1]/M0 - 1, (h['L_top']/h['L_in'])[m].mean(),
          (h['L_top']/h['L_in'])[m].std()))
fs = sorted(glob.glob(d + '/bin/*.hydro_w.*.bin'))


def prof(fn):
    H = bc.read_binary(fn); Mb = bc.read_binary(fn.replace('hydro_w', 'm1'))
    n1 = H['nx1_out_mb']; geo = np.array(H['mb_geometry'])
    r0, r1 = geo[:, 0].min(), geo[:, 1].max(); nb1 = len(set(geo[:, 0]))
    N = nb1*n1; r = r0 + (np.arange(N) + 0.5)*(r1 - r0)/N
    acc = {k: np.zeros(N) for k in ['rho', 'mom', 'E', 'F']}; cnt = np.zeros(N)
    for b in range(H['n_mbs']):
        i0 = sorted(set(geo[:, 0])).index(geo[b, 0])*n1; s = slice(i0, i0 + n1)
        dd = H['mb_data']['dens'][b]
        acc['rho'][s] += dd.mean(axis=(0, 1)); acc['mom'][s] += (dd*H['mb_data']['velx'][b]).mean(axis=(0, 1))
        acc['E'][s] += Mb['mb_data']['m1_e'][b].mean(axis=(0, 1)); acc['F'][s] += Mb['mb_data']['m1_f1'][b].mean(axis=(0, 1))
        cnt[s] += 1
    p = {k: v/cnt for k, v in acc.items()}
    p['vr'] = p['mom']/p['rho']; p['T'] = (p['E']/arad)**0.25
    p['Gam'] = -np.gradient(p['E']/3, r)/(p['rho']*GM/r**2)
    kap = np.maximum(p['Gam'], 0)*c*(GM/r**2)/np.maximum(p['F'], 1e-30)
    kr = kap*p['rho']; p['tau'] = np.concatenate([np.cumsum((0.5*(kr[1:] + kr[:-1])*np.diff(r))[::-1])[::-1], [0]])
    return H['time'], r, p


tl, r, pl = prof(fs[-1])
ie = max(0, len(fs) - 1 - int(round(win/200.0)))
te, _, pe = prof(fs[ie])
print('  profile drift %.0f -> %.0f s (rel. L1 over the domain): ' % (te, tl) +
      ' '.join('%s %.2e' % (k, np.mean(abs(pl[k] - pe[k]))/np.mean(abs(pl[k]))) for k in ['rho', 'T', 'vr', 'Gam']))
iph = np.argmin(abs(pl['tau'] - 1.0)); sup = np.where(pl['Gam'] > 1)[0]
cs = np.sqrt(5/3*(2/3)*1.0)  # placeholder, not used
print('  photosphere (tau=1) z = %.3e cm; Gamma>1 cells: %s' % (r[iph] - R0,
      ('z %.3e..%.3e (tau %.1f..%.1f)' % (r[sup[0]] - R0, r[sup[-1]] - R0, pl['tau'][sup[0]], pl['tau'][sup[-1]])) if len(sup) else 'none'))
sel = np.linspace(0, len(r) - 1, 8).astype(int)
print('  z[1e7 cm] ' + ' '.join('%8.2f' % ((r[i] - R0)/1e7) for i in sel))
for k, fmt in (('rho', '%8.1e'), ('T', '%8.2e'), ('vr', '%8.1e'), ('Gam', '%8.3f'), ('tau', '%8.1e')):
    print('  %-9s ' % k + ' '.join(fmt % pl[k][i] for i in sel))
