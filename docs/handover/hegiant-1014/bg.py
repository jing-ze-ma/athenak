"""bg.py RPH R0SPONGE run [every]: background check over the full domain.
 per dump: per radial band <v_r> (shell-mean, km/s) and max|rho/rho0-1| (shell means), mass flux through 1.0 R_ph,
 sponge-edge jumps (v_r and eint across r0), floor cells; hst: L_top/L0, Mdot_top, Min_top, P_sponge/L_top,
 (L_bot-L_in)/L_in, v1sq_wall, Mdot_bot, Picard."""
import os, sys, glob, re, numpy as np
sys.path.insert(0, os.environ.get('VISPY', 'vis/python'))
import bin_convert as bc
import os; BASE = os.environ.get('BASE', 'hegiant'); LOGF = os.environ.get('LOGF', 'run.log')
R = float(sys.argv[1]); r0 = float(sys.argv[2]); d = sys.argv[3]; ev = int(sys.argv[4]) if len(sys.argv) > 4 else 50
f = d + '/' + BASE + '.user.hst'; x = np.loadtxt(f); hdr = [l for l in open(f) if l.startswith('#  [')][0]
H = {n: x[:, i] for i, n in enumerate(re.findall(r'\[\d+\]=(\S+)', hdr))}
log = ''.join(open(x).read() for x in sorted(glob.glob(d + '/' + LOGF)))
m = re.search(r'Picard iterations mean=(\S+) max=(\S+) NON-CONVERGED=(\S+)', log)
fl = re.findall(r'fired cell-stages: at-floor=(\d+) other=(\d+)', log)
print('#### %s  rows %d  t_end %.4g  Picard %s  fofc at-floor max %s other max %s' % (d.split('/')[-1], len(H['time']), H['time'][-1],
      ('%.2f/%d NONCONV %d' % (float(m.group(1)), float(m.group(2)), float(m.group(3)))) if m else 'running',
      max([int(a) for a, b in fl]) if fl else '-', max([int(b) for a, b in fl]) if fl else '-'))
n = len(H['time']); idx = sorted(set([1, n//8, n//4, n//2, 3*n//4, n-1]))
L0 = H['L_top'][0]
for i in idx:
    print('  hst t %.3e L_top/L0 %.5f Mdot_top %.3g Min_top %.3g P_sponge/L %.2e (L_bot-L_in)/L_in %.2e v1sq_wall %.2e Mdot_bot %.2g' % (
        H['time'][i], H['L_top'][i]/L0, H['Mdot_top'][i], H['Min_top'][i], H['P_sponge'][i]/H['L_top'][i],
        (H['L_bot'][i]-H['L_in'][i])/H['L_in'][i], H['v1sq_wall'][i], H['Mdot_bot'][i]))
ts = H['time']; Isp = np.sum(0.5*(H['P_sponge'][1:]+H['P_sponge'][:-1])*np.diff(ts)); IL = np.sum(0.5*(H['L_top'][1:]+H['L_top'][:-1])*np.diff(ts))
print('  int P_sponge dt / int L_top dt = %.2e' % (Isp/IL))
fs = sorted(set(glob.glob(d + '/bin/' + BASE + '.hydro_w.*.bin')))
h0 = bc.read_binary_as_athdf(fs[0]); r = np.asarray(h0['x1v'], float); xf = np.asarray(h0['x1f'], float)
rho0 = np.asarray(h0['dens'], float).mean(axis=(0, 1)); e0 = np.asarray(h0['eint'], float).mean(axis=(0, 1))
xr = r/R; rtop = xf[-1]/R
bands = [(0.0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.01), (1.01, 1.03), (1.03, 1.1), (1.1, 1.5), (1.5, r0/R), (r0/R, 9)]
print('  bands (r/R_ph): ' + ' | '.join('%.2f-%.2f' % (a, min(b, rtop)) for a, b in bands))
ip = np.argmin(np.abs(xr - 1.0)); i0 = min(np.searchsorted(r, r0), len(r) - 1)
seen = set()
for f in fs:
    h = bc.read_binary_as_athdf(f); c = int(h['NumCycles'])
    if c in seen or (c % ev and f != fs[-1]): continue
    seen.add(c)
    rho = np.asarray(h['dens'], float).mean(axis=(0, 1)); vr = np.asarray(h['velx'], float).mean(axis=(0, 1))/1e5
    e = np.asarray(h['eint'], float).mean(axis=(0, 1)); dmin = np.asarray(h['dens'], float).min()
    s = []
    for a, b in bands:
        w = (xr >= a) & (xr < b)
        if not w.any(): s.append('    -      - '); continue
        k = np.argmax(np.abs(vr[w])); dr = np.abs(rho[w]/rho0[w] - 1)
        s.append('%+7.2f %6.3f' % (vr[w][k], dr.max()))
    jump_v = vr[i0] - vr[i0-1]; jump_e = e[i0]/e[i0-1] - e0[i0]/e0[i0-1]
    mph = rho[ip]*vr[ip]*1e5*r[ip]**2
    print('  c%-5d t %.3e  ' % (c, float(h['Time'])) + ' | '.join(s) + '  | Mdot(R_ph)/sr %+.2e  sponge-edge dv %+.2f km/s d(e ratio) %+.1e  min rho %.1e' % (mph, jump_v, jump_e, dmin))
