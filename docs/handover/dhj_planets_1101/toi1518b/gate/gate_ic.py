"""IC-fill + hydrostatic gate on a t = 0 hydro_w bin written WITH ghost zones (ng = 2)
(dhj_planets_1101; from kelt20_0925/gate/gate_ic.py: 3x EOS dump, bump terms via
dhj_remap.edges, hydrostatic residual added).
Usage: gate_ic.py <bin> <ic_profile> <gridkeys file> <p_ref bar> <r_pole> <r_eq> <grav> <omega>
- T and p of every cell from (rho, e) through the 3x EOS dump (valid for rho >= 1e-14, the
  dump's first row; the code itself extends its table to eos_logd_min = -19), compared with
  ic_profile T(p) (linear in log p);
- r(p_ref) per column vs planet_setup's pole / equator radii;
- hydrostatic residual of the active cells, |dp/dr + rho g_eff| / (rho g_eff) with
  g_eff = GM/r^2 - omega^2 r sin^2(theta) (r component; centred differences, cell centroids);
- minimum density, density at x1max."""
import sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
sys.path.insert(0, '/viper/ptmp2/jinma/dhj_planets_1101/toi1518b/scripts')
import bin_convert  # noqa: E402,F401
import dclib  # noqa: E402
import dhj_remap as dr_  # noqa: E402

BAR = 1e6
fn, icf, gk = sys.argv[1], sys.argv[2], sys.argv[3]
pref, rpole, req = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
grav, om = float(sys.argv[7]), float(sys.argv[8])
kv = dict(t.split('=') for t in open(gk).read().split())
nx1 = int(kv['mesh/nx1'])
x1min, x1max = float(kv['mesh/x1min']), float(kv['mesh/x1max'])
names = ['f_stretch_r_c%d' % k for k in range(1, 9)] + \
    ['f_stretch_r_b1_amp', 'f_stretch_r_b1_x', 'f_stretch_r_b1_w',
     'f_stretch_r_b2_amp', 'f_stretch_r_b2_x', 'f_stretch_r_b2_w']
c = np.array([float(kv['mesh/' + n]) for n in names] + [0, 0, 0, 0])
import os  # noqa: E402
eos = dclib.EOS(os.environ.get('GATE_EOS_GRID', '/viper/ptmp2/jinma/w121prod_0929/met3/dump/eos_grid.txt'))
raw, rho, T, p = dclib.load(fn, eos)
tab = np.loadtxt(icf)
Ttab = np.interp(np.log10(p), tab[:, 0], tab[:, 1])
err = T/Ttab - 1
ng = 2
ok = rho >= 1e-14
act = err[..., ng:ng + nx1]
oka = ok[..., ng:ng + nx1]
print('cells with rho >= 1e-14 (EOS dump range): %.1f %% of active' % (100*oka.mean()))
print('max |T/T_table - 1|: all cells incl. ghosts %.2e ; active %.2e' %
      (np.abs(err[ok]).max(), np.abs(act[oka]).max()))
pa = p[..., ng:ng + nx1]/BAR
for lo, hi in ((1e-6, 1e-3), (1e-3, 1), (1, 300)):
    m = oka & (pa >= lo) & (pa < hi)
    if m.any():
        print('  active %7.0e..%-7.0e bar: max |dT/T| %.2e' % (lo, hi, np.abs(act[m]).max()))
pin, Tin = p[..., ng]/BAR, T[..., ng]
j = np.unravel_index(np.argmax(pin), pin.shape)
k = np.unravel_index(np.argmin(pin), pin.shape)
print('deepest active cells: pole-like %.1f bar / %.0f K (table %.0f K), equator-like %.1f bar'
      ' / %.0f K (table %.0f K)' % (pin[k], Tin[k], Ttab[..., ng][k], pin[j], Tin[j],
                                    Ttab[..., ng][j]))
pg, Tg = p[..., 0]/BAR, T[..., 0]
m = np.unravel_index(np.argmax(pg), pg.shape)
print('deepest inner ghost: %.1f bar / %.0f K (table %.0f K, rel %.1e)'
      % (pg[m], Tg[m], Ttab[..., 0][m], err[..., 0][m]))
print('table range %.2e .. %.2e bar' % (10**tab[0, 0]/BAR, 10**tab[-1, 0]/BAR))
rf = dr_.edges(nx1, x1min, x1max, c)
rc = dr_.centroids(rf)
lpa = np.log10(pa).reshape(-1, nx1)
rr = np.array([np.interp(np.log10(pref), lq[::-1], rc[::-1]) for lq in lpa])
print('r(p_ref = %.3e bar) over %d columns: %.6e .. %.6e ; tool pole %.6e eq %.6e -> '
      'rel %.1e / %.1e' % (pref, len(rr), rr.min(), rr.max(), rpole, req,
                           rr.min()/rpole - 1, rr.max()/req - 1))
# hydrostatic residual along every column (sin^2 theta from the column's r(p_ref) bulge is
# not needed: bracket it by the pole (0) and the equator (1))
gm = grav*x1min**2
P = p[..., ng:ng + nx1].reshape(-1, nx1)
R = rho[..., ng:ng + nx1].reshape(-1, nx1)
dpdr = (P[:, 2:] - P[:, :-2])/(rc[2:] - rc[:-2])
rm = rc[1:-1]
for s2, lab in ((0.0, 'pole g'), (1.0, 'equator g_eff')):
    ge = gm/rm**2 - om**2*rm*s2
    res = np.abs(dpdr + R[:, 1:-1]*ge)/(R[:, 1:-1]*ge)
    pm = P[:, 1:-1]/BAR
    best = np.minimum(res, 9e9)
    for lo, hi in ((1e-6, 1e-3), (1e-3, 1), (1, 300)):
        mm = (pm >= lo) & (pm < hi)
        print('  hydrostatic |dp/dr + rho g|/(rho g) with %-13s %7.0e..%-7.0e bar: median '
              '%.2e  95%% %.2e' % (lab, lo, hi, np.median(best[mm]), np.quantile(best[mm], 0.95)))
print('min density %.3e g/cc; density in the top active cells %.3e .. %.3e'
      % (rho.min(), rho[..., ng + nx1 - 1].min(), rho[..., ng + nx1 - 1].max()))
