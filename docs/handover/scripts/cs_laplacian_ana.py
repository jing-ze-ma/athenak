"""Scalar-diffusion check: measured energy source per cell vs -6 kappa (T - T0) r_fac.

usage: python3 ana.py <rundir> [gamma]
The field has no radial dependence, so the radial fluxes vanish exactly and the measured
source is the discrete ANGULAR operator alone; r_fac is the discrete 1/r^2 the
area-over-volume weighting produces, 1.5 (rr^2 - rl^2)/((rr^3 - rl^3) r_c), so the
comparison is purely angular.
"""
import sys, glob, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert

d = sys.argv[1]
gamma = float(sys.argv[2]) if len(sys.argv) > 2 else 1.666667
gm1 = gamma - 1.0
files = sorted(glob.glob(d + '/bin/lap.*_w*.*.bin'))
r0 = bin_convert.read_binary(files[0]); r1 = bin_convert.read_binary(files[-1])
t0, t1 = r0['time'], r1['time']
e0 = np.asarray(r0['mb_data']['eint']); e1 = np.asarray(r1['mb_data']['eint'])
rho = np.asarray(r0['mb_data']['dens'])
src = (e1 - e0) / (t1 - t0)
T = e0 * gm1 / rho
T0 = 1.0    # p0/d0
# units -> kappa_rad in code units at (T0, p0, rho0)
L, M, tt, mu = 1.0e8, 1.0e18, 527.0, 2.3
mH, kB, sig = 1.6726e-24, 1.380649e-16, 5.670374419e-5
v = L / tt; dens = M / L**3; pres = dens * v**2; temp = mu * mH * v**2 / kB
def freedman(Tk, p, met=0.0):
    T1 = min(max(Tk, 75), 4000); p1 = min(max(p, 1), 3e8); lT = np.log10(T1); lp = np.log10(p1)
    c1, c2, c3, c4, c5, c6, c7 = 10.602, 2.882, 6.09e-15, 2.954, -2.526, 0.843, -5.490
    if T1 < 800: c8, c9, c10, c11, c12, c13 = -14.051, 3.055, 0.024, 1.877, -0.445, 0.8321
    else: c8, c9, c10, c11, c12, c13 = 82.241, -55.456, 8.754, 0.7048, -0.0414, 0.8321
    lkl = c1*np.arctan(lT-c2) - c3/(lp+c4)*np.exp((lT-c5)**2) + c6*met + c7
    lkh = c8 + c9*lT + c10*lT**2 + lp*(c11+c12*lT) + c13*met*(0.5+np.arctan((lT-2.5)/0.2)/np.pi)
    return 10**lkl + 10**lkh
import re as _re
_txt=open(d+'/out.txt').read()
_m=_re.search(r'rad_kappa_fac\s*=\s*([0-9.eE+-]+)',_txt)
kfac = float(_m.group(1)) if _m else 0.01
amp = 0.01
_m=_re.search(r'\bamp\s*=\s*([0-9.eE+-]+)',_txt)
amp = float(_m.group(1)) if _m else amp
Tk = T0 * temp; pc = 1.0 * pres; rc = 1.0 * dens
kap_cgs = 16*sig*Tk**3/(3*kfac*freedman(Tk, pc)*rc)
kap_unit = pres * v * L / temp
kap = kap_cgs / kap_unit
# discrete radial factor per cell
g = np.asarray(r0['mb_geometry']); ni = e0.shape[-1]
rf = g[0, 0] + (g[0, 1] - g[0, 0]) * np.arange(ni + 1) / ni
rl, rr = rf[:-1], rf[1:]; rcn = 0.5 * (rl + rr)
rfac = 1.5 * (rr**2 - rl**2) / ((rr**3 - rl**3) * rcn)
expect = -6.0 * kap * (T - T0) * rfac[None, None, None, :]
m = np.abs(expect) > 0.05 * np.abs(expect).max()
slope = np.sum(src * expect) / np.sum(expect * expect)
resid = src - slope * expect
print('%s: t=%.3e (%d dumps)  kappa_code=%.4g  fitted amplitude ratio %.5f  '
      'L1(resid)/L1(expect) = %.4e  Linf = %.4e  (raw L1 err vs expect %.4e)' % (
      d, t1 - t0, len(files), kap, slope,
      np.abs(resid).sum()/np.abs(expect).sum(), np.abs(resid).max()/np.abs(expect).max(),
      np.abs(src - expect).sum()/np.abs(expect).sum()))
# breakdown by distance from the panel edges (cs) -- cells (k,j) within nb of an edge
nmb, nk, nj, ni = e0.shape
kk, jj = np.meshgrid(np.arange(nk), np.arange(nj), indexing='ij')
dedge = np.minimum(np.minimum(kk, nk-1-kk), np.minimum(jj, nj-1-jj))
nedge = (np.minimum(kk, nk-1-kk) < 2).astype(int) + (np.minimum(jj, nj-1-jj) < 2).astype(int)
for name, sel in [('interior (>=2 from edge)', dedge >= 2), ('edge band (1 edge)', nedge == 1), ('vertex band (2 edges)', nedge == 2)]:
    s3 = np.broadcast_to(sel[None, :, :, None], e0.shape)
    print('   %-26s L1(resid)/L1(expect,all) = %.4e   Linf = %.4e  (n=%d)' % (name, np.abs(resid[s3]).sum()/np.abs(expect).sum(), np.abs(resid[s3]).max()/np.abs(expect).max(), sel.sum()))
