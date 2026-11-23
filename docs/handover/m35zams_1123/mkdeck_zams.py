#!/usr/bin/env python3
"""IA9_m35zams.athinput from Caltech's as-run IA9_m35mms deck: only the star and the grid change.
Grid = S22 T42L5.0 (384 log r, 6.798-8.200 Rsun) with the lateral wedge doubled in both angles
(USER 11-23: 2x lateral extent, 256 x 256 cells = the paper's 5.0e-4 rad cells)."""
import math, re, sys
import numpy as np
src, out = sys.argv[1], sys.argv[2]
ic = '/viper/ptmp2/jinma/m35zams_1115/ic/ic_m35zams_ideal.txt'
CAL = '/resnick/groups/carnegie_poc/jingze/fresh_prod_1120'

mu_mms, tu_mms = 0.6215, 7.4749268633e-09
mu = 0.6220
tu = tu_mms*mu/mu_mms                      # T_unit = mu m_u / k_B
a_cgs = 7.5657332e-15
L, GM = 6.618281e38, 4.6449357021e27
x1min, x1max = 4.7294e11, 5.7031e11        # GRID_REFERENCE (S22 T42L5.0)
f_str = -math.log(x1max/x1min)
th_a, th_b, ph_b = 1.53875208172828071, 1.60284057186151241, 0.064   # S22 wedge
half = 2.0*(math.pi/2 - th_a)
x2min, x2max, x3max = math.pi/2 - half, math.pi/2 + half, 2.0*ph_b
d = np.loadtxt(ic, comments='#'); r, rho, T = d[:, 0], d[:, 1], d[:, 5]
# WB anchors: M35MMS has he_ic_balance_rmax at the IC radius where T = 2.137e5 K and wb_rmax = 1.25x that
# (the same 1.25 ratio as T19L6.4); same rule here.
j = np.argmin(abs(T - 2.137e5)); rbal = r[j]; wbr = 1.25*rbal
rfl = r[np.argmax(rho <= 1.0001e-17)]
fin = L/(4*math.pi*x1min**2)
tfl = 1.0e3/tu
vesc = math.sqrt(2*GM/x1min)
dlnr = -f_str/384; dth = (x2max - x2min)/256; dph = x3max/256

sub = {   # (block, key): (value, comment)
 ('job','basename'): ('m35zams', 'M35ZAMS = T42L5.0 (Schultz+2022)'),
 ('mesh','nx1'): ('384', '[PAPER] S22 sect 2.1 T42L5.0: 384 logarithmically spaced radial bins'),
 ('mesh','x1min'): ('%.10e' % x1min, '[PAPER] S22: 6.798 Rsun'),
 ('mesh','x1max'): ('%.10e' % x1max, '[PAPER] S22: 8.200 Rsun (1.0378 R_ph)'),
 ('mesh','nx2'): ('256', '[USER 11-23] 2 x S22 128: theta wedge doubled at the paper cell size'),
 ('mesh','x2min'): ('%.17f' % x2min, '[USER 11-23] pi/2 - 2 x S22 half-width (S22: 0.4898 pi)'),
 ('mesh','x2max'): ('%.17f' % x2max, '[USER 11-23] pi/2 + 2 x S22 half-width (S22: 0.5102 pi)'),
 ('mesh','nx3'): ('256', '[USER 11-23] 2 x S22 128'),
 ('mesh','x3max'): ('%.17f' % x3max, '[USER 11-23] 2 x S22 0.064'),
 ('mesh','f_stretch_r'): ('%.15f' % f_str, '[PAPER] log bins: = -ln(x1max/x1min) (dlnr = %.6e; dtheta = %.6e, dphi = %.6e)' % (dlnr, dth, dph)),
 ('meshblock','nx1'): ('384', 'vet_gd: one x1 block'),
 ('meshblock','nx2'): ('128', '4 blocks = 1 per GPU on a 4-GPU node; 128 >= vet_gd near-tangent reach 127'),
 ('meshblock','nx3'): ('128', ''),
 ('units','mu'): ('%.4f' % mu, '[IC] ic_m35zams_ideal.txt (build_ic_ideal.py, total-pressure HSE fit); NOT stated in S22'),
 ('hydro','dfloor'): ('1.0e-17', '[NS] = the IC density floor above the photosphere (r > %.3f Rsun)' % (rfl/6.957e10)),
 ('hydro','vceil'): ('3.0e8', '[NS] > escape speed at the base (%.3g)' % vesc),
 ('hydro','tfloor'): (repr(tfl), '1000 K / temp_unit_kelvin (as M35MMS)'),
 ('hydro','wb_rmax'): ('%.6e' % wbr, '1.25 x he_ic_balance_rmax (M35MMS / T19L6.4 FRESHPREP ratio); > x1max'),
 ('rad_m1','arad'): ('%.10e' % (a_cgs*tu**4), 'a_cgs (7.5657332e-15) T_unit^4, T_unit = mu m_u/k_B = %.10e K' % tu),
 ('rad_m1','temp_unit_kelvin'): ('%.10e' % tu, 'mu %.4f' % mu),
 ('rad_m1','implicit_flux_x1min'): ('%.10e' % fin, '[IC] L/(4 pi r_in^2), L = MESA L 6.618281e38'),
 ('problem','he_ic_file'): (CAL + '/m35zams/ic_m35zams_ideal.txt',
     '[IC] IDEAL radiative-only HSE+thermal-eq. envelope matched to the MESA ZAMS surface row (build_ic_ideal.py), F_r = L/(4 pi r^2)'),
 ('problem','he_gm'): ('4.6449357021e+27', '[PAPER] 35 Msun point mass (G 6.6743e-8, Msun 1.98841e33)'),
 ('problem','he_ic_balance_rmax'): ('%.6e' % rbal, 'IC radius at T = 2.137e5 K (same anchor temperature as M35MMS 8.776579e11)'),
}
blk, seen, lines = None, set(), []
for ln in open(src):
    m = re.match(r'^<(\w+)>', ln)
    if m: blk = m.group(1)
    k = re.match(r'^(\w+)\s*=', ln)
    if k and (blk, k.group(1)) in sub:
        key = k.group(1); v, c = sub[(blk, key)]; seen.add((blk, key))
        ln = '%-24s = %-22s # %s\n' % (key, v, c) if c else '%-24s = %s\n' % (key, v)
    lines.append(ln)
miss = set(sub) - seen
if miss: sys.exit('keys not found: %s' % sorted(miss))
hdr = ['# M35ZAMS (= T42L5.0, Schultz, Bildsten & Jiang 2022 ApJL 924 L11) fresh IA9 deck, made 11-23 by\n',
       '#   /viper/ptmp2/jinma/m35zams_1123/mkdeck_zams.py from Caltech\'s as-run IA9_m35mms.athinput (md5 52da1397):\n',
       '#   every scheme / numerics key identical; star + grid keys replaced (S22 T42L5.0 grid, lateral wedge doubled\n',
       '#   to 256 x 256 at the paper cell size, USER 11-23); ideal gas, IC fresh_prod_1120/m35zams/ic_m35zams_ideal.txt.\n']
open(out, 'w').writelines(hdr + lines)
print('T_unit', tu, 'arad', a_cgs*tu**4, 'tfloor', tfl, 'Fin', fin, 'rbal', rbal, 'wb_rmax', wbr, 'r_floor', rfl,
      'x2', x2min, x2max, 'x3max', x3max, 'dlnr', dlnr, 'dth', dth, 'dph', dph, 'cells', 384*256*256)
