"""Loader for the entropy analysis (bundle copy of page/entropy/eload.py, 10-20): one hydro_w dump (+ the m1 dump at the
same time if it exists) on the TRUE stretched grid (common.load: the format-1.2 per-MeshBlock faces mb_x1f/mb_x2f/mb_x3f;
linspace faces are wrong for this poly + log-tail r grid). s in k_B per baryon: s_gas from the run's own EOS table by
thermodynamic integration (sthermo.STable; matches Saha H/He to <= 0.03), T_gas from the table inversion (tgas.py),
s_rad = 4 a T^3/(3 rho) (LTE, T_gas) and from the M1 E_r. tau = radial Rosseland depth from r_out (OPAC table)."""
import os
import numpy as np
from common import load as cload, dvol, tgas, m1_for, OPAC, RHO_U, E_U
from sthermo import STable, s_rad_T, s_rad_E
RSUN = 6.957e10
_ST = None
_KR = None


def tools():
    global _ST, _KR
    if _ST is None:
        _ST = STable()
        if os.path.exists(OPAC):
            hdr = [ln for ln in open(OPAC) if ln.startswith('#')]
            g = [x for x in hdr if x.startswith('# grid:')][0]; vals = hdr[hdr.index(g) + 1][1:].split()
            nT, nD = int(vals[0]), int(vals[1]); lT0, dlT, lD0, dlD = map(float, vals[2:])
            kap = np.loadtxt(OPAC).reshape(nT, nD); _KR = (kap, lT0, dlT, lD0, dlD, nT, nD)
        else:
            print('WARNING no opacity table OPAC=%s: no optical depths' % OPAC)
    return tgas(), _ST, _KR


def kappa_R(rho, T):
    kap, lT0, dlT, lD0, dlD, nT, nD = tools()[2]
    fy = np.clip((np.log10(T) - lT0) / dlT, 0, nT - 1.000001); fx = np.clip((np.log10(rho) - lD0) / dlD, 0, nD - 1.000001)
    j = fy.astype(np.int64); i = fx.astype(np.int64); wy = fy - j; wx = fx - i
    lk = (1 - wy) * ((1 - wx) * kap[j, i] + wx * kap[j, i + 1]) + wy * ((1 - wx) * kap[j + 1, i] + wx * kap[j + 1, i + 1])
    return 10 ** lk


def load(fn, m1=None, tau=True):
    TG, ST, KR = tools()
    h = cload(fn)
    D = dict(t=h['Time'], rf=h['x1f'], r=h['x1v'], thf=h['x2f'], th=h['x2v'], phf=h['x3f'], ph=h['x3v'])
    rho = h['dens'].astype(np.float64); D['rho'] = rho; D['vr'] = h['velx']; D['vth'] = h['vely']; D['vph'] = h['velz']
    D['T'] = 10 ** TG.logT(rho, h['eint'].astype(np.float64))
    rc = rho * RHO_U; D['rho_cgs'] = rc
    D['s'] = ST(rc, D['T']); D['srad'] = s_rad_T(rc, D['T'])
    D['dm'] = rho * dvol(h)
    for k in h:
        if k.startswith('s_'): D[k] = h[k].astype(np.float64)
    if m1 is None:
        m1 = m1_for(D['t']) or ''
    if m1:
        hm = cload(m1, ['m1_e']); assert abs(hm['Time'] - D['t']) < 1e-4, (m1, hm['Time'], D['t'])
        Er = hm['m1_e'].astype(np.float64) * RHO_U * E_U; D['Er'] = Er; D['Trad'] = (np.maximum(Er, 0) / 7.565733e-15) ** 0.25
        D['sradE'] = s_rad_E(rc, Er); D['m1'] = m1
    if tau and KR is not None:
        dr = np.diff(D['rf'])[None, None, :] * RSUN
        dt = kappa_R(rc, D['T']) * rc * dr
        D['tau'] = np.cumsum(dt[:, :, ::-1], axis=2)[:, :, ::-1] - 0.5 * dt
    return D
