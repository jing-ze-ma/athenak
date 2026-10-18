#!/usr/bin/env python3
"""HGPAGE 10-18: per-dump reduction of the hg145 production run (READ ONLY on the run).
usage: reduce.py <run_dir> <out_dir>   writes out_dir/red_NNNNN.npz for every hydro_w+m1 dump pair not yet reduced.
Contents: shell means, photosphere map R_ph(theta,phi) and T_eff map, equatorial + meridional slices, the 4 wedge side faces."""
import os
import sys
import glob
import time
import numpy as np
sys.path.insert(0, os.environ.get('VISPY', '/raven/u/jinma/ATHENAK/wt-accbig/vis/python'))  # AthenaK vis/python (bin_convert)
import bin_convert as bc  # noqa: E402

RUN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
A = 7.565733e-15
SIG = 5.670374e-5
RS = 6.957e10
RMAX_SL = 1.6       # slices/faces kept to r < 1.6 R_ph(IC)
RPH0 = 7.94653e12

fh = sorted(glob.glob(os.path.join(os.environ.get('BINDIR') or os.path.join(RUN, 'bin'), '*.hydro_w.*.bin')))
if os.environ.get('ONLY'):
    fh = [x for x in fh if x.split('.')[-2] in os.environ['ONLY'].split(',')]
for f in fh:
    idx = f.split('.')[-2]
    fm = f.replace('hydro_w', 'm1')
    o = os.path.join(OUT, 'red_%s.npz' % idx)
    if os.path.exists(o) or not os.path.exists(fm):
        continue
    # skip a dump that may still be being written (younger than 2 min)
    if time.time() - max(os.path.getmtime(f), os.path.getmtime(fm)) < 120:
        continue
    try:
        h = bc.read_binary_as_athdf(f)
        m = bc.read_binary_as_athdf(fm)
    except Exception as e:  # incomplete file
        print('skip', idx, e)
        continue
    r = np.asarray(h['x1v'], float); th = np.asarray(h['x2v'], float); ph = np.asarray(h['x3v'], float)
    rf = np.asarray(h['x1f'], float)
    t = float(h['Time'])
    rho = np.asarray(h['dens'], float); vr = np.asarray(h['velx'], float); vt = np.asarray(h['vely'], float)
    vp = np.asarray(h['velz'], float); ei = np.asarray(h['eint'], float)
    E = np.asarray(m['m1_e'], float); F1 = np.asarray(m['m1_f1'], float)
    w = np.sin(th)[None, :, None] * np.ones((len(ph), 1, 1))  # solid-angle weight (k, j, i)
    W = w.sum(axis=(0, 1))
    sm = lambda q: (q * w).sum(axis=(0, 1)) / W  # noqa: E731
    rhom = sm(rho)
    mvr = sm(rho * vr) / rhom
    vrrms = np.sqrt(sm(rho * (vr - mvr[None, None, :]) ** 2) / rhom)
    vrms = np.sqrt(sm(rho * (vr ** 2 + vt ** 2 + vp ** 2)) / rhom)
    Em = sm(E); Fr = sm(F1)
    hent = (5. / 3.) * ei / rho          # ideal monatomic gas enthalpy per mass (approximation, no ionisation terms)
    Fenth = sm(rho * vr * hent) - sm(rho * vr) * sm(hent)   # enthalpy flux minus the net mass-flux part
    Fke = sm(0.5 * rho * (vr ** 2 + vt ** 2 + vp ** 2) * vr)
    # photosphere: Eddington grey T(tau=2/3) = T_eff, i.e. a T_rad^4 c/4 ... use sigma T_rad^4 R^2 = L_loc/(4 pi) per column,
    # L_loc/(4pi) = <F_r r^2> over 1.2-1.4 R_ph(IC) (flux outside the photosphere, ~r^-2)
    Trad = (np.maximum(E, 1e-30) / A) ** 0.25
    sel = (r > 1.2 * RPH0) & (r < 1.4 * RPH0)
    Lloc = (F1[:, :, sel] * r[sel] ** 2).mean(axis=2)             # (k, j)
    G = SIG * Trad ** 4 * r[None, None, :] ** 2 - Lloc[:, :, None]  # >0 inside, <0 outside
    i0 = np.searchsorted(r, 0.7 * RPH0)
    Rph = np.full(Lloc.shape, np.nan)
    for k in range(G.shape[0]):
        for j in range(G.shape[1]):
            g = G[k, j, i0:]
            neg = np.where(g < 0)[0]
            if len(neg) == 0 or neg[0] == 0:
                continue
            n = neg[0]
            a, b = g[n - 1], g[n]
            Rph[k, j] = r[i0 + n - 1] + (r[i0 + n] - r[i0 + n - 1]) * a / (a - b)
    Teff = (np.maximum(Lloc, 0) / Rph ** 2 / SIG) ** 0.25
    ns = np.searchsorted(r, RMAX_SL * RPH0)
    jm = len(th) // 2; km = len(ph) // 2
    f16 = lambda q: q[..., :ns].astype(np.float16)  # noqa: E731
    eq = lambda q: 0.5 * (q[:, jm - 1, :] + q[:, jm, :])  # noqa: E731   (phi, r) at theta = pi/2
    me = lambda q: 0.5 * (q[km - 1, :, :] + q[km, :, :])  # noqa: E731   (theta, r) at the mid phi
    out = dict(t=t, r=r, th=th, ph=ph, rf=rf, ns=ns,
               rho=rhom, mvr=mvr, vrrms=vrrms, vrms=vrms, E=Em, Trad=(Em / A) ** 0.25, Fr=Fr, Fenth=Fenth, Fke=Fke,
               Rph=Rph.astype(np.float32), Teff=Teff.astype(np.float32), Lloc=Lloc.astype(np.float32),
               eq_lrho=f16(np.log10(eq(rho))), eq_T=f16(eq(Trad) / 1e3), eq_vr=f16(eq(vr) / 1e5),
               me_lrho=f16(np.log10(me(rho))), me_T=f16(me(Trad) / 1e3), me_vr=f16(me(vr) / 1e5),
               f_p0_vr=f16(vr[0] / 1e5), f_p1_vr=f16(vr[-1] / 1e5),          # phi faces (theta, r)
               f_t0_vr=f16(vr[:, 0] / 1e5), f_t1_vr=f16(vr[:, -1] / 1e5),    # theta faces (phi, r)
               f_p0_T=f16(Trad[0] / 1e3), f_p1_T=f16(Trad[-1] / 1e3), f_t0_T=f16(Trad[:, 0] / 1e3), f_t1_T=f16(Trad[:, -1] / 1e3))  # kK (float16)
    np.savez(o + '.tmp.npz', **out)
    os.replace(o + '.tmp.npz', o)
    print('reduced', idx, 't %.3f d' % (t / 86400), 'R_ph med %.2f Rsun' % (np.nanmedian(Rph) / RS),
          'Teff med %.0f' % np.nanmedian(Teff), 'vr rms max %.3f km/s' % (vrrms.max() / 1e5), flush=True)
