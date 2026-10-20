#!/usr/bin/env python3
"""ACC3D page 10-20: orbital-plane movie frames (READ ONLY on the run). usage: mk_movie.py <framedir>
One frame per hydro_w dump that has an m1 dump at the same time (m1 every 0.01 code = every 2nd hydro dump), static + stream;
only new frames are drawn (FORCE=1 redraws). Columns: log rho, log T_gas (EOS table), log P_rad/P_gas (P_rad = E/3 from m1,
P_gas from the EOS table), v_phi (rotating frame). Row 1: orbital plane (mean of the two theta rows next to it), r <= d_L1;
row 2: edge-on meridional plane through the donor axis (phi 0 on +x, phi 180 on -x). Frame files f_<hydro dump no.>.png."""
import os
import sys
import numpy as np
from multiprocessing import Pool
import common as C

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
RMAX = C.DL1
SPEC = [('rho', r'log$_{10}\rho$ [code]', 'magma', -8, 2.5),
        ('lT', r'log$_{10}T_{\rm gas}$ [K]', 'inferno', 3.6, 6.0),
        ('lbeta', r'log$_{10}P_{\rm rad}/P_{\rm gas}$', 'RdBu_r', -2, 5),
        ('vphi', r'v$_\phi$ [km/s] (rotating frame)', 'coolwarm', -300, 300)]


def frame(arg):
    th_, f, mf, ph = arg
    fn = f'{OUT}/f_{C.dnum(f):05d}.png'
    if os.path.exists(fn) and not C.FORCE:
        return None
    plt = C.mpl()
    h = C.load(f, ['dens', 'eint', 'velz'])
    tq = C.thermo(h, mf)
    nj = h['dens'].shape[1]; js = [nj // 2 - 1, nj // 2]
    q = {'rho': np.log10(np.maximum(h['dens'], 1e-30)), 'lT': tq['lT'], 'lbeta': tq['lbeta'], 'vphi': h['velz']}
    rf, pf, tf = h['x1f'], h['x3f'], h['x2f']
    ir = int(np.searchsorted(rf, RMAX))
    PP, RR = np.meshgrid(pf, rf[:ir + 1], indexing='ij'); X, Y = RR * np.cos(PP), RR * np.sin(PP)
    TT, RR2 = np.meshgrid(tf, rf[:ir + 1], indexing='ij')
    k0 = 0; k1 = int(np.argmin(np.abs(h['x3v'] - np.pi)))
    fig, axs = plt.subplots(2, 4, figsize=(20, 8.6), gridspec_kw=dict(height_ratios=[1.0, 0.5]))
    for a, (key, lab, cm, lo, hi) in zip(axs[1], SPEC):
        for kk, sg in ((k0, 1), (k1, -1)):
            pc = a.pcolormesh(sg * RR2 * np.sin(TT), RR2 * np.cos(TT), q[key][kk, :, :ir], cmap=cm, **C.nrm(key, lo, hi),
                              shading='flat', rasterized=True)
        a.set_xlim(-RMAX, RMAX); a.set_ylim(-7, 7); a.set_aspect('equal')
        a.set_title(r'edge-on ($\phi$ 0 / 180): ' + lab, fontsize=9); a.set_xlabel(r'x [R$_\odot$]'); a.set_ylabel(r'z [R$_\odot$]')
        fig.colorbar(pc, ax=a, shrink=.85)
    for a, (key, lab, cm, lo, hi) in zip(axs[0], SPEC):
        pc = a.pcolormesh(X, Y, q[key][:, js, :ir].mean(axis=1), cmap=cm, **C.nrm(key, lo, hi), shading='flat', rasterized=True)
        a.add_patch(plt.Circle((0, 0), rf[0], color='#3a424c'))
        a.add_patch(plt.Circle((0, 0), RMAX, fill=False, ec=C.C2, lw=0.8, ls='--'))
        a.set_xlim(-RMAX, RMAX); a.set_ylim(-RMAX, RMAX); a.set_aspect('equal')
        a.set_title(lab, fontsize=10); a.set_xlabel(r'x [R$_\odot$]'); a.set_ylabel(r'y [R$_\odot$]')
        a.annotate(r'L1 radius (dashed); donor $\rightarrow$', xy=(RMAX * .97, -RMAX * .97), ha='right', fontsize=8, color=C.C2)
        fig.colorbar(pc, ax=a, shrink=.85)
    t = h['Time']
    phs = 'static phase' if t < C.TS - 1e-9 else 'stream on for %.3f P_orb' % ((t - C.TS) / C.PORB)
    fig.suptitle(f'RY Per accretor 3-D, {C.RUN_LABEL} ({C.SITE}): t = {t:.4f} code = {t / C.PORB:.3f} P_orb ({phs})', fontsize=12)
    fig.tight_layout(); fig.savefig(fn + '.tmp.png', dpi=72); plt.close(fig)
    os.replace(fn + '.tmp.png', fn)
    return fn


if __name__ == '__main__':
    ML = C.dumps('m1')
    L = [(t, f, C.m1_for(t, ML), ph) for t, f, ph in C.dumps('hydro_w')]
    L = [x for x in L if x[2]]
    C.tgas()                      # parse the EOS table once (cached) before forking
    n = 0
    with Pool(C.NP) as p:
        for fn in p.imap(frame, L):
            if fn:
                n += 1
    print('movie2d frames: %d paired dumps, %d new' % (len(L), n))
