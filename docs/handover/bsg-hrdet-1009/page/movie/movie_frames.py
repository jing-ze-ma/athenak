"""Movie frames of the BSG page's Fig 3 (phi = 0 plane: ram/total pressure, kappa_R) and
Fig 4 (theta-phi shell at r ~ R_SLICE (bsg_figs.py, 40.6 Rsun): v_r, rho and T deviations) for every dump of an
arm (default arm2s2 = R2+B, arm2r2b_4n; read only).  Reuses bsg_figs.py's cached
per-dump reductions (figs/cache/<arm>/red_NNNNN.npz).  Colour ranges are FIXED for all
frames: chosen once from the newest dump at the first run and stored in ranges.json
(delete it to re-choose).  Only frames that do not exist yet are rendered.  Also writes
cadence.json: v_r rms, plume size and turnover time at the shell, and the frame-to-frame
correlation of v_r.

usage: python3 movie_frames.py [--arm arm2s2]
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import bsg_figs as B  # noqa: E402
plt = B.plt
matplotlib = B.matplotlib
MOV = os.path.dirname(os.path.abspath(__file__))


def shell_fields(red):
    """relative fluctuations [%] (used to choose the ranges and for cadence stats)"""
    vr = red['s36_vr'] / 1e5
    rho = (red['s36_rho'] / red['s36_rho'].mean() - 1) * 100
    T = (red['s36_T'] / red['s36_T'].mean() - 1) * 100
    return vr, rho, T


def choose_ranges(red):
    vr, rho, T = shell_fields(red)
    out = {}
    for nm, z in (('vr', vr), ('rho', rho), ('T', T)):
        lo, hi = np.percentile(z, [1, 99])
        m = max(float(max(abs(lo), abs(hi))), {"vr": 0.1}.get(nm, B.REL_MIN_PCT))
        out[nm] = [-m, m]
    out['ram'] = [0.0, 0.5]
    out['kap'] = [0.05, 2.5]
    out['from_dump'] = int(red['n'])
    out['t_d'] = float(red['time'] / B.DAY)
    return out


def frame3(red, rng, rmax_plot, out):
    rf = red['x1f'] / B.RSUN
    x2f = red['x2f']
    ni = int((rf <= rmax_plot + 1e-9).sum()) - 1
    R, TH = np.meshgrid(rf[:ni + 1], x2f)
    X, Y = R * np.cos(TH), R * np.sin(TH)
    pram = 0.5 * red['p0_rho'] * red['p0_v2']
    prr, psrc = B.prad_rr(red)
    ratio = (pram / (red['p0_pg'] + prr))[:, :ni]
    kp = red['p0_kap'][:, :ni]
    rph = B.r_photo(red['r'], red['tau_mean']) / B.RSUN
    th = np.linspace(np.pi / 3, 2 * np.pi / 3, 200)
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 6.0))
    blues = B.LinearSegmentedColormap.from_list('ramp', ['white', 'blue', 'navy'])
    plab = {'m1_vet': r'$P_{\rm rad}=D_{rr}E$', 'E/3': r'$P_{\rm rad}=E/3$'}.get(psrc, psrc)
    for ax, z, cm, (lo, hi), ext, tl, cl in (
            (axs[0], ratio, blues, rng['ram'], 'max', 'Ram Pressure vs Total Pressure',
             r'$\frac{1}{2}\rho v^2/(P_{\rm rad}+P_{\rm gas})$,  ' + plab),
            (axs[1], kp, 'GnBu', rng['kap'], 'both', 'Rosseland Opacity',
             r'$\kappa_{\rm Rosseland}$ [cm$^2$g$^{-1}$]  (Ma+2026 TOPS tables)')):
        pc = ax.pcolormesh(X, Y, z, cmap=cm, vmin=lo, vmax=hi, shading='flat',
                           rasterized=True)
        B.wedge_axes(ax, rmax_plot)
        if rph == rph:
            ax.plot(rph * np.cos(th), rph * np.sin(th), 'b--', lw=2)
        ax.set_title(tl, fontsize=12)
        cb = fig.colorbar(pc, ax=ax, orientation='horizontal', extend=ext, shrink=0.8,
                          pad=0.02)
        cb.set_label(cl)
    fig.suptitle(r'$\phi$ = 0 plane    t = %6.2f d' % (red['time'] / B.DAY), fontsize=20,
                 family='monospace', y=0.97)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.04, wspace=0.05)
    fig.savefig(out, dpi=100)
    plt.close(fig)


def frame4(red, rng, out):
    x2f, x3f = np.degrees(red['x2f']), np.degrees(red['x3f'])
    r36 = red['r'][int(red['i36'])] / B.RSUN
    vr = red['s36_vr'] / 1e5
    rho = red['s36_rho'] / 1e-8
    T = red['s36_T'] / 1e5
    # rho, T absolute; fixed range = shell mean of the newest dump x (1 -+ a)
    lims = dict(vr=rng['vr'], rho=[x / 1e-8 for x in rng['rho_abs']],
                T=[x / 1e5 for x in rng['T_abs']])
    fig, axs = plt.subplots(1, 3, figsize=(14.4, 6.4), constrained_layout=True)
    for ax, z, cm, key, tl, cl in (
            (axs[0], vr, 'seismic', 'vr', 'Radial Velocity', r'$v_r$ [km/s]'),
            (axs[1], rho, 'Blues', 'rho', 'Gas Density',
             r'$\rho$ [$10^{-8}$ g cm$^{-3}$]  ($\bar\rho\,(1\pm%.3g\%%)$, '
             r'$\bar\rho$ = %.3f)' % (rng['rho'][1], rng['rho_ref'] / 1e-8)),
            (axs[2], T, 'Blues_r', 'T', 'Gas Temperature',
             r'$T$ [$10^5$ K]  ($\bar T\,(1\pm%.3g\%%)$, $\bar T$ = %.4f)'
             % (rng['T'][1], rng['T_ref'] / 1e5))):
        lo, hi = lims[key]
        # z is (k = phi, j = theta): transpose for (x = phi, y = theta)
        pc = ax.pcolormesh(x3f, x2f, z.T, cmap=cm, vmin=lo, vmax=hi, shading='flat',
                           rasterized=True)
        ax.set_xlim(0, 60)
        ax.set_ylim(120, 60)
        ax.set_xticks([0, 15, 30, 45, 60])
        ax.set_yticks([60, 75, 90, 105, 120])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\phi$ [deg]')
        ax.set_title(tl, fontsize=12)
        cb = fig.colorbar(pc, ax=ax, orientation='horizontal', extend='both', shrink=0.9)
        cb.set_label(cl)
        cb.formatter.set_useOffset(False)
        cb.locator = matplotlib.ticker.MaxNLocator(5)
        cb.update_ticks()
    axs[0].set_ylabel(r'$\theta$ [deg]')
    fig.suptitle(r'shell r = %.2f $R_\odot$    t = %6.2f d' % (r36, red['time'] / B.DAY),
                 fontsize=20, family='monospace')
    fig.savefig(out, dpi=100)
    plt.close(fig)


def plume_size(vr, dth, dph):
    """first zero crossing of the azimuthally averaged 2-D autocorrelation of v_r [rad]."""
    z = vr - vr.mean()
    ny, nx = z.shape
    f = np.fft.rfft2(z, s=(2 * ny, 2 * nx))
    ac = np.fft.irfft2(np.abs(f) ** 2, s=(2 * ny, 2 * nx))
    cnt = np.fft.irfft2(np.abs(np.fft.rfft2(np.ones_like(z), s=(2 * ny, 2 * nx))) ** 2,
                        s=(2 * ny, 2 * nx))
    ac = ac / np.maximum(cnt, 1)
    ac = np.fft.fftshift(ac / ac[0, 0])
    yy, xx = np.indices(ac.shape)
    d = np.hypot((yy - ny) * dth, (xx - nx) * dph)
    bins = np.arange(0, 0.5, dth)
    k = np.digitize(d.ravel(), bins)
    prof = np.bincount(k, ac.ravel(), len(bins) + 1) / np.maximum(
        np.bincount(k, minlength=len(bins) + 1), 1)
    i0 = np.where(prof[1:] < 0)[0]
    return float(bins[i0[0]]) if len(i0) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--arm', default='repro')
    ap.add_argument('--rmax-plot', type=float, default=65.0)
    a = ap.parse_args()
    B.ARM = a.arm
    run = B.ARMS[a.arm]['run']
    B.IC, B.INP = B.ARMS[a.arm]['ic'], B.ARMS[a.arm]['inp']
    B.CACHE = B.cache_dir(run)
    os.makedirs(B.CACHE, exist_ok=True)
    x1f, _ = B.radial_faces(run)
    kap = B.Kappa()
    dumps = B.dump_list(run)
    for d in ('f3', 'f4', 'fint'):
        os.makedirs(os.path.join(MOV, 'frames', d), exist_ok=True)
    rf = os.path.join(MOV, 'ranges.json')
    # repro: bsg_figs.py (run first, BSG_RECHOOSE=1) writes the rho/T/intensity ranges;
    # the keys only this script uses (vr, ram, kap) are added from the newest dump.  With
    # BSG_RECHOOSE=1 every frame is redrawn (overwritten) with the new ranges.
    redraw = B.RECHOOSE
    B.RECHOOSE = False
    rng = json.load(open(rf)) if os.path.exists(rf) else {}
    if redraw or 'vr' not in rng:
        n, h, m, v = dumps[-1]
        ch = choose_ranges(B.reduce_dump(n, h, m, x1f, kap, v))
        for k in ('vr', 'ram', 'kap', 'from_dump', 't_d'):
            rng[k] = ch[k]
        for k in ('rho', 'T'):
            rng.setdefault(k, ch[k])
        json.dump(rng, open(rf, 'w'), indent=1)
    if 'rho_abs' not in rng or 'T_abs' not in rng:
        n, h, m, v = dumps[-1]
        rng = B.abs_ranges(B.reduce_dump(n, h, m, x1f, kap, v))
    nl = dumps[-1]
    irng = B.int_ranges(B.int_reduce(nl[0], nl[1], nl[2], x1f, kap,
                                     B.reduce_dump(nl[0], nl[1], nl[2], x1f, kap, nl[3])))
    rng.update({k: irng[k] for k in irng if k.startswith(('int_', 'teff_'))})
    stats, prev = [], None
    for n, h, m, v in dumps:
        f3 = os.path.join(MOV, 'frames', 'f3', 'f3_%05d.png' % n)
        f4 = os.path.join(MOV, 'frames', 'f4', 'f4_%05d.png' % n)
        red = B.reduce_dump(n, h, m, x1f, kap, v)
        if redraw or not os.path.exists(f3):
            frame3(red, rng, a.rmax_plot, f3)
        if redraw or not os.path.exists(f4):
            frame4(red, rng, f4)
        fi = os.path.join(MOV, 'frames', 'fint', 'fint_%05d.png' % n)
        if redraw or not os.path.exists(fi):
            B.fig_int(B.int_reduce(n, h, m, x1f, kap, red), rng, fi)
        vr = red['s36_vr'].astype(np.float64)
        r36 = float(red['r'][int(red['i36'])])
        dth = float(np.diff(red['x2f'])[0])
        dph = float(np.diff(red['x3f'])[0])
        lp = plume_size(vr, dth, dph)
        vrms = float(np.sqrt(np.mean((vr - vr.mean()) ** 2)))
        s = dict(n=int(n), t_d=float(red['time'] / B.DAY), vr_rms_kms=vrms / 1e5,
                 plume_corr_zero_rad=lp, plume_size_Rsun=lp * r36 / B.RSUN,
                 turnover_d=(lp * r36 / vrms / B.DAY) if vrms > 0 else np.nan)
        if prev is not None:
            s['corr_vr_prev'] = float(np.corrcoef(vr.ravel(), prev.ravel())[0, 1])
        prev = vr
        stats.append(s)
    json.dump(dict(ranges=rng, frames=stats), open(os.path.join(MOV, 'cadence.json'),
                                                    'w'), indent=1)
    for s in stats[-6:]:
        print(json.dumps(s))


if __name__ == '__main__':
    main()
