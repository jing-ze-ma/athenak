"""Movie frames of the BSG page's Fig 3 (phi = 0 plane: ram/total pressure, kappa_R) and
Fig 4 (theta-phi shell at r ~ R_SLICE (bsg_figs.py, 36 Rsun): v_r, rho, T) for every dump of an
arm (read only).  Reuses bsg_figs.py's cached per-dump reductions.  The frames are the static
figures of bsg_figs.py with a time stamp: his layout and his fixed colour limits.  Only frames
that do not exist yet are rendered (BSG_RECHOOSE=1: all).  Also writes
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
    ap.add_argument('--rmax-plot', type=float, default=52.0)
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
    # colour limits are fixed constants (his, bsg_figs.py); BSG_RECHOOSE=1 redraws every frame
    redraw = B.RECHOOSE
    rng = B.int_ranges()
    stats, prev = [], None
    for n, h, m, v in dumps:
        f3 = os.path.join(MOV, 'frames', 'f3', 'f3_%05d.png' % n)
        f4 = os.path.join(MOV, 'frames', 'f4', 'f4_%05d.png' % n)
        red = B.reduce_dump(n, h, m, x1f, kap, v)
        if redraw or not os.path.exists(f3):
            B.fig3(red, B.F3_RMAX, kap, f3, stamp=True)
        if redraw or not os.path.exists(f4):
            B.fig4(red, f4, stamp=True)
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
