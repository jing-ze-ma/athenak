"""cmp57_dumps.py N [N ...]: per-dump reductions for the 57.4 d comparison with Ma+2026
(10-08), read only from the merged run view $BSG_REPRO_RUN; cached in
scripts_repro/cache/cmp57/c_NNNNN.npz.  Per dump:
  - tau(r) per column from the top (A9, cell centres, Ma's kappa_R); its 2.5/50/97.5
    percentiles over columns per radius (paper Fig. 5); own photosphere per column = first
    cell outward with tau < 1 (as the intensity map); r_own, |v|, v_r there (paper Fig. 7).
  - whole-star-equivalent luminosities per radius (4 pi / Omega * sum r^2 dOmega F):
    L_rad (lab-frame m1_f1), L_enth_gas = 5/2 P_g v_r, L_enth_rad = 4/3 E v_r, L_kin.
  - shell moments: <rho>, <rho^2>/<rho>^2 (clumping), rms(rho)/<rho>, vrms, vr rms.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsg_figs as B  # noqa: E402

RUN = os.environ.get('BSG_REPRO_RUN', '/viper/ptmp2/jinma/bsg_viper_1006/merged_repro')
OUT = os.path.join(B.FIGS, 'cache', 'cmp57')
VB = np.arange(0, 151, 1.0)     # |v| histogram bins [km/s]


def one(n, x1f, kap):
    cf = os.path.join(OUT, 'c_%05d.npz' % n)
    if os.path.exists(cf):
        return
    hf = os.path.join(RUN, 'bin/bsg3d.hydro_w.%05d.bin' % n)
    mf = os.path.join(RUN, 'bin/bsg3d.m1.%05d.bin' % n)
    dh = B.bc.read_binary(hf)
    t = dh['time']
    w = B.assemble(dh)
    del dh
    nx3, nx2, nx1 = w['dens'].shape
    x2f = np.linspace(np.pi / 3, 2 * np.pi / 3, nx2 + 1)
    x3f = np.linspace(0.0, np.pi / 3, nx3 + 1)
    dom = np.diff(x3f)[:, None] * (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :]
    omega = dom.sum()
    r = 0.5 * (x1f[1:] + x1f[:-1])
    dr = np.diff(x1f)
    rho = w['dens']
    pg = (B.GAMMA - 1.0) * w['eint']
    T = pg / rho * np.float32(B.MU_MU_OVER_K)
    vr = w['velx']
    v2 = vr ** 2 + w['vely'] ** 2 + w['velz'] ** 2
    del w
    kr = B.Kappa.__call__(kap, rho, T)
    tau = np.flip(np.cumsum(np.flip((kr * rho).astype(np.float64) * dr[None, None, :],
                                    axis=2), axis=2), axis=2)
    del kr, T
    res = dict(n=n, time=t, r=r, omega=omega)
    tf = tau.reshape(-1, nx1)
    res['tau_p'] = np.percentile(tf, [2.5, 50, 97.5], axis=0)
    res['tau_mean'] = np.einsum('kji,kj->i', tau, dom) / omega
    thin = tau < 1.0
    iown = np.where(thin.any(axis=2), np.argmax(thin, axis=2), nx1 - 1)
    del tau, thin, tf
    res['r_own'] = r[iown].astype(np.float32)
    vown = np.sqrt(np.take_along_axis(v2, iown[..., None], 2)[..., 0]) / 1e5
    res['v_own'] = vown.astype(np.float32)
    res['vr_own'] = (np.take_along_axis(vr, iown[..., None], 2)[..., 0] / 1e5).astype(
        np.float32)
    res['v_own_hist'] = np.histogram(vown, VB, weights=dom)[0] / omega
    sh = lambda a: np.einsum('kji,kj->i', a.astype(np.float64), dom) / omega  # noqa: E731
    res['rho_sh'] = sh(rho)
    res['rho2_sh'] = sh(rho * rho)
    res['vrms_sh'] = np.sqrt(sh(v2))
    res['vr_rms_sh'] = np.sqrt(sh(vr * vr))
    fac = 4 * np.pi * r ** 2       # whole-star equivalent: 4 pi r^2 <F>_Omega
    res['L_enth_gas'] = fac * sh(np.float32(B.GAMMA / (B.GAMMA - 1.0)) * pg * vr)
    res['L_kin'] = fac * sh(0.5 * rho * v2 * vr)
    del pg, v2
    dm = B.bc.read_binary(mf)
    assert abs(dm['time'] - t) < 1e-6 * max(1.0, t)
    a = B.assemble(dm)
    del dm
    res['L_rad'] = fac * sh(a['m1_f1'])
    res['L_enth_rad'] = fac * sh(np.float32(4.0 / 3.0) * a['m1_e'] * vr)
    res['mdot_whole'] = fac * sh(rho * vr)      # g/s, whole-star equivalent
    del a, rho, vr
    np.savez(cf, **res)
    print('cmp57 %05d t = %.2f d' % (n, t / B.DAY), flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    B.CACHE = os.path.join(B.FIGS, 'cache', 'cmp57')
    x1f, _ = B.radial_faces(RUN)
    kap = B.Kappa()
    for a in sys.argv[1:]:
        one(int(a), x1f, kap)


if __name__ == '__main__':
    main()
