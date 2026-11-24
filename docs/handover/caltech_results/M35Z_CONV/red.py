#!/usr/bin/env python3
"""M35ZAMS convection diagnostic (orion NOTE 2026-11-24zzzo): reduce one hydro_w + vet dump (READ ONLY) to shell statistics.
usage: red.py OUTDIR N [N...]   (dump numbers, e.g. 43)   D = prod_rc4v_1124/m35zams (T2 run)
Units cgs.  T = (gamma-1) eint/rho * 7.4809404811e-9 K (deck temp_unit_kelvin); gas P = 2/3 eint; rad P = E/3 (P_rr = D11 E on
vet_vet dumps, else E/3).  Weights: solid angle (theta-phi); velocity statistics mass weighted by rho."""
import sys, os, glob
import numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
D = '/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/m35zams/bin/m35zams.'
TAB = '/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
C, A, TU, GM = 2.99792458e10, 7.5657332e-15, 7.4809404811e-9, 4.6449357021e27
L_ = open(TAB).readlines(); hd = [x for x in L_ if x.startswith('#')]
g = [x.split() for x in hd if x.startswith('# ') and len(x.split()) == 7 and x.split()[1].isdigit()][-1]
nT, nD, lT0, dlT, lD0, dlD = int(g[1]), int(g[2]), float(g[3]), float(g[4]), float(g[5]), float(g[6])
K = np.array([float(x) for x in L_ if not x.startswith('#')]).reshape(nT, nD)


def kap(T, r):
    x = np.clip((np.log10(np.maximum(T, 1)) - lT0) / dlT, 0, nT - 1.001); y = np.clip((np.log10(np.maximum(r, 1e-30)) - lD0) / dlD, 0, nD - 1.001)
    i, j = x.astype(int), y.astype(int); fx, fy = x - i, y - j
    return 10 ** (K[i, j] * (1 - fx) * (1 - fy) + K[i + 1, j] * fx * (1 - fy) + K[i, j + 1] * (1 - fx) * fy + K[i + 1, j + 1] * fx * fy)


def grid(h):
    """block order -> (nph, nth, nr) arrays via the coordinate arrays"""
    nb = h['n_mbs']
    x1f = np.array(h['mb_x1f'][0]); x1v = np.array(h['mb_x1v'][0])
    ths = sorted({round(float(v), 12) for b in range(nb) for v in h['mb_x2v'][b]})
    phs = sorted({round(float(v), 12) for b in range(nb) for v in h['mb_x3v'][b]})
    it = {v: i for i, v in enumerate(ths)}; ip = {v: i for i, v in enumerate(phs)}
    nr, nt, npp = len(x1v), len(ths), len(phs)
    dth = np.empty(nt); dph = np.empty(npp); sinth = np.empty(nt)
    for b in range(nb):
        for v, e0, e1 in zip(h['mb_x2v'][b], h['mb_x2f'][b][:-1], h['mb_x2f'][b][1:]):
            i = it[round(float(v), 12)]; dth[i] = np.cos(float(e0)) - np.cos(float(e1))   # solid-angle weight in theta
            sinth[i] = np.sin(float(v))
        for v, e0, e1 in zip(h['mb_x3v'][b], h['mb_x3f'][b][:-1], h['mb_x3f'][b][1:]):
            dph[ip[round(float(v), 12)]] = float(e1) - float(e0)
    def asm(var, name=h):
        out = np.empty((npp, nt, nr), np.float64)
        for b in range(nb):
            blk = np.asarray(name['mb_data'][var][b], dtype=np.float64)
            js = [it[round(float(v), 12)] for v in h['mb_x2v'][b]]; ks = [ip[round(float(v), 12)] for v in h['mb_x3v'][b]]
            out[np.ix_(ks, js)] = blk
        return out
    W = dth[None, :] * dph[:, None]
    return x1f, x1v, ths, phs, W, asm


def go(n, OUT):
    h = bc.read_binary(D + 'hydro_w.%05d.bin' % n)
    v = bc.read_binary(D + 'vet.%05d.bin' % n)
    x1f, x1v, ths, phs, W, asm = grid(h)
    rho, vr, vt, vp, ei = (asm(k, h) for k in ('dens', 'velx', 'vely', 'velz', 'eint'))
    E = asm('vet_e', v); F1 = asm('vet_f1', v)
    P_rr = None
    for fv in sorted(glob.glob(D + 'vet_vet.*.bin')):
        hv = bc.read_binary(fv)
        if abs(hv['time'] - h['time']) < 1.0:
            P_rr = asm('vet_d11', hv) * E
    hasD = P_rr is not None
    if P_rr is None:
        P_rr = E / 3.0
    nr = len(x1v); W3 = W[:, :, None]; WT = W.sum()
    sm = lambda q: (q * W3).sum((0, 1)) / WT                       # solid-angle shell mean
    Pg = (2. / 3.) * ei; T = Pg / rho * TU                         # = T_gas [K]
    kapR = kap(T, rho); chi = kapR * rho; dr = np.diff(x1f)
    tau = np.cumsum((chi * dr[None, None, :])[:, :, ::-1], axis=2)[:, :, ::-1]   # from the top face inward, at each cell's inner face
    Ptot = Pg + E / 3.0
    Csg2 = 5. / 3. * Pg / rho
    beta = Pg / Ptot; G1 = beta + (4 - 3 * beta) ** 2 * (2. / 3.) / (beta + 12 * (2. / 3.) * (1 - beta))
    Cst2 = G1 * Ptot / rho
    cm = lambda q, w: (q * w * rho * W3).sum((0, 1))               # rho*solid-angle weighted sums
    mom = lambda q, w: cm(q, w) / cm(np.ones_like(q), w)
    out = dict(n=n, t=float(h['time']), cyc=int(h['cycle']), r=x1v, rf=x1f, nth=len(ths), nph=len(phs), hasD=hasD,
               E_med_over_aT4=float(np.median(E[:, :, 50] / (A * T[:, :, 50] ** 4))))
    out['rho'] = sm(rho); out['T'] = sm(T); out['Pg'] = sm(Pg); out['Er'] = sm(E); out['F1'] = sm(F1); out['Prr'] = sm(P_rr)
    out['kapR'] = sm(chi) / sm(rho); out['tau_sh'] = np.cumsum((sm(chi) * dr)[::-1])[::-1]
    # per-column radii of tau = 1, 2/3 (linear in r inside the cell)
    for lev, nm in ((1.0, 'tau1'), (2. / 3., 'tau23')):
        k = nr - 1 - np.argmax(tau[:, :, ::-1] >= lev, axis=2)
        rows = np.ogrid[:rho.shape[0], :rho.shape[1]]
        tin = tau[rows[0], rows[1], k]; tout = tin - chi[rows[0], rows[1], k] * dr[k]
        rr = x1f[k + 1] - (lev - tout) / (tin - tout) * dr[k]
        out['rc_' + nm] = np.array([np.average(rr.ravel(), weights=W.ravel()), rr.min(), rr.max()])
    out['tau_col_med'] = np.array([np.median(tau[:, :, i]) for i in range(nr)])
    Vph = [vt, vp]; vars_all = {'vr': vr, 'vt': vt, 'vp': vp}
    for tag, mask in (('', np.ones_like(rho)), ('_m1', (tau >= 1.0).astype(float))):
        w = mask; wm = np.maximum(w, 0)
        den = cm(np.ones_like(rho), wm)
        ok = den > 0
        one = np.where(ok, 1, np.nan)
        mvr = np.where(ok, cm(vr, wm) / np.where(ok, den, 1), np.nan)
        out['mvr' + tag] = mvr
        dvr = vr - np.nan_to_num(mvr)[None, None, :]
        out['vr_rms_tot' + tag] = np.sqrt(mom(vr ** 2, wm)) * one
        out['vr_rms_fl' + tag] = np.sqrt(mom(dvr ** 2, wm)) * one
        out['vt_rms' + tag] = np.sqrt(mom(vt ** 2, wm)) * one
        out['vp_rms' + tag] = np.sqrt(mom(vp ** 2, wm)) * one
        mtv = np.where(ok, cm(vt, wm) / np.where(ok, den, 1), 0); mpv = np.where(ok, cm(vp, wm) / np.where(ok, den, 1), 0)
        out['vt_rms_fl' + tag] = np.sqrt(mom((vt - mtv[None, None, :]) ** 2, wm)) * one
        out['vp_rms_fl' + tag] = np.sqrt(mom((vp - mpv[None, None, :]) ** 2, wm)) * one
        sp2 = vr ** 2 + vt ** 2 + vp ** 2
        out['mach_g_rms' + tag] = np.sqrt(mom(sp2 / Csg2, wm)) * one
        out['mach_g_max' + tag] = np.array([np.max(np.where(mask[:, :, i] > 0, np.sqrt(sp2 / Csg2)[:, :, i], 0)) for i in range(nr)])
        out['mach_t_rms' + tag] = np.sqrt(mom(sp2 / Cst2, wm)) * one
        out['mach_t_max' + tag] = np.array([np.max(np.where(mask[:, :, i] > 0, np.sqrt(sp2 / Cst2)[:, :, i], 0)) for i in range(nr)])
        Tm = mom(T, wm)
        out['dT_rms' + tag] = np.sqrt(mom((T / Tm[None, None, :] - 1) ** 2, wm)) * one
    # KE per unit radius pieces (radial mean flow / fluctuation / lateral): sums over the wedge, dV = r^2 dr dOmega
    dV = (x1f[1:] ** 3 - x1f[:-1] ** 3) / 3.0
    mv = sm(rho * vr) / sm(rho)
    KEr_mean = 0.5 * sm(rho) * mv ** 2 * 4 * np.pi * dV            # scaled to 4 pi from the wedge mean
    KEr_fl = 0.5 * 4 * np.pi * dV * sm(rho * (vr - mv[None, None, :]) ** 2)
    KElat = 0.5 * 4 * np.pi * dV * sm(rho * (vt ** 2 + vp ** 2))
    mt_ = sm(rho * vt) / sm(rho); mp_ = sm(rho * vp) / sm(rho)     # shell-mean lateral flow (shear / ray-effect dipole) split off
    out['mvt'] = mt_; out['mvp'] = mp_
    out['KElat_fl'] = 0.5 * 4 * np.pi * dV * sm(rho * ((vt - mt_[None, None, :]) ** 2 + (vp - mp_[None, None, :]) ** 2))
    out.update(KEr_mean=KEr_mean, KEr_fl=KEr_fl, KElat=KElat, dV=dV)
    # fluxes, scaled to the full sphere [erg/s]: L = 4 pi r^2 <.>_Omega
    r2 = 4 * np.pi * x1v ** 2
    h_tot = (ei + Pg + E + P_rr) / rho                              # total (gas + radiation) specific enthalpy
    out['Lrad_lab'] = r2 * sm(F1)
    Lcom = r2 * sm(F1 - (E + P_rr) * vr)                            # comoving radiative flux (cell F1, E, P_rr; estimate)
    out['Lrad_com'] = Lcom
    rvr = rho * vr
    out['Lenth'] = r2 * (sm(rvr * h_tot) - sm(rvr) * sm(h_tot))      # resolved enthalpy flux, mean mass flux removed
    ke = 0.5 * (vr ** 2 + vt ** 2 + vp ** 2)
    out['Lkin'] = r2 * (sm(rvr * ke) - sm(rvr) * sm(ke))
    out['mdot'] = r2 * sm(rvr)
    # store compact slice data
    iFe = None
    for nm, arr in (('Tprof', out['T']), ('rhop', out['rho'])):
        out[nm] = arr
    np.savez_compressed(os.path.join(OUT, 'm35_%05d.npz' % n), **{k_: v_ for k_, v_ in out.items()})
    print('done', n, 't %.1f' % out['t'], 'cyc', out['cyc'], 'E/aT4 med(deep)', round(out['E_med_over_aT4'], 4),
          'tau_col top/mid', out['tau_col_med'][-1], flush=True)


if __name__ == '__main__':
    os.makedirs(sys.argv[1], exist_ok=True)
    for n_ in sys.argv[2:]:
        go(int(n_), sys.argv[1])
