#!/usr/bin/env python3
"""tests_r21: is the horizontal isothermality that drives the porosity PHYSICAL?

Per-shell channel/clump sizes, transverse optical depths, timescales, EOS
consistency of the density contrast, and the effective-opacity reduction.
READ-ONLY on the run.  usage: chan.py [dumpindex ...]
"""
import os
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r19')
import bin_convert as bc                                          # noqa: E402
from poros import (RUN, RS, LSTAR, TURN, ARAD, CLIGHT, kappa,     # noqa: E402
                   grid, temp_of, TB, bil, cellidx, weight)

EOSTBL = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
          'eos_table_box_w5.txt')
TARGETS = [0.60, 0.72, 0.80, 0.88, 0.93, 0.96]

# ------------------------------------------------------------- pressure table
_d = np.loadtxt(EOSTBL, comments='#')
PTB = _d[:, 1].reshape(TB['ny'], TB['nx'])
del _d


def egas(rho, T):
    ix, iy, u, v = cellidx(np.log10(rho), np.log10(T))
    return rho*10.0**bil(TB['E'], ix, iy, u, v)


def pgas(rho, T):
    ix, iy, u, v = cellidx(np.log10(rho), np.log10(T))
    return rho*10.0**bil(PTB, ix, iy, u, v)


def ptot(rho, T):
    return pgas(rho, T) + weight(np.log10(rho))*ARAD*T**4/3.0


def utot(rho, T):
    return egas(rho, T) + weight(np.log10(rho))*ARAD*T**4


def thermo(rho, T):
    """p, u, beta, Gamma1, c_s, Cv(volumetric)."""
    p = ptot(rho, T)
    u = utot(rho, T)
    beta = pgas(rho, T)/p
    dl = 0.004
    chiT = (np.log(ptot(rho, T*10**dl)) - np.log(ptot(rho, T/10**dl)))/(2*dl*np.log(10))
    chir = (np.log(ptot(rho*10**dl, T)) - np.log(ptot(rho/10**dl, T)))/(2*dl*np.log(10))
    cv = (utot(rho, T*10**dl) - utot(rho, T/10**dl))/(T*(10**dl - 10**-dl))
    g1 = chir + chiT**2*p/(T*cv)
    return p, u, beta, g1, np.sqrt(g1*p/rho), cv


def dlnk_dlnrho(T, rho):
    d = 0.02
    return (np.log(kappa(T, rho*10**d)) - np.log(kappa(T, rho/10**d)))/(2*d*np.log(10))


def dlnk_dlnT(T, rho):
    d = 0.01
    return (np.log(kappa(T*10**d, rho)) - np.log(kappa(T/10**d, rho)))/(2*d*np.log(10))


# ------------------------------------------------------------- loading
def load(idx, fields=('dens', 'eint', 'velx', 'vely', 'velz')):
    d = bc.read_binary(os.path.join(RUN, 'bin', 'he4.hydro_w.%05d.bin' % idx))
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    out = {}
    nj, nk = d['nx2_mb'], d['nx3_mb']
    for nm in fields:
        g = np.empty((n3, n2, n1), dtype=np.float64)
        src = d['mb_data'][nm]
        for m in range(d['n_mbs']):
            _, lx2, lx3, _ = d['mb_logical'][m]
            g[lx3*nk:(lx3+1)*nk, lx2*nj:(lx2+1)*nj, :] = src[m]
        out[nm] = g
    del d['mb_data']
    return d, out


def amean(X, wj):
    """area-weighted mean of a shell map X (n3,n2)."""
    return float(np.dot(wj, X.mean(axis=0))/wj.sum())


def acf_len(X, axis):
    """periodic autocorrelation length (lags to 1/e) of shell map X along axis,
    averaged over the other axis.  Returns lag in cells (float)."""
    Y = X - X.mean(axis=axis, keepdims=True)
    F = np.fft.rfft(Y, axis=axis)
    C = np.fft.irfft(F*np.conj(F), n=Y.shape[axis], axis=axis)
    C = C.mean(axis=1-axis)
    C = C/C[0]
    n = len(C)//2
    for k in range(1, n):
        if C[k] <= np.exp(-1.0):
            return k - 1 + (C[k-1] - np.exp(-1.0))/(C[k-1] - C[k])
    return float(n)


# ------------------------------------------------------------- per-dump work
def run(idx, fh):
    d, g = load(idx)
    rc, dr, thc, wj = grid(d)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    r0 = rc/RS
    dth = (d['x2max'] - d['x2min'])/n2
    dph = (d['x3max'] - d['x3min'])/n3
    sinm = float(np.dot(wj, np.sin(thc))/wj.sum())
    ks = [int(np.argmin(np.abs(r0 - t))) for t in TARGETS]
    band = sorted(set(sum([[k-1, k, k+1] for k in ks], [])))
    sl = {k: i for i, k in enumerate(band)}

    rho = np.ascontiguousarray(g['dens'][:, :, band])
    ein = np.ascontiguousarray(g['eint'][:, :, band])
    vr = np.ascontiguousarray(g['velx'][:, :, band])
    vt = np.ascontiguousarray(g['vely'][:, :, band])
    vp = np.ascontiguousarray(g['velz'][:, :, band])
    del g
    T = temp_of(rho, ein)
    kap = kappa(T, rho)
    p, u, beta, g1, cs, cv = thermo(rho, T)

    t = d['time']
    print('\n' + '=' * 108, file=fh)
    print('dump %d   t = %.1f s = %.4f turnover   (grid %dx%dx%d, dtheta=%.5f '
          'dphi=%.5f rad)' % (idx, t, t/TURN, n1, n2, n3, dth, dph), file=fh)

    hdr1 = ('%6s %6s | %9s %9s %9s %9s | %8s %8s %8s | %8s %8s'
            % ('r/R', 'i', '<rho>', '<T>', 'beta', 'Gam1', 'Lth/cm', 'Lph/cm',
               'L/cells', 'rmsrho', 'rmsT'))
    rows = []
    for k in ks:
        i = sl[k]
        rh, Tk = rho[:, :, i], T[:, :, i]
        rm = amean(rh, wj)
        Tm = amean(Tk, wj)
        x = rh/rm
        lth_c = acf_len(x, 1)
        lph_c = acf_len(x, 0)
        dl_th = rc[k]*dth
        dl_ph = rc[k]*sinm*dph
        Lth, Lph = lth_c*dl_th, lph_c*dl_ph
        L = np.sqrt(Lth*Lph)
        Lcell = np.sqrt(lth_c*lph_c)
        # channel / clump populations
        q = np.percentile(rh, [10, 50, 90])
        lo = rh <= q[0]
        hi = rh >= q[2]
        rlo, rhi = rh[lo].mean(), rh[hi].mean()
        Tlo, Thi = Tk[lo].mean(), Tk[hi].mean()
        klo, khi = kappa(Tlo, rlo), kappa(Thi, rhi)
        tau_lo = klo*rlo*L
        tau_hi = khi*rhi*L
        tau_m = amean(kap[:, :, i], wj)*rm*L
        rrms = np.sqrt(amean((x - 1.0)**2, wj))
        trms = np.sqrt(amean((Tk/Tm - 1.0)**2, wj))
        # extreme channels
        q1 = np.percentile(rh, 1.0)
        m1 = rh.min()
        r1 = rh[rh <= q1].mean()
        T1 = Tk[rh <= q1].mean()
        k1 = kappa(T1, r1)
        kmn = kappa(Tk.ravel()[rh.argmin()], m1)
        krc = kap[:, :, i]*rh
        f_thin_cell = float((krc*dl_th < 1.0).mean())
        f_thin_L = float((krc*L < 1.0).mean())
        rows.append(dict(k=k, r=r0[k], rm=rm, Tm=Tm, L=L, Lth=Lth, Lph=Lph,
                         r1=r1, k1=k1, tau1=k1*r1*L, tau1c=k1*r1*dl_th,
                         rmin=m1, taumin=kmn*m1*L, mfp1=1.0/(k1*r1),
                         f_thin_cell=f_thin_cell, f_thin_L=f_thin_L,
                         Lcell=Lcell, lth_c=lth_c, lph_c=lph_c,
                         rlo=rlo, rhi=rhi, Tlo=Tlo, Thi=Thi, klo=klo, khi=khi,
                         tau_lo=tau_lo, tau_hi=tau_hi, tau_m=tau_m,
                         rrms=rrms, trms=trms, dl_th=dl_th, dl_ph=dl_ph))
    print('\n--- 1. channel/clump size and TRANSVERSE optical depth ---', file=fh)
    print(hdr1, file=fh)
    for R in rows:
        print('%6.3f %6d | %9.3e %9.4e %9.3e %9.4f | %8.3e %8.3e %8.2f | %8.4f %8.5f'
              % (R['r'], R['k'], R['rm'], R['Tm'],
                 amean(beta[:, :, sl[R['k']]], wj), amean(g1[:, :, sl[R['k']]], wj),
                 R['Lth'], R['Lph'], R['Lcell'], R['rrms'], R['trms']), file=fh)
    print('\n%6s | %10s %10s %9s | %10s %10s %9s | %9s'
          % ('r/R', 'rho_ch', 'kap_ch', 'tau_perp', 'rho_cl', 'kap_cl', 'tau_perp',
             'tau_mean'), file=fh)
    for R in rows:
        print('%6.3f | %10.3e %10.4f %9.3f | %10.3e %10.4f %9.3f | %9.3f'
              % (R['r'], R['rlo'], R['klo'], R['tau_lo'],
                 R['rhi'], R['khi'], R['tau_hi'], R['tau_m']), file=fh)

    print('\n%6s | %10s %10s %9s %9s | %10s %9s | %9s %9s'
          % ('r/R', 'rho_1%', 'rho_1%/<r>', 'tau_L', 'tau_cell', 'rho_min',
             'tau_L_min', 'f(dtau_c<1)', 'f(tauL<1)'), file=fh)
    for R in rows:
        print('%6.3f | %10.3e %10.4f %9.3f %9.3f | %10.3e %9.3f | %9.4f %9.4f'
              % (R['r'], R['r1'], R['r1']/R['rm'], R['tau1'], R['tau1c'],
                 R['rmin'], R['taumin'], R['f_thin_cell'], R['f_thin_L']), file=fh)
    print('  (mfp of the 1%% channel 1/(kappa rho) in cm, vs L: ' +
          '  '.join('%.3f:%.2e/%.2e' % (R['r'], R['mfp1'], R['L']) for R in rows),
          file=fh)

    # --- 2. timescales -------------------------------------------------------
    print('\n--- 2. timescales on the channel scale L (s) ---', file=fh)
    print('%6s | %10s %10s %8s | %10s %10s %10s %10s | %8s %8s'
          % ('r/R', 't_diff_ch', 't_diff_mn', 'Cv/4aT3', 't_sound', 't_conv_h',
             't_conv_tot', 't_coolrad', 'Mach_h', 'Mach_tot'), file=fh)
    trow = []
    for R in rows:
        i = sl[R['k']]
        L = R['L']
        cvm = amean(cv[:, :, i], wj)
        fcorr = cvm/(4.0*ARAD*R['Tm']**3)
        tdif_ch = 3.0*R['klo']*R['rlo']*L*L/CLIGHT*fcorr
        tdif_mn = 3.0*amean(kap[:, :, i]*rho[:, :, i], wj)*L*L/CLIGHT*fcorr
        csm = amean(cs[:, :, i], wj)
        vh = np.sqrt(amean(vt[:, :, i]**2 + vp[:, :, i]**2, wj))
        vt3 = np.sqrt(amean(vr[:, :, i]**2 + vt[:, :, i]**2 + vp[:, :, i]**2, wj))
        # radial radiative cooling time: U_tot / (F / H_T)
        Tprof = np.array([amean(T[:, :, sl[R['k']+dd]], wj) for dd in (-1, 0, 1)])
        dTdr = (Tprof[2] - Tprof[0])/(rc[R['k']+1] - rc[R['k']-1])
        HT = R['Tm']/abs(dTdr)
        F = LSTAR/(4*np.pi*rc[R['k']]**2)
        tcool = amean(u[:, :, i], wj)*HT/F
        tsnd, tcnv = L/csm, L/vh
        trow.append((R, tdif_ch, tdif_mn, tsnd, tcnv, tcool, vh/csm, vt3/csm,
                     vh, csm, fcorr, HT, dTdr, cvm))
        print('%6.3f | %10.3e %10.3e %8.3f | %10.3e %10.3e %10.3e %10.3e | '
              '%8.4f %8.4f'
              % (R['r'], tdif_ch, tdif_mn, fcorr, tsnd, tcnv, L/max(vt3, 1e-30),
                 tcool, vh/csm, vt3/csm), file=fh)
    print('\n%6s | %10s %10s %10s | ratios' % ('r/R', 'v_h(cm/s)', 'c_s', 'H_T/cm'),
          file=fh)
    print('%6s | %10s %10s %10s | %10s %10s %10s'
          % ('', '', '', '', 'tdif/tsnd', 'tdif/tcnv', 'tdif/tcool'), file=fh)
    for (R, td, tdm, ts, tc, tcl, mh, mt, vh, csm, fc, HT, dTdr, cvm) in trow:
        print('%6.3f | %10.3e %10.3e %10.3e | %10.4f %10.4f %10.4f'
              % (R['r'], vh, csm, HT, td/ts, td/tc, td/tcl), file=fh)

    # --- 3. EOS consistency --------------------------------------------------
    print('\n--- 3. EOS decomposition of the density contrast ---', file=fh)
    print('  dlnrho = (1/beta) dlnp - (1 + 4(1-beta)/beta) dlnT  (exact for '
          'p=pg+aT^4/3, pg ~ rho T)', file=fh)
    print('%6s | %9s %9s %9s %9s | %9s %9s %9s | %8s %8s %8s'
          % ('r/R', 'rms dlnr', 'rms dlnp', 'rms dlnT', 'beta', 'p-term',
             'T-term', 'pred rms', 'M_h', 'M^2/beta', 'rhov2/p'), file=fh)
    for R in rows:
        i = sl[R['k']]
        rh, Tk, pk = rho[:, :, i], T[:, :, i], p[:, :, i]
        rm, Tm = R['rm'], R['Tm']
        pm = amean(pk, wj)
        bm = amean(beta[:, :, i], wj)
        dlr = np.log(rh/rm)
        dlp = np.log(pk/pm)
        dlt = np.log(Tk/Tm)
        cT = 1.0 + 4.0*(1.0 - bm)/bm
        pterm = dlp/bm
        tterm = -cT*dlt
        pred = pterm + tterm
        csm = amean(cs[:, :, i], wj)
        vh = np.sqrt(amean(vt[:, :, i]**2 + vp[:, :, i]**2, wj))
        M = vh/csm
        rv2 = amean(rh*(vt[:, :, i]**2 + vp[:, :, i]**2), wj)/pm
        print('%6.3f | %9.4f %9.5f %9.5f %9.3e | %9.4f %9.4f %9.4f | '
              '%8.4f %8.4f %8.4f'
              % (R['r'], np.sqrt(amean(dlr**2, wj)), np.sqrt(amean(dlp**2, wj)),
                 np.sqrt(amean(dlt**2, wj)), bm,
                 np.sqrt(amean(pterm**2, wj)), np.sqrt(amean(tterm**2, wj)),
                 np.sqrt(amean(pred**2, wj)), M, M*M/bm, rv2), file=fh)
        # correlation check of the identity
        cc = (amean((pred - amean(pred, wj))*(dlr - amean(dlr, wj)), wj)
              / np.sqrt(amean((pred - amean(pred, wj))**2, wj)
                        * amean((dlr - amean(dlr, wj))**2, wj)))
        print('        identity check: corr(pred, dlnrho) = %.5f   slope = %.4f'
              % (cc, amean(pred*dlr, wj)/amean(dlr*dlr, wj)), file=fh)

    # --- 4. opacity slopes, resolution, kappa_eff ----------------------------
    print('\n--- 4. opacity slopes, resolution, effective opacity ---', file=fh)
    print('%6s | %9s %9s %9s %9s | %8s %8s %8s | %9s %9s %9s'
          % ('r/R', 'dlnk/dlnr', 'dlnk/dlnT', 'k(ch)/k(m)', 'k(cl)/k(m)',
             'w_ch/cel', 'dtau_cel', 'dtau_ch', '<Fcol>/Fr', 'Fmean/Fr',
             'keff/kap'), file=fh)
    for R in rows:
        k = R['k']
        i = sl[k]
        km = kappa(R['Tm'], R['rm'])
        # per-column flux with the local radial gradient
        Tm3 = [T[:, :, sl[k+dd]] for dd in (-1, 0, 1)]
        dTdr3 = (Tm3[2] - Tm3[0])/(rc[k+1] - rc[k-1])
        Fcol = -(4.0*ARAD*CLIGHT/3.0)*T[:, :, i]**3/(kap[:, :, i]*rho[:, :, i])*dTdr3
        Fcm = amean(Fcol, wj)
        Tms = [amean(x, wj) for x in Tm3]
        dTmdr = (Tms[2] - Tms[0])/(rc[k+1] - rc[k-1])
        Fmean = -(4.0*ARAD*CLIGHT/3.0)*R['Tm']**3/(km*R['rm'])*dTmdr
        Freq = LSTAR/(4*np.pi*rc[k]**2)
        dtau_cel = amean(kap[:, :, i]*rho[:, :, i], wj)*dr[k]
        print('%6.3f | %9.4f %9.4f %9.4f %9.4f | %8.2f %8.3f %8.3f | '
              '%9.4f %9.4f %9.4f'
              % (R['r'], dlnk_dlnrho(R['Tm'], R['rm']), dlnk_dlnT(R['Tm'], R['rm']),
                 R['klo']/km, R['khi']/km,
                 R['Lcell'], dtau_cel, R['klo']*R['rlo']*R['dl_th'],
                 Fcm/Freq, Fmean/Freq, Fmean/Fcm), file=fh)
    # --- 5. lateral (transverse-diffusion) flux vs radial flux ---------------
    print('\n--- 5. lateral radiative flux carried by the transverse diffusion '
          'operator ---', file=fh)
    print('%6s | %10s %10s %9s | %10s %10s %9s | %9s'
          % ('r/R', '<|Fperp|>', 'Freq', 'ratio', '<|Fperp|>ch', 'F_free',
             'Fp/Ffree', 'tsmooth/s'), file=fh)
    for R in rows:
        k, i = R['k'], sl[R['k']]
        Tk, rh, kk = T[:, :, i], rho[:, :, i], kap[:, :, i]
        gth = (np.roll(Tk, -1, 1) - np.roll(Tk, 1, 1))/(2*R['dl_th'])
        gph = (np.roll(Tk, -1, 0) - np.roll(Tk, 1, 0))/(2*R['dl_ph'])
        gmag = np.sqrt(gth**2 + gph**2)
        Fp = (4.0*ARAD*CLIGHT/3.0)*Tk**3/(kk*rh)*gmag
        Freq = LSTAR/(4*np.pi*rc[k]**2)
        Ffree = CLIGHT*ARAD*Tk**4        # c E_rad (free-streaming cap)
        low = rh <= np.percentile(rh, 10)
        print('%6.3f | %10.3e %10.3e %9.4f | %10.3e %10.3e %9.4f | %9.3e'
              % (R['r'], amean(Fp, wj), Freq, amean(Fp, wj)/Freq,
                 Fp[low].mean(), amean(Ffree, wj), amean(Fp, wj)/amean(Ffree, wj),
                 3.0*R['klo']*R['rlo']*R['L']**2/CLIGHT), file=fh)
    fh.flush()


if __name__ == '__main__':
    ids = [int(a) for a in sys.argv[1:]] or [14, 20]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'out_chan.txt'), 'w') as fh:
        for i in ids:
            run(i, fh)
    print(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'out_chan.txt')).read())
