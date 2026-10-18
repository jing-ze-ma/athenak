#!/usr/bin/env python3
"""bsg_puls.py IC EOSDUMP ROSS GM_IN L RIN RTOP [RUN ...]: cheap physics-vs-closure diagnostics
(10-16).
(b) linear ADIABATIC radial pulsation modes of the IC column between the closed wall r_in
    (xi = 0) and RTOP (outer bound: delta P/P + (4 + w^2 r^3/(G m)) xi = 0), Gamma_1 of the
    code EOS + radiation, gravity G m(r) of the column (Unno et al. radial form, y1 = xi =
    dr/r, y2 = Lagrangian dP/P):
        dy1/dlnr = -3 y1 - y2/Gamma1,   dy2/dlnr = V (y2 + (4 + w^2 r^3/(G m)) y1),  V = rho g r/P;
    the thermal time of the layer above tau 30 (int c_v rho T dV / L) and its KH ratio;
(c) Gamma_rad = kappa F_r/(c g) and Gamma_L = kappa L/(4 pi c G m) along the column vs the radii
    where the runs' shell-mean |v_r| peaks; the oscillation period from <v_r>(t) zero crossings."""
import glob
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from icstd.micro import Star, CL, RSUN, G_CODE, A_RAD_  # noqa: E402
from icstd.seam import read_col, tau_from_top  # noqa: E402


def main():
    ic, eosf, ross, gmin, L, rin, rtop = sys.argv[1:8]
    gmin, L, rin, rtop = float(gmin), float(L), float(rin), float(rtop)
    runs = sys.argv[8:]
    r, rho, eint, Fr, E, T, fm = read_col(ic)
    st = Star(gmin, L, 1.0, eosf, ross, -21.0, 2.0, 0.0, 1e-4, 1e-3, verbose=False)
    k = (r >= rin*0.999) & (r <= rtop)
    r, rho, T, Fr, fm = r[k], rho[k], T[k], Fr[k], fm[k]
    P, pg, chir, chit, cp, dlt, nad = st.thermo(T, rho)
    cv = cp - P*chit**2/(rho*T*chir)
    g1 = chir + chit**2*P/(rho*T*cv)
    m = gmin/G_CODE + np.concatenate([[0], np.cumsum(0.5*np.diff(r)*4*np.pi*(
        r[1:]**2*rho[1:] + r[:-1]**2*rho[:-1]))])
    g = G_CODE*m/r**2
    V = rho*g*r/P
    lr = np.log(r)
    fV = lambda x: np.interp(x, lr, V)          # noqa: E731
    fg1 = lambda x: np.interp(x, lr, g1)        # noqa: E731
    fs = lambda x: np.interp(x, lr, r**3/(G_CODE*m))   # noqa: E731

    def resid(w):
        def rhs(x, y):
            return [-3*y[0] - y[1]/fg1(x), fV(x)*(y[1] + (4 + w*w*fs(x))*y[0])]
        s = solve_ivp(rhs, [lr[0], lr[-1]], [0.0, 1.0], method='LSODA', rtol=1e-8, atol=1e-12)
        y = s.y[:, -1]
        return (y[1] + (4 + w*w*fs(lr[-1]))*y[0])/np.hypot(*y)
    tdyn = np.sqrt(r[-1]**3/(G_CODE*m[-1]))
    ws = np.linspace(0.2, 12.0, 120)/tdyn
    fsv = [resid(w) for w in ws]
    modes = []
    for i in range(len(ws) - 1):
        if np.sign(fsv[i]) != np.sign(fsv[i+1]):
            modes.append(brentq(resid, ws[i], ws[i+1], xtol=1e-6*ws[i]))
    print('# IC %s: %.3f-%.3f Rsun (closed wall at r_in, outer bound at %.3f Rsun), t_dyn %.3e s'
          % (ic, r[0]/RSUN, r[-1]/RSUN, r[-1]/RSUN, tdyn))
    for n, w in enumerate(modes[:4]):
        print('  adiabatic radial mode n=%d: period %.4e s = %.3f d' % (n, 2*np.pi/w, 2*np.pi/w/86400))
    # sound crossing, thermal time above tau 30
    cs = np.sqrt(g1*P/rho)
    print('  sound crossing r_in..top: %.3e s = %.3f d' % (np.trapz(1/cs, r), np.trapz(1/cs, r)/86400))
    kap = st.kap(T, rho)
    tau = tau_from_top(r, rho, kap, rtop)
    for tc in (30.0, 100.0, 1000.0):
        j = tau <= tc
        U = np.trapz((rho*(cv*T) + A_RAD_*T**4)[j]*4*np.pi*r[j]**2, r[j])
        print('  layer above tau %g (r > %.3f Rsun): thermal time U/L %.3e s = %.3f d'
              % (tc, r[j].min()/RSUN, U/L, U/L/86400))
    GamL = kap*L/(4*np.pi*CL*G_CODE*m)
    Gamr = kap*Fr*4*np.pi*r**2/(4*np.pi*CL*G_CODE*m)
    # acoustic cut-off of the photosphere/atmosphere: w_c = c_s/(2 H_P), H_P = P/(rho g_eff)
    iph = int(np.argmin(abs(tau - 2/3)))
    for lab, ge in (('g', g), ('g (1 - Gamma_rad)', g*(1 - np.minimum(Gamr, 0.99)))):
        Hp = P/(rho*ge)
        wc = cs/(2*Hp)
        print('  acoustic cut-off period 2 pi/w_c at R_ph (%s): %.3e s = %.3f d; range over tau 0.1-3: '
              '%.3f-%.3f d' % (lab, 2*np.pi/wc[iph], 2*np.pi/wc[iph]/86400,
                               (2*np.pi/wc[(tau > 0.1) & (tau < 3)]).min()/86400,
                               (2*np.pi/wc[(tau > 0.1) & (tau < 3)]).max()/86400))
    print('  r[Rsun]  tau      Gamma_L  Gamma_rad  F_MLT/F  T[K]')
    for x in np.linspace(0.85, 1.0, 16)*r[np.argmin(abs(tau - 2/3))]:
        i = np.argmin(abs(r - x))
        print('  %7.3f  %8.2e  %6.3f   %6.3f    %5.3f  %.3e' % (r[i]/RSUN, tau[i], GamL[i], Gamr[i],
                                                               fm[i], T[i]))
    for run in runs:
        import bin_convert as bc
        fs_ = sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin')))
        ts, vr = [], []
        for f in fs_:
            h = bc.read_binary_as_athdf(f)
            x2f = np.asarray(h['x2f'], float)
            dom = (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :, None] * \
                np.diff(np.asarray(h['x3f'], float))[:, None, None]
            d = np.asarray(h['dens'], float)
            w = np.broadcast_to(dom, d.shape)
            v = (d*np.asarray(h['velx'], float)*w).sum(axis=(0, 1))/(d*w).sum(axis=(0, 1))/1e5
            ts.append(float(h['Time']))
            vr.append(v)
            rr = np.asarray(h['x1v'], float)
        ts, vr = np.array(ts), np.array(vr)
        rmsr = np.sqrt((vr**2).mean(axis=0))
        ipk = int(np.argmax(rmsr))
        print('# run %s: rms_t <v_r> peaks at r %.3f Rsun (%.2f km/s); profile:' % (run, rr[ipk]/RSUN,
                                                                                rmsr[ipk]))
        for x in np.linspace(0.8, 1.15, 15)*rr[ipk]:
            i = np.argmin(abs(rr - x))
            j = np.argmin(abs(r - rr[i])) if rr[i] <= r[-1] else None
            print('   r %7.3f  rms<v_r> %6.3f km/s  Gamma_L %s  Gamma_rad %s'
                  % (rr[i]/RSUN, rmsr[i], '%.3f' % GamL[j] if j is not None else '-',
                     '%.3f' % Gamr[j] if j is not None else '-'))
        # period from zero crossings at the peak radius, per third of the run
        v = vr[:, ipk]
        zc = ts[1:][np.sign(v[1:]) != np.sign(v[:-1])]
        nseg = max(3, int(ts[-1]/86400/1.0 + 0.5))
        for a, b in [(q/nseg, (q + 1)/nseg if q < nseg - 1 else 1.01) for q in range(nseg)]:
            z = zc[(zc >= a*ts[-1]) & (zc < b*ts[-1])]
            per = 2*np.mean(np.diff(z))/86400 if z.size > 2 else np.nan
            amp = np.abs(v[(ts >= a*ts[-1]) & (ts < b*ts[-1])]).max()
            print('   t %.2f-%.2f d: period %.3f d, max |<v_r>| %.2f km/s at %.3f Rsun'
                  % (a*ts[-1]/86400, min(b, 1)*ts[-1]/86400, per, amp, rr[ipk]/RSUN))


if __name__ == '__main__':
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                                    'vis', 'python'))
    main()
