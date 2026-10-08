"""sanity checks + T_d(a,r), kappa_pr(a), S=1 temperatures. usage: venv/bin/python checks.py"""
import numpy as np, time
import grainlib as gl
ch = gl.Chem()
for sn in gl.STARS:
    st = gl.Star(sn)
    print('== %s  R=%.3e vesc=%.1f km/s v_con=%.2f kE=%.3f Mdot0=%.2e Msun/yr' %
          (sn, st.R, st.vesc/1e5, st.v_con/1e5, st.kE, st.mdot0/gl.Msun*gl.YR))
    rr = np.array([2, 3, 4, 5, 7, 10]) * st.R
    print(' r/R', rr/st.R, '\n W', st.W(rr), '\n rho_FT', st.rho_ft(rr), '\n BB Td', st.Teff*st.W(rr)**.25)
    for var in ['lowk', 'pure', 'Fe3e-4', 'Fe1e-3']:
        m = gl.Model(sn, var, ch)
        for a in [1e-7, 1e-5, 3e-5, 1e-4]:
            Td = [m.Td_at(a, x) for x in rr]
            # S at Td with FT density (nH), fcond=0
            lS = [ch.lnS(t, st.rho_ft(x)/ch.mH, 0.) for t, x in zip(Td, rr)]
            print(' %-6s a=%5.0f nm Td=%s  lnS=%s' % (var, a*1e7, np.round(Td).astype(int), np.round(lS, 1)))
        if sn == 'golden16':
            aa = np.array([1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4])
            print('   kpr_dust [cm2/g dust] a(um)=', aa*1e4, np.round([m.kpr_dust(x) for x in aa]),
                  ' Qpr', np.round(np.interp(np.log(aa), m.op.la, m.qpr), 4))
# T(S=1) vs nH for a grain: solve lnS(T, nH)=0
for nH in [1e8, 1e9, 1e10, 1e11]:
    T = np.arange(700, 1500, 1.)
    l = np.array([ch.lnS(t, nH, 0.) for t in T])
    print('nH=%.0e  T(S=1)=%.0f K' % (nH, T[np.argmin(abs(l))]))

# ---- growth time scale: analytic a = a0 + V_mon J t vs numerical (fixed r, S>>1, tiny f)
from scipy.integrate import solve_ivp
m = gl.Model('golden16', 'pure', ch)
st = m.st
r = 4 * st.R; nH = st.rho_ft(r) / ch.mH; Tg = 800.; al = 0.1; xs = 1e-16
vbar = np.sqrt(8 * gl.kB * Tg / (np.pi * gl.M_MG))
rate = gl.V_MON * al * nH * ch.eMg * vbar / 8
sol = solve_ivp(lambda t, y: [m.dadt(y[0], r, nH, Tg, xs, al)[0]], (0, 1e9), [gl.A_SEED],
                rtol=1e-8, t_eval=[1e8, 1e9])
print('growth @4R golden16 nH=%.2e: analytic a(1e8,1e9 s)=%s  numerical=%s  t(0.3um)=%.2e s' %
      (nH, gl.A_SEED + rate * np.array([1e8, 1e9]), sol.y[0], 3e-5 / rate))
# ---- gas-grain heating vs radiative heating, worst case (m20 2R, Tg=1700, 1 nm)
for sn in gl.STARS:
    mm = gl.Model(sn, 'pure', ch)
    s = mm.st
    for rr_ in [2, 3]:
        r = rr_ * s.R; nH = s.rho_ft(r) / ch.mH
        for a in [1e-7, 1e-5]:
            ia = np.argmin(abs(mm.op.la - np.log(a)))
            Td = mm.Td_at(a, r)
            qT = np.interp(np.log(s.Teff), mm.op.lT, mm.op.QaT[ia])
            Prad = 4 * s.W(r) * 5.6704e-5 * s.Teff**4 * qT
            vH2 = np.sqrt(8 * gl.kB * 1700 / (np.pi * 2 * 1.008 * gl.mu_))
            Pcol = 0.586 * nH * vH2 / 4 * 2 * gl.kB * (1700 - Td)
            print('%s r=%dR a=%.0fnm Td=%.0f: gas-grain/radiative heating = %.1e' %
                  (sn, rr_, a * 1e7, Td, Pcol / Prad))
# ---- Mg photoionisation by the diluted Planck field (radiative Saha, Te = 1000 K)
me, hP, chi = 9.109e-28, 6.626e-27, 7.646 * 1.602e-12
for sn in gl.STARS:
    s = gl.Star(sn)
    for rr_ in [3, 4]:
        r = rr_ * s.R; nH = s.rho_ft(r) / ch.mH
        TR = s.Teff
        sah = s.W(r) * np.sqrt(1000 / TR) * 2 * (2 / 1) * (2 * np.pi * me * gl.kB * TR / hP**2)**1.5 \
            * np.exp(-chi / (gl.kB * TR))
        ne = 1e-4 * nH
        print('%s r=%dR nH=%.1e: Mg+/Mg (Planck Teff, ne=1e-4 nH) = %.1e' % (sn, rr_, nH, sah / ne))
