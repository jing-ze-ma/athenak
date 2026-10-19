#!/usr/bin/env python3
"""S3b/S4 of the accretor RHD port: the PHYSICAL Plaskett L1 stream (estimates, no run).
Usage: stream_physics.py EOS_DUMP ROSS_TABLE PLANCK_TABLE
(1) env13's isothermal stream (26.2 kK, Mdot 1e-4, pressure-supported widths) with the real
    Z 0.02 opacities: optical depths, diffusion time vs flight time, thin wings.
(2) the L1 state from the donor's own atmosphere (optically thick overflow, Kolb & Ritter
    1990 form calibrated to Ryu et al. 2025): Mdot = Q 2 pi/(Omega^2 sqrt(BC)) int c_T dP over
    the donor's grey atmosphere (EOS + Rosseland table, gas + radiation pressure).
(3) the core carried ADIABATICALLY (gas + trapped radiation) from L1 to r_out, with the
    pressure-supported transverse widths there, solved self-consistently with Mdot.
(4) donor irradiation of the gainer's facing side.
"""
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.special import erfc

kB, mH, G = 1.380649e-16, 1.6726e-24, 6.674e-8
A_RAD, C_L, SIG = 7.5657e-15, 2.99792458e10, 5.6704e-5
Msun, Rsun, Lsun = 1.989e33, 6.957e10, 3.828e33
yr = 3.15576e7


def read_dump(fn):
    with open(fn) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    d = np.loadtxt(fn, comments='#')
    x = xmin + dx*np.arange(int(nx)); y = ymin + dy*np.arange(int(ny))
    return x, y, d[:, 0].reshape(int(ny), int(nx)), d[:, 1].reshape(int(ny), int(nx))


def read_tab(fn):
    hdr = None
    with open(fn) as fh:
        for ln in fh:
            if not ln.startswith('#'):
                break
            w = ln.split()
            if len(w) == 7 and w[1].isdigit():
                hdr = [float(v) for v in w[1:]]
    nT, nD, lT0, dlT, lD0, dlD = hdr
    v = np.loadtxt(fn, comments='#').reshape(int(nT), int(nD))
    return RectBivariateSpline(lT0 + dlT*np.arange(int(nT)), lD0 + dlD*np.arange(int(nD)), v, kx=1, ky=1)


xe, ye, le, lp = read_dump(sys.argv[1])
sp_e = RectBivariateSpline(ye, xe, le, kx=3, ky=3)
sp_p = RectBivariateSpline(ye, xe, lp, kx=3, ky=3)
kR_s, kP_s = read_tab(sys.argv[2]), read_tab(sys.argv[3])
kR = lambda r, T: 10**kR_s(np.log10(T), np.log10(r), grid=False)
kP = lambda r, T: 10**kP_s(np.log10(T), np.log10(r), grid=False)
pg = lambda r, T: r*10**sp_p(np.log10(T), np.log10(r), grid=False)
eg = lambda r, T: r*10**sp_e(np.log10(T), np.log10(r), grid=False)   # erg/cc
pr = lambda T: A_RAD*T**4/3.0


def rho_of(p, T):
    return 10**brentq(lambda ld: ld + sp_p(np.log10(T), ld, grid=False) - np.log10(p), xe[0], xe[-1])


# binary (code constants of ry_per_accretor / plaskett_setup_1007.py)
Ma, Md, a = 16.0, 18.0, 32.5575
Om = 1.9708e-5                          # 1/s
B, C = 6.994, 7.994                     # L1 curvature (Omega^2 units)
Pnn, Pzz = 4.135, 9.218                 # transverse curvature at the r_out crossing (Omega^2)
vout = 134.87e5                         # |v| at r_out, cm/s
t_fl = 0.0507*6.957e5                   # L1 -> r_out flight, s
t_imp = 4.6*Rsun/(0.5*(134.87e5 + 395.5e5))   # r_out -> photosphere, ~ s
Mdot = 1.0e-4*Msun/yr
Td, Rd, Ld = 26200.0, 12.8, 10**4.84
Tg, Rg = 28300.0, 9.0
print('# Mdot %.4e g/s, flight L1->r_out %.3e s, r_out->impact ~%.3e s' % (Mdot, t_fl, t_imp))


def report(tag, rho, T, sy, sz, tflight):
    k, kp = kR(rho, T), kP(rho, T)
    u = eg(rho, T) + A_RAD*T**4
    ty, tz = k*rho*sy*Rsun, k*rho*sz*Rsun
    s = min(sy, sz)*Rsun
    tdiff = 3.0*k*rho*s*s*u/(C_L*A_RAD*T**4)
    tthin = eg(rho, T)/(4.0*kp*rho*SIG*T**4)      # optically thin cooling time of the gas
    print('%-34s rho %.3e T %7.0f  kR %.3f kP %.3f  Prad/Pgas %.4f  tau_perp %.3g tau_z %.3g  '
          't_diff %.3g s = %.3g flights; thin t_cool %.3g s' % (tag, rho, T, k, kp, pr(T)/pg(rho, T),
          ty, tz, tdiff, tdiff/tflight, tthin))
    return tdiff


# ---------------------------------------------------------------- (1) env13 stream
print('\n## (1) env13 stream (isothermal 26.2 kK) with the Z 0.02 tables')
report('L1 (Ryu analytic x0.721)', 6.20e-7, Td, 0.4801, 0.4819, t_fl)
report('r_out peak (env13 widths)', 5.11e-8, Td, 0.670, 0.449, t_fl)
# thin wing: where the column outside offset x (in the plane) reaches tau 1
rho0, sg = 5.11e-8, 0.670
for xs in np.arange(0.0, 6.01, 0.5):
    rr = rho0*np.exp(-0.5*xs*xs)
    tau = kR(rr, Td)*rho0*sg*Rsun*np.sqrt(np.pi/2)*erfc(xs/np.sqrt(2))
    print('  r_out wing x = %.1f sigma: rho %.3e  tau(>x) %.3g  Prad/Pgas %.3f' % (xs, rr, tau, pr(Td)/pg(rr, Td)))

# ---------------------------------------------------------------- (2) L1 from the donor atmosphere
print('\n## (2) L1 state from the donor atmosphere (grey, Eddington T(tau), EOS + Rosseland, gas + rad)')
Q_iso, Q_ad = 0.721, 0.649
gd = G*Md*Msun/(Rd*Rsun)**2
fac = 2*np.pi/(Om**2*np.sqrt(B*C))
for gfac in (1.0, 0.3):
    g = gd*gfac
    # integrate the atmosphere in tau: dP_gas/dtau = g/kappa - dP_rad/dtau, P_rad = a T^4/3
    taus = np.logspace(-4, 6, 4000)
    Tt = lambda t: (0.75*Td**4*(t + 2.0/3.0))**0.25
    def rhs(lt, y):
        t = np.exp(lt); T = Tt(t)
        p = np.exp(y[0]); r = rho_of(p, T)
        dprad = A_RAD*0.75*Td**4/3.0                # dP_rad/dtau (Eddington)
        return [t*(g/kR(r, T) - dprad)/p]
    # start: P_gas small at tau 1e-4
    sol = solve_ivp(rhs, (np.log(1e-4), np.log(1e6)), [np.log(1e-2)], t_eval=np.log(taus),
                    rtol=1e-7, atol=1e-10)
    tt = np.exp(sol.t); P = np.exp(sol.y[0]); T = Tt(tt)
    rho = np.array([rho_of(p, x) for p, x in zip(P, T)])
    cT = np.sqrt(P/rho)
    # Mdot(P_L1) = Q fac int c_T dP (gas pressure)
    cum = np.concatenate([[0.0], np.cumsum(0.5*(cT[1:] + cT[:-1])*np.diff(P))])
    for Q, nm in ((Q_iso, 'iso 0.721'), (Q_ad, 'adiab 0.649')):
        md = Q*fac*cum
        i = np.searchsorted(md, Mdot)
        f = (Mdot - md[i-1])/(md[i] - md[i-1])
        lt = np.log(tt[i-1]) + f*np.log(tt[i]/tt[i-1]); tL = np.exp(lt)
        TL = Tt(tL); PL = np.exp(np.log(P[i-1]) + f*np.log(P[i]/P[i-1])); rL = rho_of(PL, TL)
        cL = np.sqrt(PL/rL)
        k = kR(rL, TL)
        print('g %.3g (x%.1f), Q %s: L1 at tau %.4g: T %.0f K, rho %.3e, P_gas %.3e, Prad/Pgas %.3f, '
              'c_T %.2f km/s, sigma_y %.3f sigma_z %.3f Rsun, kappa %.3f' % (g, gfac, nm, tL, TL, rL, PL,
              pr(TL)/PL, cL/1e5, 0.932*cL/(Om*np.sqrt(B))/Rsun, cL/(Om*np.sqrt(C))/Rsun, k))
        if gfac == 1.0 and Q == Q_iso:
            L1 = (rL, TL, cL)
    i = np.searchsorted(tt, 2.0/3.0)
    print('   donor photosphere (tau 2/3, g x%.1f): rho %.3e g/cc, P_gas %.3e' % (gfac, rho[i], P[i]))

# ---------------------------------------------------------------- (3) adiabatic core to r_out
print('\n## (3) core carried adiabatically (gas + trapped radiation) from L1 to r_out')
rL, TL, cL = L1
sy, sz = 0.932*cL/(Om*np.sqrt(B))/Rsun, cL/(Om*np.sqrt(C))/Rsun
report('L1 (donor atmosphere, iso Q)', rL, TL, sy, sz, t_fl)


def adiabat(r0, T0, r1):
    # d(e_tot/rho) = -P_tot d(1/rho): integrate T(rho) along ln rho
    def dT(lr, y):
        r = np.exp(lr); T = y[0]
        h = 1e-4
        etot = lambda rr, TT: (eg(rr, TT) + A_RAD*TT**4)/rr
        dedT = (etot(r, T*(1 + h)) - etot(r, T*(1 - h)))/(2*h*T)
        dedr = (etot(r*(1 + h), T) - etot(r*(1 - h), T))/(2*h)     # d e/d ln rho
        P = pg(r, T) + pr(T)
        return [(P/r - dedr)/dedT]
    s = solve_ivp(dT, (np.log(r0), np.log(r1)), [T0], rtol=1e-8)
    return s.y[0][-1]


rho_pk = 5.11e-8
for it in range(60):
    T1 = adiabat(rL, TL, rho_pk)
    c1 = np.sqrt((pg(rho_pk, T1) + pr(T1))/rho_pk)
    sp_, sz_ = c1/(Om*np.sqrt(Pnn)), c1/(Om*np.sqrt(Pzz))
    new = Mdot/(vout*2*np.pi*sp_*sz_)
    if abs(new/rho_pk - 1) < 1e-8:
        break
    rho_pk = np.sqrt(new*rho_pk)
print('adiabatic r_out: rho_peak %.4e g/cc (%.3f code), T %.0f K, c_iso(P_tot) %.2f km/s, sigma_perp %.3f, '
      'sigma_z %.3f Rsun, arc width (/cos 24.3) %.3f' % (rho_pk, rho_pk/5.11e-8, T1, c1/1e5, sp_/Rsun, sz_/Rsun,
      sp_/Rsun/np.cos(np.radians(24.3))))
report('r_out peak (adiabatic core)', rho_pk, T1, sp_/Rsun, sz_/Rsun, t_fl)
# radiative equilibrium temperature of an optically thin parcel at r_out (both stars, dilution)
def W(R, d):
    return 0.5*(1 - np.sqrt(1 - (R/d)**2))
rg, rd_ = 13.5015, np.hypot(a - 13.5015*np.cos(np.radians(3.8)), 13.5015*np.sin(np.radians(3.8)))
Tthin = (W(Rg, rg)*Tg**4 + W(Rd, rd_)*Td**4)**0.25
Tthin_g = (W(Rg, rg)*Tg**4)**0.25
print('thin-parcel radiative equilibrium T at r_out (grey, both stars) %.0f K; gainer only (vacuum r_out) %.0f K'
      % (Tthin, Tthin_g))

Teq = Tthin   # the core cannot radiate below the irradiated surface temperature
# one-zone core along the flight: adiabatic work + radiative diffusion loss of the trapped
# radiation, d(e_tot/rho)/dt = -P_tot d(1/rho)/dt - (a T^4/rho)/t_diff, t_diff = 3 kR rho s^2/c
# (s = the smaller transverse sigma); rho(t) = Mdot/(v 2 pi s_perp s_z) with pressure-supported
# widths s = c_iso/(Omega sqrt(curvature)), curvature and |v| linear in t between L1 (B, C,
# c_L) and r_out (Pnn, Pzz, 134.9 km/s).  A model: brackets, not a prediction.
def onezone(r0, T0, nstep=20000):
    dt = t_fl/nstep
    T, r = T0, r0
    etot = lambda rr, TT: (eg(rr, TT) + A_RAD*TT**4)/rr
    for n in range(nstep):
        f = (n + 1)/nstep
        v = (1 - f)*np.sqrt((pg(r, T) + pr(T))/r) + f*vout
        cp, cz = (1 - f)*B/0.932**2 + f*Pnn, (1 - f)*C + f*Pzz
        c2 = (pg(r, T) + pr(T))/r
        rn = Mdot/(v*2*np.pi*c2/(Om**2*np.sqrt(cp*cz)))
        rn = r*min(max(rn/r, 0.995), 1.005)    # relax rho toward the width-consistent value
        s = np.sqrt(c2)/(Om*np.sqrt(max(cp, cz)))
        td = 3.0*kR(r, T)*r*s*s/C_L
        loss = min(A_RAD*(T**4 - Teq**4)/r*dt/td, 0.5*A_RAD*T**4/r)
        e1 = etot(r, T) - (pg(r, T) + pr(T))*(1/rn - 1/r) - loss
        T = brentq(lambda x: etot(rn, x) - e1, 1.0e3, 3.0e6)
        r = rn
    return r, T, np.sqrt((pg(r, T) + pr(T))/r)
r2, T2, c2_ = onezone(L1[0], L1[1])
print('one-zone (adiabatic + diffusion) r_out: rho_peak %.4e g/cc (%.3f code), T %.0f K, c_iso %.2f km/s, '
      'sigma_perp %.3f sigma_z %.3f Rsun, arc width %.3f, Prad/Pgas %.3f' % (r2, r2/5.11e-8, T2, c2_/1e5,
      c2_/(Om*np.sqrt(Pnn))/Rsun, c2_/(Om*np.sqrt(Pzz))/Rsun, c2_/(Om*np.sqrt(Pnn))/Rsun/np.cos(np.radians(24.3)),
      pr(T2)/pg(r2, T2)))
report('r_out peak (one-zone)', r2, T2, c2_/(Om*np.sqrt(Pnn))/Rsun, c2_/(Om*np.sqrt(Pzz))/Rsun, t_fl)
# isothermal-at-26.2kK widths but Mdot with the hot L1 is the same: env13 r_out row for contrast
# ---------------------------------------------------------------- (4) donor irradiation
print('\n## (4) donor irradiation of the gainer (point-source dilution of the donor disc)')
for ph, rp_ in ((0.0, 9.4847), (14.0, 9.4316), (45.0, 9.1466), (90.0, 9.0013)):
    x, y = rp_*np.cos(np.radians(ph)), rp_*np.sin(np.radians(ph))
    d = np.hypot(a - x, y)
    mu = ((a - x)*np.cos(np.radians(ph)) - y*np.sin(np.radians(ph)))/d    # cos of incidence (radial normal)
    Fin = SIG*Td**4*(Rd/d)**2*max(mu, 0.0)
    print('phi %5.1f: d %.2f Rsun, cos(inc) %.3f, F_irr/F_gainer(28.3 kK) = %.3f' % (ph, d, mu, Fin/(SIG*Tg**4)))
