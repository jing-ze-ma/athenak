#!/usr/bin/env python3
"""S5 (accretor-rhd-1008, user 10-09): stream_physics.py with the 3.81-d geometry
(plaskett_setup_1009.py output, rhd/s5/setup/setup_1009.out) and Ryu+2025's ADIABATIC L1 rate
factor 0.649 (optically thick overflow) for the L1 state that seeds the one-zone core.
The one-zone r_out state (section 3, last lines) is the S5 stream state; section (5) adds the
window quantities the pgen needs (code units, rho_u = 5.11e-8 g/cc) and the injected radiation.

Original: S3b/S4 of the accretor RHD port: the PHYSICAL Plaskett L1 stream (estimates, no run).
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


# binary: 3.81 d (plaskett_setup_1009.py, NX1=640, PDAY=3.81; env13 3.69-d values in comments)
Ma, Md, a = 16.0, 18.0, 33.2596         # env13 32.5575
Om = 1.908715e-5                        # 1/s (env13 1.9708e-5)
B, C = 6.994, 7.994                     # L1 curvature (Omega^2 units; depend on q only)
# The ballistic launch speed is the L1 overflow speed of the hot L1 state of section (2)
# (c_T = sqrt(P_gas/rho) = 45.45 km/s, the speed in the Mdot integral), not the donor-Teff
# 18.68 km/s of env13: problem/t_don = 155115 K (mu 0.62) in the pgen.  Central trajectory
# for that launch (rhd/s5/setup/traj_launch.py; env13-launch values in brackets):
Pnn, Pzz = 3.6132, 9.222                # transverse curvature at the r_out crossing [4.128, 9.218]
vout = 140.27e5                         # |v| at r_out, cm/s [133.47]
cosa = 126.46/140.27                    # |v_r|/|v| at r_out [121.64/133.47]
r_out = 13.7926                         # Rsun (0.85 d_L1, d_L1 16.22664; env13 13.5015)
phi_s = 3.410                           # deg, central trajectory at r_out [3.798]
t_fl = 0.0348*6.957e5                   # L1 -> r_out flight, s [0.0522 code]
t_imp = (0.0550 - 0.0348)*6.957e5       # r_out -> photosphere (r 9.39, phi 15.07 deg, |v| 410.7), s
Mdot = 1.0e-4*Msun/yr
Td, Rd = 26200.0, 12.9433               # donor Teff (HR figure, kept) and R = R_L (Eggleton, 3.81 d)
Ld = 4*np.pi*(Rd*Rsun)**2*SIG*Td**4/Lsun
Tg, Rg = 28300.0, 9.0
RHO_U = 5.11e-8                         # code density unit (g/cc), unchanged from env13/S3a
Q_L1 = 0.649                            # Ryu+2025 adiabatic (optically thick overflow), user 10-09
print('# donor R_L %.4f Rsun, L %.4g Lsun (log %.3f) at Teff %.0f K' % (Rd, Ld, np.log10(Ld), Td))
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
        if gfac == 1.0 and Q == Q_L1:
            L1 = (rL, TL, cL)
    i = np.searchsorted(tt, 2.0/3.0)
    print('   donor photosphere (tau 2/3, g x%.1f): rho %.3e g/cc, P_gas %.3e' % (gfac, rho[i], P[i]))


# ======================================================================================
# stream_path_1015 (accretor-bigdom-1015, user 10-15): the stream state ALONG the ballistic
# path from L1 (launch c_T(L1) = 45.45 km/s, as the pgen), for an injection nozzle between
# L1 and 13.79: the same one-zone core as section (3) of stream_physics_1009 (adiabatic gas +
# trapped radiation, diffusion loss toward T_eq), but with the true |v|(t) and the true
# transverse curvatures Phi_nn, Phi_zz along the path (not the linear L1 -> r_out ramp).
# Output per radius: |v|, v_r, v_phi, phi, rho_peak, T, Prad/Pgas, sound speeds (isothermal
# gas c_T, adiabatic gas-only c_g from the EOS, adiabatic gas + radiation c_tot), Mach numbers,
# sigma_perp and the arc width sigma_perp/cos(angle from radial) for problem/stream_width.
# Usage: python3 stream_path_1015.py EOS_DUMP ROSS PLANCK
print('\n## stream_path_1015: one-zone core along the true ballistic path from L1')
GU = G*Msun/(Rsun*1e10)
om = np.sqrt(GU*(Ma + Md)/a**3); aom = a*om
mu = Md/(Ma + Md); xa, xd = -mu, 1 - mu
gx = lambda x: -(1-mu)*(x-xa)/abs(x-xa)**3 - mu*(x-xd)/abs(x-xd)**3 + x
xl1 = brentq(gx, xa + 1e-3, xd - 1e-3)
def rhs(t, s):
    x, y, vx, vy = s; r1 = np.hypot(x-xa, y); r2 = np.hypot(x-xd, y)
    return [vx, vy, -(1-mu)*(x-xa)/r1**3 - mu*(x-xd)/r2**3 + x + 2*vy,
            -(1-mu)*y/r1**3 - mu*y/r2**3 + y - 2*vx]
rL, TL, cL = L1
ev = lambda t, s: np.hypot(s[0]-xa, s[1]) - 13.0/a
ev.terminal = True; ev.direction = -1
so = solve_ivp(rhs, [0, 6], [xl1 - 1e-4, 0, -cL/1e5/aom, 0], events=ev, rtol=1e-11, atol=1e-13,
               dense_output=True, max_step=1e-4)
tend = so.t_events[0][0]
GMA, GMD, XCM = GU*Ma, GU*Md, a*Md/(Ma + Md)
def roche3(x, y, z):
    return -GMA/np.sqrt(x*x+y*y+z*z) - GMD/np.sqrt((x-a)**2+y*y+z*z) - 0.5*om**2*((x-XCM)**2+y*y)
def path(t):
    s = so.sol(t); dx, y = s[0]-xa, s[1]; r = np.hypot(dx, y)
    vr = (s[2]*dx + s[3]*y)/r*aom; vp = (-s[2]*y + s[3]*dx)/r*aom
    xo, yo = dx*a, y*a
    vh = np.array([s[2], s[3]]); vh /= np.hypot(*vh); nh = np.array([-vh[1], vh[0]]); h = 1e-3
    f0 = roche3(xo, yo, 0)
    pnn = (roche3(xo+h*nh[0], yo+h*nh[1], 0) - 2*f0 + roche3(xo-h*nh[0], yo-h*nh[1], 0))/h**2/om**2
    pzz = (roche3(xo, yo, h) - 2*f0 + roche3(xo, yo, -h))/h**2/om**2
    return r*a, vr, vp, np.degrees(np.arctan2(y, dx)), pnn, pzz
etot = lambda rr, TT: (eg(rr, TT) + A_RAD*TT**4)/rr
def cad(r, T, rad=True):
    # adiabatic sound speed: dP/drho along the adiabat (gas only: e_gas, P_gas)
    h = 1e-4
    if rad:
        e = etot; P = lambda rr, TT: pg(rr, TT) + pr(TT)
    else:
        e = lambda rr, TT: eg(rr, TT)/rr; P = pg
    dedT = (e(r, T*(1+h)) - e(r, T*(1-h)))/(2*h*T)
    dedr = (e(r*(1+h), T) - e(r*(1-h), T))/(2*h*r)
    dTdr = (P(r, T)/r**2 - dedr)/dedT          # de = (P/rho^2) drho along the adiabat
    dPdr = ((P(r*(1+h), T) - P(r*(1-h), T))/(2*h*r) + (P(r, T*(1+h)) - P(r, T*(1-h)))/(2*h*T)*dTdr)
    return np.sqrt(dPdr)
# one-zone along the path: time steps of the flight; transverse curvature along the stream's
# normal (pnn) and z; widths s = c_iso/(Om sqrt(curv)); rho = Mdot/(|v| 2 pi s_perp s_z)
nst = 40000
T, r = TL, rL
rows = []
targets = [16.20, 16.10, 16.00, 15.75, 15.50, 15.25, 15.00, 14.75, 14.50, 14.25, 14.00, 13.7926]
ti = 0
tprev = 0.0
for n in range(1, nst + 1):
    t = tend*n/nst
    rr_, vr, vp, ph, pnn, pzz = path(t)
    v = np.hypot(vr, vp)*1e5
    # near L1 the in-plane normal curvature is B (Lubow-Shu); use max(pnn, small) for widths
    cpn = max(pnn, 0.5); cpz = max(pzz, 0.5)
    c2 = (pg(r, T) + pr(T))/r
    vv = max(v, np.sqrt(pg(r, T)/r))           # at the launch the flow speed is c_T
    rn = Mdot/(vv*2*np.pi*c2/(Om**2*np.sqrt(cpn*cpz)))
    rn = r*min(max(rn/r, 0.995), 1.005)
    s = np.sqrt(c2)/(Om*np.sqrt(max(cpn, cpz)))
    dt = (t - tprev)/om*6.957e5               # s (om is 1/code time)
    tprev = t
    Teq = 19239.0
    td = 3.0*kR(r, T)*r*s*s/C_L
    loss = min(A_RAD*(T**4 - Teq**4)/r*dt/td, 0.5*A_RAD*T**4/r)
    e1 = etot(r, T) - (pg(r, T) + pr(T))*(1/rn - 1/r) - loss
    T = brentq(lambda x: etot(rn, x) - e1, 1.0e3, 3.0e6)
    r = rn
    while ti < len(targets) and rr_ <= targets[ti]:
        c2 = (pg(r, T) + pr(T))/r
        sp = np.sqrt(c2)/(Om*np.sqrt(cpn))/Rsun
        cosa = abs(vr)/np.hypot(vr, vp)
        cg, ct = cad(r, T, False), cad(r, T, True)
        rows.append((targets[ti], t/om, ph, vr, vp, np.hypot(vr, vp), r, r/RHO_U, T, pr(T)/pg(r, T),
                     np.sqrt(pg(r, T)/r)/1e5, cg/1e5, ct/1e5, v/cg, v/ct, sp, sp/cosa, pnn))
        ti += 1
print('%8s %7s %6s %8s %7s %7s %10s %7s %7s %6s %6s %6s %6s %6s %6s %6s %6s %6s' % ('r', 't_code', 'phi',
      'v_r', 'v_phi', '|v|', 'rho_g/cc', 'code', 'T_K', 'Pr/Pg', 'c_T', 'c_gad', 'c_tot', 'M_g', 'M_tot',
      'sperp', 'arcw', 'Pnn'))
for w in rows:
    print('%8.4f %7.4f %6.2f %8.2f %7.2f %7.2f %10.3e %7.4f %7.0f %6.3f %6.2f %6.2f %6.2f %6.2f %6.2f %6.3f %6.3f %6.2f' % w)
