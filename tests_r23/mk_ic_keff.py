#!/usr/bin/env python3
"""tests_r23: build a PURELY RADIATIVE hydrostatic ic_profile for the He4 presupernova
envelope, with an EFFECTIVE (porosity-reduced) opacity kappa_eff = f(r) kappa(T,rho).

WHY.  The present IC (tests_r14/ic_he4_tall_neww.txt) descends from column_sph.py, an
MLT model in which the convective flux carries up to 12 % of L.  The 3-D run has no
subgrid convective flux: it develops real convection, and because the developed envelope
is POROUS the shell-averaged radiative flux exceeds the flux of the shell-mean state
(kappa_eff/kappa = 0.67..0.88 across the convection zone, tests_r19/out_porosity.txt).
Started from the MLT IC the 3-D envelope therefore radiates ~15 % more than L, cools and
contracts.  Literature practice (Jiang et al.) is to start from a purely RADIATIVE
hydrostatic envelope; "accelerated evolution" (Anders et al. 2018) then adjusts the mean
state with the measured fluxes -- which here is exactly what f(r) < 1 does in advance.

WHAT IS BUILT.  The structure integration is column_sph.py's, with ONE change: the MLT
closure is removed and the temperature gradient is the radiative one computed with
kappa_eff,

    nabla = nabla_rad,eff = 3 kappa_eff F p / (4 a c g T^4),   F(r) = L/(4 pi r^2),
    kappa_eff(r, T, rho) = f(r) kappa(T, rho).

Everything else is literally column_sph.py / column.py: hydrostatic dp/dr = -rho GM/r^2
with a point mass, the ideal mu = 1.3423 gas + a T^4/3, the same Rosseland table and its
bilinear interpolation, the same grey Eddington atmosphere T(tau) = (3/4 Teff^4
(tau + 2/3))^1/4 integrated in ln tau from tau_top = 1e-4 down to tau = 2/3 and SHOT so
that the tau = 2/3 surface lands on the anchor radius, the same interior integration
inward in ln p down to r_stop = 0.35 R, the same node counts (NATM = 4000 atmosphere
nodes, NINT = 12000 interior nodes).  Above the convection zone f == 1, so the whole
atmosphere segment and the photosphere are bit-for-bit the original treatment.

The conversion to the run's variables is make_ic_sph.py's, unchanged in substance:
rho is kept, T is re-solved from the EOS TABLE DUMP so that the run's own table returns
the column's TOTAL pressure with UNTAPERED radiation,
    p_gas_table(rho, T) + a T^4/3 = p_column(r),
and then eint = e_gas(rho, T) + w(rho) a T^4 with w the NEW density-only taper of
tests_r14/convert_ic.py (rho_lo 1e-11 -> w = 0, rho_hi 1e-10 -> w = 1, temperature gate
off), i.e. the taper ic_he4_tall_neww.txt is written in.  Above the column's top the
floor atmosphere of tests_r13/mk_tall_ic.py is appended unchanged (rho relaxes to
3e-13 g/cm^3 over 0.005 R, T = the column's top T, rows every 2e7 cm up to 1.32 R).

THE ONE-PARAMETER FAMILY.  With L fixed and the closure fixed, the model is a
one-parameter family labelled by the photospheric radius R_s; a single condition can be
imposed.  --anchor selects it:
    rphot  R_s = R = 2.3717e11 cm, the original builder's own choice (default)
    rhoin  R_s adjusted so that rho(x1min) equals the present IC's rho(x1min)
    mass   R_s adjusted so that the mass in the mesh range equals the present IC's
Matching (rho, T) at x1min simultaneously is NOT possible at fixed L and closure; see
the report printed at the end for how far each model lands from the other two targets.

usage:  mk_ic_keff.py --f0 0.85 --out ic_rad_f085.txt [--anchor rphot|rhoin|mass]
        mk_ic_keff.py --ftab   --out ic_rad_ftab.txt
        mk_ic_keff.py --all                       (builds the four deliverables)
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
PRESN = '/viper/u2/jinma/ATHENAK/bench/hestar_presn'
sys.path.insert(0, PRESN)
import column as col                                              # noqa: E402

DUMP = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
        'eos_table_box_w5.txt')
ICOLD = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r14/ic_he4_tall_neww.txt'

# star (stars.json, he4_presn) and the run's mesh
LOGL, MNOM, LOGG = 4.78, 3.15, 3.8712
XC, YC, ZC = 0.0, 0.98, 0.02
TAB = '/viper/u2/jinma/ATHENAK/bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt'
RSTAR = 2.3717e11
X1MIN, X1MAX = 1.18585e11, 2.964625e11         # 0.50 R .. 1.25 R (mesh/x1max, tests_r15)

# convection zone and the porosity window
CZ_LO, CZ_HI, CZ_BL = 0.635, 0.967, 0.02       # in units of R; cubic blend width
# tests_r19/out_porosity.txt + tests_r21/out_chan.txt, kappa_eff/kappa at 3.5 turnovers
FTAB_R = np.array([0.60, 0.72, 0.80, 0.88, 0.93, 0.96])
FTAB_F = np.array([0.73, 0.88, 0.88, 0.80, 0.67, 0.81])

TAU_TOP, R_STOP, NATM, NINT = 1.0e-4, 0.35, 4000, 12000
R_BELOW, NDOWN = 0.45, 1500          # anchor = inner: how far below x1min to extend
RHO_ATM, RTOP_ATM, DR_ATM = 3.0e-13, 1.32*RSTAR, 2.0e7   # mk_tall_ic.py

# new EOS radiation taper (tests_r14/convert_ic.py)
RHO_LO_NEW, RHO_HI_NEW = 1e-11, 1e-10
A_RAD = 7.5657332503e-15


# ----------------------------------------------------------------- EOS table dump
def read_dump(fn):
    with open(fn) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    return (xmin + dx*np.arange(nx), ymin + dy*np.arange(ny),
            d[:, 0].reshape(ny, nx), d[:, 1].reshape(ny, nx))


def bilin(x, y, f, xq, yq):
    ix = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    iy = np.clip(np.searchsorted(y, yq) - 1, 0, len(y) - 2)
    tx = (xq - x[ix])/(x[ix+1] - x[ix])
    ty = (yq - y[iy])/(y[iy+1] - y[iy])
    return ((1-tx)*(1-ty)*f[iy, ix] + tx*(1-ty)*f[iy, ix+1]
            + (1-tx)*ty*f[iy+1, ix] + tx*ty*f[iy+1, ix+1])


def smooth(s):
    s = np.clip(s, 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def w_new(rho):
    s = (np.log10(rho) - np.log10(RHO_LO_NEW))/(np.log10(RHO_HI_NEW)
                                                - np.log10(RHO_LO_NEW))
    return smooth(s)


# ----------------------------------------------------------------- porosity factor f(r)
def make_f(f0=None, tab=False):
    """f(r): 1 outside the convection zone, a plateau f0 (or the measured table) inside,
    cubic-blended over CZ_BL*R at each edge."""
    lo, hi, bl = CZ_LO*RSTAR, CZ_HI*RSTAR, CZ_BL*RSTAR

    def win(r):
        return (smooth((r - (lo - 0.5*bl))/bl)*smooth(((hi + 0.5*bl) - r)/bl))

    if tab:
        rt = FTAB_R*RSTAR

        def f(r):
            fi = np.interp(r, rt, FTAB_F)      # clamped at the ends
            return 1.0 + (fi - 1.0)*win(r)
    else:
        def f(r):
            return 1.0 + (f0 - 1.0)*win(r)
    return f


# ----------------------------------------------------------------- structure
def rho_of(p, T):
    """column.rho_of with a floor: near/above Eddington the radiative model drives
    p_gas -> 0 and we want to see where, not crash."""
    pg = max(p - col.a_rad*T**4/3.0, 1e-10*p)
    return pg*col.MU*col.mu_amu/(col.kB*T)


def build(fofr, r_anchor, tau_top=TAU_TOP):
    """Integrate the radiative-equilibrium envelope whose tau = 2/3 surface sits at
    r_anchor.  Returns arrays sorted by ascending r plus scalars."""
    L = 10.0**LOGL*col.Lsun
    M = MNOM*col.Msun
    GM = col.G*M
    Teff = (L/(4.0*np.pi*r_anchor**2*col.sigma_sb))**0.25

    def gof(r):
        return GM/(r*r)

    def Fof(r):
        return L/(4.0*np.pi*r*r)

    def Tof(tau):
        return (0.75*Teff**4*(tau + 2.0/3.0))**0.25

    # ---- grey Eddington atmosphere, identical to column_sph.run_sph -----------------
    def seed(r0):
        g0, F0 = gof(r0), Fof(r0)
        T0 = Tof(tau_top)
        pg = 1.0
        for _ in range(500):
            rho = pg*col.MU*col.mu_amu/(col.kB*T0)
            k = col.kappa(T0, rho)
            pn = max((g0/k - F0/col.c_l), 1.0e-6*g0/k)*tau_top
            if abs(pn - pg)/pg < 1.0e-13:
                pg = pn
                break
            pg = 0.5*(pg + pn)
        return pg + col.a_rad*T0**4/3.0

    def atm(lt, y):
        tau = np.exp(lt)
        T = Tof(tau)
        p, r = y
        rho = rho_of(p, T)
        k = col.kappa(T, rho)
        return [tau*gof(r)/k, -tau/(k*rho)]

    def shoot(r0):
        return solve_ivp(atm, [np.log(tau_top), np.log(2.0/3.0)], [seed(r0), r0],
                         rtol=1.0e-11, atol=[1.0e-10, 1.0e-3], dense_output=True,
                         max_step=0.01, first_step=1.0e-4)

    r0 = r_anchor*1.02
    s = shoot(r0)
    for _ in range(40):
        f = s.y[1, -1] - r_anchor
        if abs(f)/r_anchor < 1.0e-13:
            break
        r0b = r0 - f
        sb = shoot(r0b)
        fb = sb.y[1, -1] - r_anchor
        r0n = r0b - fb*(r0b - r0)/(fb - f) if fb != f else r0b
        r0, s = r0n, shoot(r0n)
    p_ph, r_ph = s.y[0, -1], s.y[1, -1]

    R_, P_, T_, TAU_, F_ = [], [], [], [], []
    for lt in np.linspace(np.log(tau_top), np.log(2.0/3.0), NATM):
        tau = np.exp(lt)
        pp, rr = s.sol(lt)
        R_.append(rr)
        P_.append(pp)
        T_.append(Tof(tau))
        TAU_.append(tau)
        F_.append(float(fofr(np.array([rr]))[0]))

    # ---- interior, inward in ln p, RADIATIVE with kappa_eff ------------------------
    def deriv(lnp, y):
        p = np.exp(lnp)
        T, tau, r = y
        g, F = gof(r), Fof(r)
        rho = rho_of(p, T)
        k = col.kappa(T, rho)
        keff = float(fofr(np.array([r]))[0])*k
        nrad = 3.0*keff*F*p/(4.0*col.a_rad*col.c_l*g*T**4)
        return [T*nrad, p*k/g, -p/(rho*g)]

    r_stop = R_STOP*RSTAR

    def stopr(lnp, y):
        return y[2] - r_stop
    stopr.terminal = True
    stopr.direction = -1

    si = solve_ivp(deriv, [np.log(p_ph), np.log(p_ph) + 60.0], [Teff, 2.0/3.0, r_ph],
                   rtol=1.0e-10, atol=[1.0e-4, 1.0e-10, 1.0e-2], max_step=0.003,
                   events=stopr, dense_output=True, first_step=1.0e-4)
    for lnp in np.linspace(np.log(p_ph), si.t[-1], NINT)[1:]:
        T, tau, r = si.sol(lnp)
        R_.append(r)
        P_.append(np.exp(lnp))
        T_.append(T)
        TAU_.append(tau)
        F_.append(float(fofr(np.array([r]))[0]))

    r = np.array(R_)
    o = np.argsort(r)
    r, p, T, tau, fr = (r[o], np.array(P_)[o], np.array(T_)[o], np.array(TAU_)[o],
                        np.array(F_)[o])
    rho = np.array([rho_of(pp, tt) for pp, tt in zip(p, T)])
    kap = np.array([col.kappa(tt, dd) for tt, dd in zip(T, rho)])
    reached = si.t_events[0].size > 0
    return dict(r=r, p=p, T=T, tau=tau, rho=rho, kappa=kap, f=fr, Teff=Teff, L=L, M=M,
                GM=GM, r_ph=r_ph, p_ph=p_ph, reached=reached, r_end=r[0])


def T_ideal_of(rho, p):
    """invert p = rho k T/(mu m_u) + a T^4/3 for the COLUMN's ideal-gas + radiation T."""
    lo, hi = 1.0, 1.0e8
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if rho*col.kB*mid/(col.MU*col.mu_amu) + col.a_rad*mid**4/3.0 < p:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def build_up(fofr, r_in, rho_in, p_in):
    """THE DEFAULT PATH.  March the radiative-equilibrium envelope OUTWARD in ln p from
    the inner boundary state (r_in, rho_in, p_in), stop at the photosphere -- defined
    exactly as the grey Eddington tau = 2/3 level, sigma T^4 = F(r) -- and attach the
    original builder's grey atmosphere above it.  Also extends a short way BELOW r_in
    (to R_BELOW R) so the pgen has profile under the ghost zones of x1min.

    Direction matters.  In this envelope Prad/Pgas ~ 50 and d p_gas/dp = 1 - Gamma_eff,
    so p_gas is a tiny residual of two large terms and the INWARD march amplifies any
    change in the closure exponentially (the surface-anchored build() below gives 0.63 /
    3.06 / 687 times the present envelope mass for f0 = 1 / 0.85 / 0.75).  The OUTWARD
    march is the contracting direction: p_gas relaxes onto its near-Eddington attractor
    and the model forgets its start within a few beta H_p.  Anchoring at r_in also keeps
    the deep envelope -- where the thermal time is longest and the run cannot adjust --
    identical to the present IC by construction.
    """
    L = 10.0**LOGL*col.Lsun
    M = MNOM*col.Msun
    GM = col.G*M
    T_in = T_ideal_of(rho_in, p_in)

    def gof(r):
        return GM/(r*r)

    def Fof(r):
        return L/(4.0*np.pi*r*r)

    def deriv(lnp, y):
        p = np.exp(lnp)
        T, r = y
        g, F = gof(r), Fof(r)
        rho = rho_of(p, T)
        k = col.kappa(T, rho)
        keff = float(fofr(np.array([r]))[0])*k
        nrad = 3.0*keff*F*p/(4.0*col.a_rad*col.c_l*g*T**4)
        return [T*nrad, -p/(rho*g)]

    def phot(lnp, y):
        T, r = y
        return T - (Fof(r)/col.sigma_sb)**0.25
    phot.terminal = True
    phot.direction = -1

    lp0 = np.log(p_in)
    su = solve_ivp(deriv, [lp0, lp0 - 60.0], [T_in, r_in], rtol=1.0e-10,
                   atol=[1.0e-4, 1.0e-2], max_step=0.003, events=phot,
                   dense_output=True, first_step=1.0e-4)
    if su.t_events[0].size == 0:
        raise RuntimeError('outward march never reached the photosphere (ended at '
                           'r = %.4e, T = %.4e)' % (su.y[1, -1], su.y[0, -1]))
    lp_ph = float(su.t_events[0][0])
    Teff, r_ph = float(su.y_events[0][0][0]), float(su.y_events[0][0][1])
    p_ph = np.exp(lp_ph)

    # ---- below the anchor: the same ODE inward, a short way only -------------------
    def stopr(lnp, y):
        return y[1] - R_BELOW*RSTAR
    stopr.terminal = True
    stopr.direction = -1
    sd = solve_ivp(deriv, [lp0, lp0 + 20.0], [T_in, r_in], rtol=1.0e-10,
                   atol=[1.0e-4, 1.0e-2], max_step=0.003, events=stopr,
                   dense_output=True, first_step=1.0e-4)
    lp_bot = float(sd.t_events[0][0]) if sd.t_events[0].size else float(sd.t[-1])

    # ---- the grey Eddington atmosphere, outward from tau = 2/3 ---------------------
    def Tof(tau):
        return (0.75*Teff**4*(tau + 2.0/3.0))**0.25

    def atm(lt, y):
        tau = np.exp(lt)
        T = Tof(tau)
        p, r = y
        rho = rho_of(p, T)
        k = col.kappa(T, rho)
        return [tau*gof(r)/k, -tau/(k*rho)]
    sa = solve_ivp(atm, [np.log(2.0/3.0), np.log(TAU_TOP)], [p_ph, r_ph], rtol=1.0e-11,
                   atol=[1.0e-10, 1.0e-3], dense_output=True, max_step=0.01,
                   first_step=1.0e-4)

    R_, P_, T_ = [], [], []
    ndn = int(NDOWN*max(lp_bot - lp0, 0.0)/max(lp0 - lp_ph, 1e-9)) + 2
    for lnp in np.linspace(lp_bot, lp0, max(ndn, 2))[:-1]:
        T, r = sd.sol(lnp)
        R_.append(r)
        P_.append(np.exp(lnp))
        T_.append(T)
    for lnp in np.linspace(lp0, lp_ph, NINT):
        T, r = su.sol(lnp)
        R_.append(r)
        P_.append(np.exp(lnp))
        T_.append(T)
    for lt in np.linspace(np.log(2.0/3.0), np.log(TAU_TOP), NATM)[1:]:
        pp, rr = sa.sol(lt)
        R_.append(rr)
        P_.append(pp)
        T_.append(Tof(np.exp(lt)))

    r = np.array(R_)
    o = np.argsort(r)
    r, p, T = r[o], np.array(P_)[o], np.array(T_)[o]
    rho = np.array([rho_of(pp, tt) for pp, tt in zip(p, T)])
    kap = np.array([col.kappa(tt, dd) for tt, dd in zip(T, rho)])
    fr = fofr(r)
    # tau, integrated inward from tau_top at the outermost node
    tau = np.zeros_like(r)
    tau[-1] = TAU_TOP
    krho = kap*rho
    for i in range(len(r) - 2, -1, -1):
        tau[i] = tau[i+1] + 0.5*(krho[i] + krho[i+1])*(r[i+1] - r[i])
    return dict(r=r, p=p, T=T, tau=tau, rho=rho, kappa=kap, f=fr, Teff=Teff, L=L, M=M,
                GM=GM, r_ph=r_ph, p_ph=p_ph, reached=True, r_end=r[0], T_in=T_in)


def mass_in_mesh(r, rho):
    m = (r >= X1MIN) & (r <= X1MAX)
    return float(np.trapezoid(4.0*np.pi*r[m]**2*rho[m], r[m]))


def solve_anchor(fofr, mode, old):
    """secant on R_s for the requested anchor condition."""
    if mode == 'rphot':
        return RSTAR, build(fofr, RSTAR)
    if mode == 'inner':
        rho_in = float(np.exp(np.interp(X1MIN, old['r'], np.log(old['rho']))))
        p_in = float(np.exp(np.interp(X1MIN, old['r'], np.log(old['p']))))
        st = build_up(fofr, X1MIN, rho_in, p_in)
        print('  anchor inner: r_in = %.6e (%.4f R), rho_in = %.6e, p_in = %.6e, '
              'T_in(ideal) = %.6e ; emergent R_phot = %.6e (%.5f R), Teff = %.1f K'
              % (X1MIN, X1MIN/RSTAR, rho_in, p_in, st['T_in'], st['r_ph'],
                 st['r_ph']/RSTAR, st['Teff']))
        return st['r_ph'], st
    ro, rhoo = old['r'], old['rho']
    tgt_rho = float(np.exp(np.interp(X1MIN, ro, np.log(rhoo))))
    tgt_m = mass_in_mesh(ro, rhoo)

    def resid(rs):
        st = build(fofr, rs)
        if mode == 'rhoin':
            v = float(np.exp(np.interp(X1MIN, st['r'], np.log(st['rho']))))
            return np.log(v/tgt_rho), st
        v = mass_in_mesh(st['r'], st['rho'])
        return np.log(v/tgt_m), st

    a, b = RSTAR, RSTAR*1.01
    fa, sa = resid(a)
    fb, sb = resid(b)
    for _ in range(30):
        if abs(fb) < 1e-6 or fb == fa:
            break
        c = b - fb*(b - a)/(fb - fa)
        c = float(np.clip(c, 0.8*RSTAR, 1.3*RSTAR))
        a, fa = b, fb
        fb, sb = resid(c)
        b = c
    print('  anchor %s: R_s = %.6e (%.5f R), residual ln = %.2e' % (mode, b, b/RSTAR,
                                                                    fb))
    return b, sb


# ----------------------------------------------------------------- write the ic_profile
def to_ic(st, x, y, le, lp):
    """make_ic_sph.py: keep rho, solve T_table from p_table(rho,T) = p_column, then
    eint = e_gas + w_new(rho) a T^4.  Then append mk_tall_ic.py's floor atmosphere."""
    r, p, rho = st['r'], st['p'], st['rho']
    lr = np.log10(rho)

    def ptab(lt):
        return rho*10.0**bilin(x, y, lp, lr, lt) + A_RAD*10.0**(4.0*lt)/3.0
    lo = np.full_like(lr, y[0] + 1e-9)
    hi = np.full_like(lr, y[-1] - 1e-9)
    for _ in range(200):
        mid = 0.5*(lo + hi)
        f = ptab(mid) - p
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    lt = 0.5*(lo + hi)
    Tt = 10.0**lt
    resid = np.abs(ptab(lt)/p - 1.0)
    inside = ((lr > x[0]) & (lr < x[-1]) & (lt > y[0] + 1e-6) & (lt < y[-1] - 1e-6))
    eg = rho*10.0**bilin(x, y, le, lr, lt)
    eint = eg + w_new(rho)*A_RAD*Tt**4

    # floor atmosphere (mk_tall_ic.py)
    Tend = Tt[-1]
    rn = np.arange(r[-1] + DR_ATM, RTOP_ATM, DR_ATM)
    xx = (rn - r[-1])/(0.005*RSTAR)
    rhon = RHO_ATM + (rho[-1] - RHO_ATM)*np.exp(-xx)
    Tn = np.full_like(rn, Tend)
    egn = rhon*10.0**bilin(x, y, le, np.log10(rhon), np.log10(Tn))
    en = egn + w_new(rhon)*A_RAD*Tn**4
    return (np.concatenate((r, rn)), np.concatenate((rho, rhon)),
            np.concatenate((eint, en)), Tt, np.concatenate((Tt, Tn)),
            resid.max(), bool(inside.all()), lr.min(), x[0])


# ----------------------------------------------------------------- report
def report(name, st, r_all, rho_all, e_all, T_all, old, extra):
    R = RSTAR
    r, rho, T, kap, fr, tau = st['r'], st['rho'], st['T'], st['kappa'], st['f'], st['tau']
    L, M = st['L'], st['M']
    edd = 4.0*np.pi*col.G*M*col.c_l
    gam = fr*kap*L/edd
    rph = float(np.interp(np.log(2.0/3.0), np.log(tau[::-1]), r[::-1]))
    mm = mass_in_mesh(r_all, rho_all)
    mold = mass_in_mesh(old['r'], old['rho'])
    ig = np.argmax(gam)
    ngam = (gam > 1.0).sum()
    lines = []
    lines.append('=== %s' % name)
    lines.append('  R_s = %.6e (%.5f R) ; R_phot(tau=2/3) = %.6e (%.5f R) ; Teff = %.0f K'
                 % (extra['rs'], extra['rs']/R, rph, rph/R, st['Teff']))
    lines.append('  column top r = %.6e (%.5f R), tau_top = %.1e ; base r = %.6e '
                 '(%.4f R) ; ODE reached the base: %s'
                 % (r[-1], r[-1]/R, tau[-1], r[0], r[0]/R, st['reached']))
    lines.append('  mass in [%.4e, %.4e] = [%.3f, %.3f] R : %.5e g  (old IC %.5e g, '
                 'ratio %.4f)' % (X1MIN, X1MAX, X1MIN/R, X1MAX/R, mm, mold, mm/mold))
    lines.append('  max Gamma_rad = kappa_eff L/(4 pi G M c) = %.4f at r/R %.4f '
                 '(T %.3e, rho %.3e, kappa %.3f, f %.3f)'
                 % (gam[ig], r[ig]/R, T[ig], rho[ig], kap[ig], fr[ig]))
    if ngam:
        m = gam > 1.0
        lines.append('  Gamma > 1 on %d of %d nodes, r/R = %.4f .. %.4f'
                     % (ngam, len(r), r[m].min()/R, r[m].max()/R))
    else:
        lines.append('  Gamma > 1 nowhere')
    # density / gas-pressure inversion
    pg = st['p'] - col.a_rad*T**4/3.0
    d = np.diff(rho)
    if np.any(d > 0):
        b = np.where(d > 0)[0]
        # contiguous run containing the strongest inversion
        lines.append('  DENSITY INVERSION over r/R %.4f .. %.4f : rho %.4e -> %.4e, '
                     'rho_max/rho_min = %.3f'
                     % (r[b.min()]/R, r[b.max()+1]/R, rho[b.min()], rho[b.max()+1],
                        rho[b.min():b.max()+2].max()/rho[b.min():b.max()+2].min()))
        lines.append('  p_gas/p at the inversion: min %.3e (at r/R %.4f) ; p_gas '
                     'floored: %s'
                     % ((pg/st['p']).min(),
                        r[int(np.argmin(pg/st['p']))]/R,
                        bool((pg <= 1.01e-10*st['p']).any())))
    else:
        lines.append('  no density inversion (rho monotone outward-decreasing)')
    lines.append('  profile vs the present IC (T from the EOS-table solve):')
    lines.append('     r/R      rho_new      rho_old  ratio       T_new       T_old'
                 '  ratio')
    for s in (0.6, 0.7, 0.8, 0.9, 0.95, 1.0):
        rr = s*R
        a1 = float(np.exp(np.interp(rr, r_all, np.log(rho_all))))
        b1 = float(np.exp(np.interp(rr, old['r'], np.log(old['rho']))))
        a2 = float(np.exp(np.interp(rr, r_all, np.log(T_all))))
        b2 = float(np.exp(np.interp(rr, old['r'], np.log(old['T']))))
        lines.append('   %5.2f  %11.4e  %11.4e  %5.3f  %10.4e  %10.4e  %5.3f'
                     % (s, a1, b1, a1/b1, a2, b2, a2/b2))
    lines.append('  fits inside x1max = %.6e : %s (file top %.6e = %.4f R)'
                 % (X1MAX, bool(r_all.max() >= X1MAX), r_all.max(), r_all.max()/R))
    lines.append('  EOS-table solve: max |p_bilin/p_col - 1| = %.2e ; every node inside '
                 'the table: %s (min log10 rho %.3f vs table edge %.3f)'
                 % (extra['resid'], extra['inside'], extra['lrmin'], extra['ledge']))
    mlow = rho_all < 10.0**extra['ledge']
    lines.append('  rows below the EOS table density edge: %d of %d%s'
                 % (mlow.sum(), len(rho_all),
                    ' (all in the floor atmosphere)' if mlow.sum() and
                    r_all[mlow].min() > r[-1] else ''))
    return '\n'.join(lines)


HDR = """# tests_r23/mk_ic_keff.py : PURELY RADIATIVE hydrostatic ic_profile for the
# 4.0 Msun He star presupernova envelope, with an EFFECTIVE opacity
#     kappa_eff(r) = f(r) kappa(T,rho),   f = %s
# f = 1 outside the convection zone (%.3f .. %.3f R), cubic-blended over %.3f R.
# NO MLT / no subgrid convective flux: radiation carries all of L at every radius.
# Structure: column_sph.py's integration (point-mass hydrostatics, ideal mu = %.4f gas
# + a T^4/3, Rosseland table %s, grey Eddington atmosphere shot so tau = 2/3 lands on
# R_s = %.6e cm, inward in ln p to 0.35 R) with nabla = nabla_rad,eff instead of the MLT
# nabla.  Conversion: make_ic_sph.py (rho kept; T solved from the EOS table dump so that
# p_gas_table(rho,T) + a T^4/3 = p_column) and the NEW density-only radiation taper of
# tests_r14/convert_ic.py, eint = e_gas(rho,T) + w(rho) a T^4, w = 0 below rho %.1e,
# 1 above %.1e, temperature gate off.  Floor atmosphere above the column top from
# tests_r13/mk_tall_ic.py (rho -> %.1e over 0.005 R, T = the top T, to %.3f R).
# anchor = %s ; R_phot(tau=2/3) = %.6e cm ; Teff = %.1f K
# r[cm]  rho[g/cm^3]  eint[erg/cm^3]
"""


def one(tag, fofr, fdesc, anchor, old, tbl, out):
    x, y, le, lp = tbl
    print('--- building %s (f = %s, anchor = %s)' % (out, fdesc, anchor))
    rs, st = solve_anchor(fofr, anchor, old)
    (r_all, rho_all, e_all, Tt, T_all, resid, inside,
     lrmin, ledge) = to_ic(st, x, y, le, lp)
    fn = os.path.join(HERE, out)
    with open(fn, 'w') as fh:
        fh.write(HDR % (fdesc, CZ_LO, CZ_HI, CZ_BL, col.MU, os.path.basename(TAB), rs,
                        RHO_LO_NEW, RHO_HI_NEW, RHO_ATM, RTOP_ATM/RSTAR, anchor,
                        float(np.interp(np.log(2.0/3.0), np.log(st['tau'][::-1]),
                                        st['r'][::-1])), st['Teff']))
        for a, b, c in zip(r_all, rho_all, e_all):
            fh.write('%.10e %.10e %.10e\n' % (a, b, c))
    txt = report(out, st, r_all, rho_all, e_all, T_all, old,
                 dict(rs=rs, resid=resid, inside=inside, lrmin=lrmin, ledge=ledge))
    print(txt)
    print('  wrote %s (%d rows)\n' % (fn, len(r_all)))
    return txt


def load_old(tbl):
    x, y, le, lp = tbl
    d = np.loadtxt(ICOLD, comments='#')
    r, rho, eint = d[:, 0], d[:, 1], d[:, 2]
    # T of the present IC under the NEW taper (bisection, poros.py/relax_ic.py)
    lrho = np.log10(rho)
    w = w_new(rho)
    lo = np.full_like(r, y[0] + 1e-9)
    hi = np.full_like(r, y[-1] - 1e-9)
    for _ in range(80):
        mid = 0.5*(lo + hi)
        f = rho*10.0**bilin(x, y, le, lrho, mid) + w*A_RAD*10.0**(4*mid) - eint
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    T = 10.0**(0.5*(lo + hi))
    p = rho*10.0**bilin(x, y, lp, lrho, np.log10(T)) + A_RAD*T**4/3.0
    return dict(r=r, rho=rho, eint=eint, T=T, p=p)


def main():
    col.set_composition(XC, YC, ZC, tab=TAB)
    tbl = read_dump(DUMP)
    old = load_old(tbl)
    av = sys.argv[1:]
    anchor = av[av.index('--anchor') + 1] if '--anchor' in av else 'inner'
    if '--all' in av:
        jobs = [('f100', make_f(1.0), '1.00 (pure radiative, real opacity)',
                 'ic_rad_f100.txt'),
                ('f085', make_f(0.85), '0.85 plateau', 'ic_rad_f085.txt'),
                ('f075', make_f(0.75), '0.75 plateau', 'ic_rad_f075.txt'),
                ('ftab', make_f(tab=True), 'measured table 0.73/0.88/0.88/0.80/0.67/0.81'
                 ' at r/R 0.60/0.72/0.80/0.88/0.93/0.96', 'ic_rad_ftab.txt')]
    elif '--fmatch' in av:
        # solve for the plateau f0 whose envelope mass in the mesh equals the present
        # IC's, at anchor = rphot.  This is the best-conditioned member of the family:
        # the run conserves mass and cannot change the deep envelope on its own time.
        tgt = mass_in_mesh(old['r'], old['rho'])

        def res(f0):
            st = build(make_f(f0), RSTAR)
            v = mass_in_mesh(st['r'], st['rho'])
            print('    f0 = %.5f -> mass %.5e (ratio %.4f)' % (f0, v, v/tgt))
            return np.log(v/tgt)
        a, b = 1.00, 0.90
        fa, fb = res(a), res(b)
        for _ in range(20):
            if abs(fb) < 3e-3 or fb == fa:
                break
            c = float(np.clip(b - fb*(b - a)/(fb - fa), 0.70, 1.05))
            a, fa, b = b, fb, c
            fb = res(c)
        print('  mass-matched plateau f0 = %.5f' % b)
        jobs = [('fmatch', make_f(b), '%.4f plateau (solved so that the envelope mass '
                 'in the mesh equals the present IC\'s)' % b, 'ic_rad_fmass.txt')]
    elif '--ftab' in av:
        jobs = [('ftab', make_f(tab=True), 'measured table',
                 av[av.index('--out') + 1])]
    else:
        f0 = float(av[av.index('--f0') + 1])
        jobs = [('f%03d' % round(100*f0), make_f(f0), '%.2f plateau' % f0,
                 av[av.index('--out') + 1])]
    out = []
    for tag, fofr, fdesc, fn in jobs:
        try:
            out.append(one(tag, fofr, fdesc, anchor, old, tbl, fn))
        except Exception as e:                                       # noqa: BLE001
            msg = '=== %s FAILED: %s: %s' % (fn, type(e).__name__, e)
            print(msg)
            out.append(msg)
    print('\n\n############ SUMMARY (anchor = %s)\n' % anchor)
    print('\n'.join(out))


if __name__ == '__main__':
    main()
