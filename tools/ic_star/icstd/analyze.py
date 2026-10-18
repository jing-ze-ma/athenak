"""Column analysis and the rules in units of the star (tools/ic_star).

Rules (each reproduces the AG Car E1 production value; constants in spec.DEFAULTS):
  r_in (wall x1min)  'hp_below_cz': n_hp pressure scale heights below the base of the Fe
                     convection zone (the contiguous F_MLT/F > cz_fmlt region containing the
                     Fe opacity peak, or the Fe peak itself if it is radiative);
  column bottom      rlo_fac x1min (AG Car A: 0.95; the code needs ghosts below x1min);
  x1max              rtop_fac R_ph (E1: 3.0); column top col_top_fac x1max (1.15);
  column floor       10^round(log10(rho_floor_rel rho_ph))      (E1: 1e-16; rho_ph 2.56e-11)
  he_sponge_dmax     dmax_rel x floor                            (E1: 2e-16)
  dfloor             dfloor_rel x floor                          (E1: 1e-20)
  he_sponge_rate     10^round(log10(sponge_rate_rel g/c_s,iso)) at R_ph   (E1: 1e-3 /s;
                     g/c_s 1.26e-5 /s for AG Car A: rate = 80 g/c_s)
  he_sponge_r0       sponge_r0_rel x1max                          (E1: 0.9)
  tfloor_kelvin      100 K floor of tfloor_rel T_rad(x1max), T_rad = (E_col/a)^1/4 of the
                     diluted IC radiation at x1max                 (E1: 3000 K; T_rad 3956 K)
  eos_logt_min       floor(10 log10 tfloor)/10                    (E1: 3.4)
  pfloor             10^round(log10(pfloor_rel P_gas(floor, tfloor)))  (E1: 1e-12)
  vceil              10^round(log10(vceil_rel v_esc(R_ph)))       (E1: 1e8; v_esc 2.3e7)
  wb_rmax            wb_rmax_rel R_ph                             (E1: 1.1)
  he_ic_balance_rmax balance_rmax_rel R_ph                        (E1: 342/388.3 = 0.8808)
  mlt_ramp_start / _time  start_turnovers / time_turnovers x t_Fe, t_Fe = H_p/v_MLT (textbook
                     alpha_MLT bubble speed) at the Fe opacity peak     (E1: 2 and 10 x 2.016e6 s)
  implicit_flux_x1min = the column's F_r at x1min, log-linear in r (= he_star_m1.cpp's fine-grid
                     interpolation; the code fatals if they differ by > 1e-3).
"""
import numpy as np

from .micro import A_RAD_, CL


def decade(x):
    return 10.0**np.round(np.log10(x))


def iso_rho_ph(st, c):
    """rho_R of the isothermal atmosphere shot (cheap; used by the column-floor rule).  The
    floor enters the shot only through tau(floor) -> a provisional floor 1e-6 rho_hi..."""
    from .column import _rhs_iso
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    from .micro import SIG, q_hopf
    R, LUM = st.RSTAR, st.LUM
    Tph = (3.0*(LUM/(4*np.pi*R**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
    rhs = _rhs_iso(st, Tph, c['atm_rad'])
    fl = c['rho_lo']*1e-3

    def ev(rr, y):
        return y[0] - np.log(fl)
    ev.terminal = True

    def resid(lnr):
        s = solve_ivp(rhs, [R, 3.0*R], [lnr, 2.0/3.0], events=[ev], method='LSODA',
                      rtol=1e-8, atol=1e-13)
        return -1.0 if s.t_events[0].size == 0 else s.y_events[0][0][1]
    return np.exp(brentq(resid, np.log(c['rho_lo']), np.log(c['rho_hi']), xtol=1e-8)), Tph


def t_rad_at(st, r, Tph):
    R, L = st.RSTAR, st.LUM
    F, Fc = L/(4*np.pi*r**2), L/(4*np.pi*R**2)
    E = (F/CL)*(1.0 + (CL*A_RAD_*Tph**4/Fc - 1.0)*(R/r)**2)
    return (E/A_RAD_)**0.25


def profile(st, col):
    """per-row diagnostics below R_ph (and the full-column tau)"""
    r, rho, T = col['r'], col['rho'], col['T']
    n = r.size
    kap = st.kap(T, rho)
    tau = np.zeros(n)
    for i in range(n - 2, -1, -1):
        tau[i] = tau[i+1] + 0.5*(r[i+1] - r[i])*(rho[i]*kap[i] + rho[i+1]*kap[i+1])
    m = r <= col['rph']*(1 + 1e-12)
    keys = ('nab', 'nad', 'nrad', 'Hp', 'v', 'P', 'pg', 'nab_mlt')
    o = {k: np.full(n, np.nan) for k in keys}
    for i in np.where(m)[0]:
        g = st.mlt_gradient(T[i], rho[i], r[i])
        g0 = st.mlt_gradient(T[i], rho[i], r[i], pp=False)
        o['nab'][i], o['v'][i], o['P'][i], o['Hp'][i] = g[0], g[3], g[5], g[6]
        o['nad'][i], o['nrad'][i], o['pg'][i], o['nab_mlt'][i] = g[7], g[8], g[9], g0[0]
    F = st.LUM/(4*np.pi*r**2)
    gr = st.GM/r**2
    lnP = np.log(np.where(m, o['P'], 1.0))
    with np.errstate(all='ignore'):
        nab_col = np.gradient(np.log(T), r)/np.gradient(lnP, r)
        Fd = 4*A_RAD_*CL*T**3/(3*kap*rho)*np.abs(np.gradient(T, r))
    return dict(kap=kap, tau=tau, below=m, F=F, g=gr, gam_L=kap*F/(CL*gr),
                gam_rad=kap*col['Fr']/(CL*gr), nab_col=nab_col, Fd=Fd,
                Lsum=4*np.pi*r**2*(col['Fr'] + col['fmlt']*F), **o)


def fe_cz(col, pr, cz_fmlt=0.01, lo=None):
    """Fe opacity peak (max kappa_R at 1e5 < T < 1e6 K below R_ph, r >= lo) and the
    contiguous convective zone around it."""
    r, T, fm = col['r'], col['T'], col['fmlt']
    sel = pr['below'] & (T > 1e5) & (T < 1e6)
    if lo is not None:
        sel &= r >= lo
    if not sel.any():
        return None
    ip = np.where(sel)[0][np.argmax(pr['kap'][sel])]
    out = dict(i_fe=ip, r_fe=r[ip], T_fe=T[ip], rho_fe=col['rho'][ip], kap_fe=pr['kap'][ip],
               fe_convective=bool(fm[ip] > cz_fmlt))
    cz = fm > cz_fmlt
    if cz[ip]:
        a = ip
        while a > 0 and cz[a-1]:
            a -= 1
        b = ip
        while b + 1 < r.size and cz[b+1]:
            b += 1
        out.update(i_base=a, r_base=r[a], i_top=b, r_top=r[b])
    else:
        out.update(i_base=ip, r_base=r[ip], i_top=ip, r_top=r[ip])
    # turnover t = H_p / v_MLT at the Fe peak (or at the CZ's max-F_MLT row)
    j = ip
    if not pr['v'][ip] > 0:
        seg = np.arange(out['i_base'], out['i_top'] + 1)
        j = seg[np.argmax(fm[seg])]
    out['i_turn'] = j
    out['t_turn'] = pr['Hp'][j]/pr['v'][j] if pr['v'][j] > 0 else np.nan
    # all convective layers below R_ph (for the report)
    zs, i = [], 0
    mm = pr['below']
    while i < r.size:
        if cz[i] and mm[i]:
            k = i
            while k + 1 < r.size and cz[k+1] and mm[k+1]:
                k += 1
            zs.append((r[i], r[k], float(fm[i:k+1].max())))
            i = k + 1
        else:
            i += 1
    out['zones'] = zs
    return out


def hp_between(r, Hp, r0, r1):
    """number of H_p from r0 up to r1 (trapezoid on the rows)"""
    m = (r >= r0) & (r <= r1) & np.isfinite(Hp)
    if m.sum() < 2:
        return 0.0
    x, h = r[m], Hp[m]
    return float(np.sum(0.5*(1/h[1:] + 1/h[:-1])*np.diff(x)))


def r_at_hp_below(r, Hp, r_base, n_hp):
    """radius n_hp scale heights below r_base (None if the column is not deep enough)"""
    m = (r <= r_base) & np.isfinite(Hp)
    x, h = r[m][::-1], Hp[m][::-1]          # descending from r_base
    acc = 0.0
    for k in range(1, x.size):
        d = 0.5*(1/h[k] + 1/h[k-1])*(x[k-1] - x[k])
        if acc + d >= n_hp:
            return x[k-1] - (n_hp - acc)/(d/(x[k-1] - x[k]))
        acc += d
    return None


def flux_at(col, r0):
    """he_star_m1 fine-grid interpolation of F_r: log-linear in r"""
    r, F = col['r'], col['Fr']
    k = max(0, min(np.searchsorted(r, r0) - 1, r.size - 2))
    w = (r0 - r[k])/(r[k+1] - r[k])
    return float(np.exp((1.0 - w)*np.log(F[k]) + w*np.log(F[k+1])))


def background(st, s, Tph, rho_ph, x1max, col_floor):
    """deck values of the floor / sponge / top layer from the rules"""
    b = s['background']
    R = st.RSTAR
    out, why = {}, {}
    pg = rho_ph*10**st.eos.lp(np.log10(Tph), np.log10(rho_ph))
    cs = np.sqrt(pg/rho_ph)
    g = st.GM/R**2
    nu = g/cs
    out['rho_floor'] = col_floor
    out['he_sponge_dmax'] = b['dmax_rel']*col_floor
    out['dfloor'] = b['dfloor_rel']*col_floor
    out['he_sponge_rate'] = float(b['sponge_rate']) if b.get('sponge_rate') else \
        decade(b['sponge_rate_rel']*nu)
    why['he_sponge_rate'] = 'rule %g x g/c_s,iso at R_ph (= %.3e /s; g %.3e, c_s %.3e)' % (
        b['sponge_rate_rel'], nu, g, cs)
    out['he_sponge_r0'] = b['sponge_r0_rel']*x1max
    trad = t_rad_at(st, x1max, Tph)
    out['T_rad_x1max'] = trad
    out['tfloor_kelvin'] = float(b['tfloor']) if b.get('tfloor') else \
        100.0*np.floor(b['tfloor_rel']*trad/100.0)
    why['tfloor_kelvin'] = 'rule floor100(%g x T_rad(x1max) = %.1f K)' % (b['tfloor_rel'], trad)
    out['eos_logt_min'] = np.floor(10*np.log10(out['tfloor_kelvin']) + 1e-9)/10.0
    pgf = col_floor*10**st.eos.lp(np.log10(out['tfloor_kelvin']), np.log10(col_floor))
    out['pfloor'] = decade(b['pfloor_rel']*pgf)
    why['pfloor'] = 'rule %g x P_gas(floor, tfloor) = %.3e' % (b['pfloor_rel'], pgf)
    vesc = np.sqrt(2*st.GM/R)
    out['vceil'] = decade(b['vceil_rel']*vesc)
    why['vceil'] = 'rule %g x v_esc(R_ph) = %.3e cm/s' % (b['vceil_rel'], vesc)
    out['wb_rmax'] = b['wb_rmax_rel']*R
    out['he_ic_balance_rmax'] = b['balance_rmax_rel']*R
    out['nu_ph'], out['g_ph'], out['cs_ph'], out['v_esc'] = nu, g, cs, vesc
    out['rho_ph'] = rho_ph
    return out, why
