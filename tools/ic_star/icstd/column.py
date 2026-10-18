"""The standard IC column (the AG Car method), tools/ic_star.

ONE set of equations (make_ic_agcar_ge.py, Viper lbv_1008/agcar/geos/scripts, = AG Car A/B and
the E1 production IC ic_agcar_A_ge_deep70.txt):
  * R_ph (tau_R = 2/3) given, or sqrt(L/(4 pi sigma Teff^4)); T(tau) grey Hopf (Mihalas fit);
  * ABOVE R_ph: isothermal atmosphere at T_ph = T_Hopf(2/3), gravity only, rho_R SHOT so that
    its tau(R_ph) = 2/3, down to the column floor rho_floor, constant floor above to r_col_top;
    E above R_ph DILUTED: E = (F/c)(1 + (c a T_ph^4/F_ph - 1)(R_ph/r)^2) (he_ic_rad = sph1d
    replaces E and re-marches the thin layers in the code: this only seeds it);
  * R_ph inward to tau = tau_j (3): Hopf T(tau) + gas HSE with the radiation force of F;
  * tau_j inward to the column bottom: dP_tot/dr = -rho GM/r^2, dlnT/dlnP_tot = nabla (local MLT,
    alpha_MLT, + MLT++ f1, f2, pp_alpha; = nabla_rad where radiative), F_r = F nabla/nabla_rad,
    F_MLT/F = 1 - nabla/nabla_rad, F = L/(4 pi r^2) (L constant: no sources in the envelope);
  * GM a point mass (he_gm).
Output rows: r rho eint F_r(radiative) E T F_MLT/F (he_ic_cols = 5, 7 columns).

Sampling schemes (same equations; they differ only in where the rows sit and the ODE
tolerances):
  'agcar_ge' (THE STANDARD; = make_ic_agcar_ge.py): interior rows exp(arange(ln r_bot, ln r_j,
       dlnr)) plus the wall radius r_wall as an exact row; Hopf rows arange(r_j, R_ph, drn);
       atmosphere rows arange(R_ph, r_floor, drn), floor rows every fl_step*drn to r_col_top
       (+ one step); one interior solve r_j -> r_bot (LSODA rtol 1e-12, atol 1e-12, dense).
  legacy_deep (only for reproducing ic_agcar_A_ge_deep70.txt = make_ic_agcar_deep.py on top of
       ic_agcar_A_ge.txt): stage 1 = the standard to r_bot1 = rlo_fac r_wall1; stage 2 restarts
       from stage 1's first row AS PRINTED (%.10e) inward on the continued log grid with rtol
       1e-10, atol 1e-13 to r_bot.
  'mltpp_inward' (compatibility: make_ic_agcar_mltpp.py --mode inward, = hg145 ic2/gen.sh): log
       grid anchored at r_j, Hopf rows r_j + k drn and the R_ph row, iso atmosphere shot with
       rtol 1e-10, interior rtol 1e-9 atol 1e-12, MLT++ in the 'f2' form.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .micro import CL, SIG, A_RAD_, q_hopf


def fmt_row(row):
    return '%.10e %.10e %.10e %.10e %.10e %.10e %.6e\n' % tuple(row)


# --------------------------------------------------------------------------- shared pieces
def _rhs_hopf(st, GM=None):
    LUM, eos = st.LUM, st.eos
    GM = st.GM if GM is None else GM

    def rhs1(rr, y):                    # Hopf layer, y = (ln rho, tau)
        lnr, tau = y
        tau = max(tau, 1.0e-8)
        rho_ = np.exp(lnr)
        F = LUM/(4*np.pi*rr**2)
        q = q_hopf(tau)
        T_ = (3.0*F/(4.0*SIG)*(tau + q))**0.25
        kap_ = st.kap(T_, rho_)
        dtau = -kap_*rho_
        dlnT = 0.25*(-2.0/rr + (1.0 + 0.1331*3.4488*np.exp(-3.4488*tau))*dtau/(tau + q))
        lt_, lr_ = np.log10(T_), np.log10(rho_)
        p_ = rho_*10.0**eos.lp(lt_, lr_)
        geff = GM/rr**2 - kap_*F/CL
        dlnp = -rho_*geff/p_
        dlnr = (dlnp - eos.lp(lt_, lr_, dx=1)*dlnT)/(1.0 + eos.lp(lt_, lr_, dy=1))
        return [dlnr, dtau]
    return rhs1


def _rhs_iso(st, Tph, rad=0, GM=None):
    LUM, eos = st.LUM, st.eos
    GM = st.GM if GM is None else GM

    def rhs_iso(rr, y):                 # isothermal at T_ph, y = (ln rho, tau)
        rho_ = np.exp(y[0])
        kap_ = st.kap(Tph, rho_)
        lt_, lr_ = np.log10(Tph), np.log10(rho_)
        p_ = rho_*10.0**eos.lp(lt_, lr_)
        if rad:
            F = LUM/(4*np.pi*rr**2)
            geff = GM/rr**2 - kap_*F/CL
        else:
            geff = GM/rr**2
        return [-rho_*geff/p_/(1.0 + eos.lp(lt_, lr_, dy=1)), -kap_*rho_]
    return rhs_iso


def _rhs_interior(st, with_tau=False):
    GM = st.GM

    def rhs2(rr, y):                    # MLT++ interior, y = (ln T, ln rho[, tau])
        T_, rho_ = np.exp(y[0]), np.exp(y[1])
        o = st.mlt_gradient(T_, rho_, rr)
        nab = o[0]
        P_, pg_, chir, chit, cp, dlt, nad_ = st.thermo(T_, rho_)
        dlnP = -rho_*(GM/rr**2)/P_
        dlnT = nab*dlnP
        if with_tau:
            return [dlnT, (dlnP - chit*dlnT)/chir, -o[4]*rho_]
        return [dlnT, (dlnP - chit*dlnT)/chir]
    return rhs2


def _mlt_cols(st, r, T, rho, isc):
    n = r.size
    ffr, fm = np.ones(n), np.zeros(n)
    for i in np.where(isc)[0]:
        o = st.mlt_gradient(T[i], rho[i], r[i])
        ffr[i], fm[i] = o[1], o[2]
    return ffr, fm


def _ecol(st, r, T, rph, tph):
    F = st.LUM/(4*np.pi*r**2)
    E = A_RAD_*T**4
    up = r > rph
    Ec, Fc = A_RAD_*tph**4, st.LUM/(4*np.pi*rph**2)
    return np.where(up, (F/CL)*(1.0 + (CL*Ec/Fc - 1.0)*(rph/r)**2), E)


# --------------------------------------------------------------------------- the standard
def column_agcar_ge(st, c):
    """make_ic_agcar_ge.py column() + sample(); c: dict of controls (see spec defaults).
    Returns dict r, rho, T, eint, Fr, E, fmlt + landmarks."""
    LUM, R = st.LUM, st.RSTAR
    eos = st.eos
    rtol = c['rtol']
    rhs1 = _rhs_hopf(st)
    Tph = (3.0*(LUM/(4*np.pi*R**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
    rhs_iso = _rhs_iso(st, Tph, c['atm_rad'])
    rho_floor = c['rho_floor']

    def top_event(rr, y):
        return y[0] - np.log(rho_floor)
    top_event.terminal = True

    def shoot(lnr):
        return solve_ivp(rhs_iso, [R, 3.0*R], [lnr, 2.0/3.0], events=[top_event],
                         method='LSODA', rtol=rtol, atol=1e-13, dense_output=True)

    def resid(lnr):
        s = shoot(lnr)
        return -1.0 if s.t_events[0].size == 0 else s.y_events[0][0][1]
    lnrho_R = brentq(resid, np.log(c['rho_lo']), np.log(c['rho_hi']), xtol=1e-13)
    su = shoot(lnrho_R)
    rfl = su.t_events[0][0]

    def evj(rr, y):
        return y[1] - c['tauj']
    evj.terminal = True
    s1 = solve_ivp(rhs1, [R, 0.5*R], [lnrho_R, 2.0/3.0], events=evj, method='LSODA',
                   rtol=rtol, atol=1e-13, dense_output=True)
    rj = s1.t_events[0][0]
    yj = s1.y_events[0][0]
    Tj = (3.0*(LUM/(4*np.pi*rj**2))/(4.0*SIG)*(yj[1] + q_hopf(yj[1])))**0.25
    rhoj = np.exp(yj[0])
    rbot = c['r_bot1'] if c.get('legacy_deep') else c['r_bot']
    rhs2 = _rhs_interior(st)
    kw = {}
    if c.get('tstop'):                  # explore only (the build passes no events)
        def evT(rr, y):
            return y[0] - np.log(c['tstop'])
        evT.terminal = True
        kw['events'] = [evT]
    s2 = solve_ivp(rhs2, [rj, rbot], [np.log(Tj), np.log(rhoj)], method='LSODA', rtol=rtol,
                   atol=1e-12, dense_output=True, **kw)
    assert s2.success, s2.message
    rend = s2.t[-1]
    # ---- sample (make_ic_agcar_ge.sample)
    rwall = c['r_wall1'] if c.get('legacy_deep') else c['r_wall']
    rg = np.exp(np.arange(np.log(rbot), np.log(rj), c['dlnr']))
    rg = np.unique(np.concatenate([rg, [rwall]]))
    drn = c['drn_atm']
    r_in = rg[(rg < rj) & (rg >= rend)]
    y2 = s2.sol(r_in)
    r_h = np.arange(rj, R, drn)
    y1 = s1.sol(r_h)
    tau_h = np.maximum(y1[1], 1e-8)
    T_h = (3.0*(LUM/(4*np.pi*r_h**2))/(4.0*SIG)*(tau_h + q_hopf(tau_h)))**0.25
    r_u = np.arange(R, rfl, drn)
    yu = su.sol(r_u)
    rtop = c['r_col_top']
    r_f = np.arange(r_u[-1] + drn, rtop + drn*c['fl_step'], drn*c['fl_step'])
    r = np.concatenate([r_in, r_h, r_u, r_f])
    T = np.concatenate([np.exp(y2[0]), T_h, np.full(r_u.size + r_f.size, Tph)])
    rho = np.concatenate([np.exp(y2[1]), np.exp(y1[0]), np.exp(yu[0]),
                          np.full(r_f.size, rho_floor)])
    ffr, fm = _mlt_cols(st, r, T, rho, r < rj)
    lt, lr = np.log10(T), np.log10(rho)
    eint = rho*10.0**eos.le(lt, lr)
    Fr = ffr*(LUM/(4*np.pi*r**2))
    E = _ecol(st, r, T, R, Tph)
    out = dict(r=r, rho=rho, T=T, eint=eint, Fr=Fr, E=E, fmlt=fm, rph=R, tph=Tph,
               rho_R=np.exp(lnrho_R), rj=rj, r_floor=rfl, scheme='agcar_ge', n_deep=0)
    if c.get('legacy_deep'):
        out = deep_continue(st, out, c)
    return out


def deep_continue(st, col, c):
    """make_ic_agcar_deep.py: rows of the column kept verbatim (as printed), new rows below
    its first row on the continued log grid (rtol 1e-10, atol 1e-13) down to c['r_bot']."""
    LUM, eos = st.LUM, st.eos
    row0 = [float('%.10e' % col['r'][0]), float('%.10e' % col['rho'][0]),
            float('%.10e' % col['T'][0])]
    r0, rho0, T0 = row0
    rhs2 = _rhs_interior(st)
    dl = c['dlnr']
    nk = int(np.floor(np.log(r0/c['r_bot'])/dl)) + 1
    rn = r0*np.exp(-dl*np.arange(nk + 1))
    s = solve_ivp(rhs2, [r0, rn[-1]], [np.log(T0), np.log(rho0)], t_eval=rn, method='LSODA',
                  rtol=1e-10, atol=1e-13)
    assert s.success, s.message
    r = s.t[1:][::-1]
    T = np.exp(s.y[0][1:])[::-1]
    rho = np.exp(s.y[1][1:])[::-1]
    ffr, fm = _mlt_cols(st, r, T, rho, np.ones(r.size, bool))
    eint = rho*10.0**eos.le(np.log10(T), np.log10(rho))
    Fr = ffr*LUM/(4*np.pi*r**2)
    E = A_RAD_*T**4
    out = dict(col)
    for k, v in (('r', r), ('rho', rho), ('T', T), ('eint', eint), ('Fr', Fr), ('E', E),
                 ('fmlt', fm)):
        out[k] = np.concatenate([v, col[k]])
    out['n_deep'] = r.size
    out['scheme'] = 'agcar_ge+legacy_deep'
    return out


# --------------------------------------------------------------------------- compatibility
def column_mltpp_inward(st, c):
    """make_ic_agcar_mltpp.py --mode inward (Raven icB2-1011 reconstruction), bit for bit."""
    LUM, RSTAR = st.LUM, st.RSTAR
    eos = st.eos
    DRN = c['dlnr']*RSTAR
    rhs1 = _rhs_hopf(st)
    Tph = (3.0*(LUM/(4*np.pi*RSTAR**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
    rhs_iso = _rhs_iso(st, Tph, c['atm_rad'])
    rho_fl = c['rho_floor']

    def top_event(rr, y):
        return y[0] - np.log(rho_fl)
    top_event.terminal = True

    def blow_event(rr, y):
        return y[0] - np.log(1.0e-6)
    blow_event.terminal = True

    def shoot(lnrho_R):
        return solve_ivp(rhs_iso, [RSTAR, RSTAR + 2.0*RSTAR], [lnrho_R, 2.0/3.0],
                         events=[top_event, blow_event], method='LSODA', rtol=1e-10,
                         atol=1e-13)

    def resid(lnrho_R):
        s = shoot(lnrho_R)
        if s.t_events[0].size == 0:
            return -1.0
        return s.y_events[0][0][1]
    lo, hi = np.log(c['rho_lo']), np.log(c['rho_hi'])
    assert resid(lo)*resid(hi) < 0
    lnrho_R = brentq(resid, lo, hi, xtol=1e-12)
    su = shoot(lnrho_R)
    rtop_atm = su.t_events[0][0]
    rx = np.arange(RSTAR, rtop_atm, DRN)
    su = solve_ivp(rhs_iso, [RSTAR, rtop_atm], [lnrho_R, 2.0/3.0], t_eval=rx,
                   method='LSODA', rtol=1e-10, atol=1e-13)
    assert su.success, su.message
    r_f = np.arange(rx[-1] + DRN, c['r_col_top'], DRN*c['fl_step'])
    r_u = np.concatenate([su.t[1:], r_f])
    rho_u = np.concatenate([np.exp(su.y[0][1:]), np.full(r_f.size, rho_fl)])
    T_u = np.full(r_u.size, Tph)
    rho_R = np.exp(lnrho_R)

    def evj(rr, y):
        return y[1] - c['tauj']
    evj.terminal = True
    rd = np.arange(RSTAR, 0.9*RSTAR, -DRN)
    s1 = solve_ivp(rhs1, [RSTAR, 0.9*RSTAR], [np.log(rho_R), 2.0/3.0], t_eval=rd,
                   events=evj, method='LSODA', rtol=1e-10, atol=1e-13, dense_output=True)
    rj = s1.t_events[0][0]
    yj = s1.y_events[0][0]
    Tj = (3.0*(LUM/(4*np.pi*rj**2))/(4.0*SIG)*(yj[1] + q_hopf(yj[1])))**0.25
    rhoj = np.exp(yj[0])
    rhs2 = _rhs_interior(st, with_tau=True)
    nk = int(np.floor(np.log(rj/c['r_bot'])/c['dlnr'] + 1e-9))
    rn = rj*np.exp(-c['dlnr']*np.arange(nk + 1))
    rn = rn[rn >= c['r_bot']]
    s2 = solve_ivp(rhs2, [rj, rn[-1]], [np.log(Tj), np.log(rhoj), yj[1]], t_eval=rn,
                   method='LSODA', rtol=1e-9, atol=1e-12)
    assert s2.success, s2.message
    r_h = rj + DRN*np.arange(1, int(np.ceil((RSTAR - rj)/DRN)) + 1)
    r_h = np.concatenate([r_h[r_h < RSTAR*(1.0 - 1e-12)], [RSTAR]])
    yh = s1.sol(r_h)
    yh[:, -1] = [np.log(rho_R), 2.0/3.0]
    tau_h = np.maximum(yh[1], 1e-8)
    T_h = (3.0*(LUM/(4*np.pi*r_h**2))/(4.0*SIG)*(tau_h + q_hopf(tau_h)))**0.25
    rho_h = np.exp(yh[0])
    r = np.concatenate([s2.t[::-1], r_h, r_u])
    T = np.concatenate([np.exp(s2.y[0])[::-1], T_h, T_u])
    rho = np.concatenate([np.exp(s2.y[1])[::-1], rho_h, rho_u])
    ffr, fm = _mlt_cols(st, r, T, rho, r < rj)
    lt, lr = np.log10(T), np.log10(rho)
    eint = rho*10.0**eos.le(lt, lr)
    Fr = ffr*LUM/(4*np.pi*r**2)
    F = LUM/(4*np.pi*r**2)                       # (make_ic_agcar_mltpp.ecol, its op order)
    E = A_RAD_*T**4
    up = r > RSTAR
    Fph = LUM/(4*np.pi*RSTAR**2)
    E[up] = F[up]/CL*(1.0 + (CL*A_RAD_*Tph**4/Fph - 1.0)*(RSTAR/r[up])**2)
    return dict(r=r, rho=rho, T=T, eint=eint, Fr=Fr, E=E, fmlt=fm, rph=RSTAR, tph=Tph,
                rho_R=rho_R, rj=rj, r_floor=rtop_atm, scheme='mltpp_inward', n_deep=0)


def explore(st, c, r_stop, tstop):
    """standard equations, one dense interior solve r_j -> r_stop or T = tstop (for the
    r_in rule); coarser rows (dlnr x 5)."""
    c2 = dict(c)
    c2['r_bot'] = r_stop
    c2['r_wall'] = r_stop
    c2['tstop'] = tstop
    c2['dlnr'] = 5*c['dlnr']
    c2['legacy_deep'] = False
    return column_agcar_ge(st, c2)


def build(st, c):
    if c['scheme'] == 'agcar_ge':
        return column_agcar_ge(st, c)
    if c['scheme'] == 'mltpp_inward':
        return column_mltpp_inward(st, c)
    raise ValueError('unknown scheme ' + c['scheme'])


# --------------------------------------------------------------------------- MESA-gated MLT++
def gating_rule(st, col, xa1c=0.0, ctrl=None):
    """MESA r15140 set_gradT_excess_alpha (mlt_info.f90:1330-1489) on a built column, over
    the rows below R_ph: beta = min P_gas/P x (1 + xa1c), lambda = min(1, max L_rad/L_Edd =
    kappa F_r/(c g) from R_ph down to the beta-minimum row) -> alpha in the (beta, lambda)
    plane.  Returns alpha, beta, r(beta min), lambda."""
    from .micro import gradT_excess_alpha, PP_CTRL
    r, T, rho = col['r'], col['T'], col['rho']
    m = r <= col['rph']*(1 + 1e-12)
    P, pg = st.thermo(T[m], rho[m])[:2]
    b = pg/P
    kb = int(np.argmin(b))
    beta = float(b[kb])*(1.0 + xa1c)
    sel = r[m] >= r[m][kb]
    kap = st.kap(T[m], rho[m])
    g = st.GM/r[m]**2
    lam = min(1.0, float(np.max((kap*col['Fr'][m]/(CL*g))[sel])))
    c = dict(PP_CTRL)
    if ctrl:
        c.update(ctrl)
    return gradT_excess_alpha(beta, lam, c), beta, float(r[m][kb]), lam


def solve_gated(st, c, r_stop, tstop, xa1c, ctrl, tol, log):
    """alpha* = rule(column(alpha*)) (bisection on [0, 1], as the AG Car A fixed point);
    coarse columns (dlnr 1e-3, rtol 1e-9) to r_stop or T = tstop."""
    def rule(a):
        st.PP_ALPHA = a
        cc = dict(c, scheme='agcar_ge', legacy_deep=False, dlnr=1e-3, rtol=1e-9,
                  r_bot=r_stop, r_wall=r_stop, tstop=tstop)
        col = column_agcar_ge(st, cc)
        al, be, rb, lam = gating_rule(st, col, xa1c, ctrl)
        log('    MLT++ gating: alpha used %.5f -> beta_min %.4f at %.3f Rsun, lambda %.4f -> '
            'rule %.5f' % (a, be, rb/6.957e10, lam, al))
        return al
    r0 = rule(0.0)
    if r0 <= 0.0:
        st.PP_ALPHA = 0.0
        return 0.0
    r1 = rule(1.0)
    if r1 >= 1.0:
        st.PP_ALPHA = 1.0
        return 1.0
    lo, hi = 0.0, 1.0
    while hi - lo > tol:
        mid = 0.5*(lo + hi)
        if rule(mid) > mid:
            lo = mid
        else:
            hi = mid
    st.PP_ALPHA = 0.5*(lo + hi)
    return st.PP_ALPHA
