"""Mode 'mesa' of the standard IC (tools/ic_star; USER decision 10-15, option A).

When a MESA profile exists the STRUCTURE is MESA's, the microphysics ours:
  * T(P): MESA's ln T as a function of ln P_tot (PCHIP), i.e. MESA's actual gradient
    nabla_MESA = dlnT/dlnP, as a function of pressure;
  * HSE re-integrated OUTWARD from the column bottom with OUR EOS (rho from
    P - a T^4/3 = P_gas,code(rho, T)) and the code's gravity: G m(r) = G m(r_bot) +
    G int 4 pi r^2 rho dr (deck: he_gm = G m(x1min) of this column, he_gm_column = true);
    start P(r_bot) = P_MESA(r_bot), m(r_bot) = m_MESA(r_bot);
  * stop where T(P) reaches the Hopf tau = tau_j temperature (3 F/(4 sigma)(tau_j + q))^1/4,
    F = L/(4 pi r^2); above: the Hopf layer integrated OUTWARD (gas HSE with the radiation
    force of F, tau from tau_j down to 2/3) -> OUR R_ph and rho_ph; above R_ph the standard
    isothermal gravity-only atmosphere at T_ph from rho_ph down to the column floor, floor
    to the column top (as mode 'global');
  * F_rad from OUR kappa_R: F_rad/F = nabla/nabla_rad(ours, g = G m/r^2, F = L/(4 pi r^2));
    F_MLT/F = 1 - F_rad/F CLIPPED to [0, 1]; where F_rad > F (our kappa lower than MESA's)
    the column's F_r = F (no negative scaffold) and the excess F_rad/F - 1 is reported:
    those layers will LOSE heat in the run (radiation alone carries more than L there).
L = L_ph of the MESA header (constant; L_MESA(r)/L_ph at the wall is checked).
MESA's own flux split for the CZ comparison: 'mesa' format F_rad,M/F = nabla_M/nabla_rad,M
(MESA opacity, L(r), m(r)); 'cgs6' format F_rad,M/F = 4 pi r^2 F_r,M/L.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from .column import _rhs_hopf, _rhs_iso, _ecol
from .micro import A_RAD_, CL, G_CODE, SIG, q_hopf


def mesa_struct(st, mesa):
    """MESA arrays at the zone centres: r, T, rho, P_tot, m, L, kappa, nabla_M, F_rad,M/F,
    H_p.  cgs6 (no P column): P from MESA's own HSE integrated inward from the top row,
    P_top = P_code(rho, T) + a T^4/3 there."""
    r = mesa['rc']
    T, rho = mesa['T'], mesa['rho']
    if mesa['fmt'] == 'mesa':
        P = mesa['P']
        mc = np.interp(r, mesa['r'], mesa['m'])
        Lc = np.interp(r, mesa['r'], mesa['L'])
    else:
        mc, Lc = mesa['m'], mesa['L']
        P = np.empty_like(r)
        lt, lr = np.log10(T[-1]), np.log10(rho[-1])
        P[-1] = rho[-1]*10**st.eos.lp(lt, lr) + A_RAD_*T[-1]**4/3.0
        gr = G_CODE*mc*rho/r**2
        for i in range(r.size - 2, -1, -1):
            P[i] = P[i+1] + 0.5*(gr[i] + gr[i+1])*(r[i+1] - r[i])
    lnP, lnT = np.log(P), np.log(T)
    nab = np.gradient(lnT, lnP)
    g = G_CODE*mc/r**2
    if mesa['fmt'] == 'mesa':
        kap = mesa['kap']
        nrad = 3*kap*P*Lc/(16*np.pi*A_RAD_*CL*G_CODE*mc*T**4)
        frad = np.minimum(nab/nrad, 1.0)
    else:
        kap = mesa['kap']
        frad = np.minimum(mesa['L']/mesa['L_ph'], 1.0)
    return dict(r=r, T=T, rho=rho, P=P, m=mc, L=Lc, kap=kap, nab=nab, frad=frad,
                Hp=P/(rho*g), fmlt=1.0 - frad)


def column_mesa(st, c, mesa, ms, log=print):
    """the column; c as mode 'global' (+ r_bot, r_wall, r_col_top, rho_floor rule)"""
    LUM = st.LUM
    eos = st.eos
    o = np.argsort(np.log(ms['P']))
    lnPs, lnTs = np.log(ms['P'])[o], np.log(ms['T'])[o]
    keep = np.concatenate([[True], np.diff(lnPs) > 0])
    fT = PchipInterpolator(lnPs[keep], lnTs[keep])
    dT = fT.derivative()
    lnPmin = lnPs[keep][0]
    r0 = c['r_bot']
    P0 = float(np.exp(np.interp(np.log(r0), np.log(ms['r']), np.log(ms['P']))))
    m0 = float(np.interp(r0, ms['r'], ms['m']))*c.get('gm_scale', 1.0)   # (diagnostic variants)
    guess = [float(np.exp(np.interp(np.log(r0), np.log(ms['r']), np.log(ms['rho']))))]

    def rho_of(P, T):
        lt = np.log10(T)
        pg = P - A_RAD_*T**4/3.0
        if pg <= 0:
            raise RuntimeError('P_gas <= 0 (radiation-dominated beyond the table)')
        lpg = np.log10(pg)
        x = np.log10(guess[0])
        for _ in range(30):
            f = x + eos.lp(lt, x) - lpg
            dx = -f/(1.0 + eos.lp(lt, x, dy=1))
            x += dx
            if abs(dx) < 1e-13:
                break
        guess[0] = 10.0**x
        return 10.0**x

    tauj = c['tauj']

    def rhs(rr, y):
        lnP, m = y
        P = np.exp(lnP)
        T = np.exp(fT(lnP))
        rho = rho_of(P, T)
        return [-G_CODE*m*rho/(rr**2*P), 4*np.pi*rr**2*rho]

    def ev_join(rr, y):
        Tj = (3.0*(LUM/(4*np.pi*rr**2))/(4.0*SIG)*(tauj + q_hopf(tauj)))**0.25
        return np.exp(fT(y[0])) - Tj
    ev_join.terminal = True

    def ev_pmin(rr, y):
        return y[0] - lnPmin
    ev_pmin.terminal = True
    rmax = 3.0*ms['r'][-1]
    s = solve_ivp(rhs, [r0, rmax], [np.log(P0), m0], method='LSODA', rtol=1e-10,
                  atol=[1e-12, m0*1e-12], dense_output=True,
                  events=[ev_pmin] if c.get('mesa_top') == 'surface' else [ev_join, ev_pmin])
    dl = c['dlnr']
    ext = None
    surf = c.get('mesa_top') == 'surface'
    if surf:
        # 10-16 mesa_top = surface: MESA's T(P) all the way to its surface zone (its
        # photosphere); R_ph = that radius, the iso atmosphere directly on top (no Hopf layer).
        # The scaffold then reaches the photosphere as in MESA (F_MLT/F of our kappa on MESA's
        # T gradient); the taper (extend) removes it above.
        if not s.t_events[0].size:
            raise RuntimeError('mode mesa (surface): the MESA profile top was not reached')
        rj = s.t_events[0][0]
        yj = s.y_events[0][0]
        Tj = float(np.exp(fT(yj[0])))
        guess[0] = float(np.exp(np.interp(np.log(min(rj, ms['r'][-1])), np.log(ms['r']),
                                          np.log(ms['rho']))))
        rhoj = rho_of(np.exp(yj[0]), Tj)
        mj = yj[1]
        r_e = rj
    elif s.t_events[0].size:
        rj = s.t_events[0][0]
        yj = s.y_events[0][0]
        Tj = float(np.exp(fT(yj[0])))
        guess[0] = float(np.exp(np.interp(np.log(min(rj, ms['r'][-1])), np.log(ms['r']),
                                          np.log(ms['rho']))))
        rhoj = rho_of(np.exp(yj[0]), Tj)
        mj = yj[1]
        r_e = rj
    elif s.t_events[1].size and c.get('top_fit', 'two_sided') == 'two_sided':
        return two_sided(st, c, mesa, ms, s, fT, dT, rho_of, guess, log)
    elif s.t_events[1].size:
        # the MESA profile ends below tau_j (e.g. L. Ma's BSG file, top row T 63 kK): continue
        # OUTWARD with the mode-'global' interior equations (local MLT++ with OUR kappa,
        # dP_tot/dr = -rho G m/r^2) to the tau_j temperature
        r_e = s.t_events[1][0]
        ye = s.y_events[1][0]
        Te = float(np.exp(fT(ye[0])))
        guess[0] = float(np.exp(np.interp(np.log(min(r_e, ms['r'][-1])), np.log(ms['r']),
                                          np.log(ms['rho']))))
        rhoe = rho_of(np.exp(ye[0]), Te)
        gm_save = st.GM

        def rhs_x(rr, y):
            T_, rho_, m_ = np.exp(y[0]), np.exp(y[1]), y[2]
            st.GM = G_CODE*m_
            nab_ = st.mlt_gradient(T_, rho_, rr)[0]
            P_, pg_, chir, chit, cp, dlt, nad_ = st.thermo(T_, rho_)
            dlnP = -rho_*(G_CODE*m_/rr**2)/P_
            dlnT = nab_*dlnP
            return [dlnT, (dlnP - chit*dlnT)/chir, 4*np.pi*rr**2*rho_]

        def ev_j2(rr, y):
            Tj_ = (3.0*(LUM/(4*np.pi*rr**2))/(4.0*SIG)*(tauj + q_hopf(tauj)))**0.25
            return y[0] - np.log(Tj_)
        ev_j2.terminal = True

        def ev_stall(rr, y):
            return y[1] - (np.log(rhoe) - 12.0)
        ev_stall.terminal = True

        def run(eps):
            """outward extension from rho_e (1 + eps); -> (sol, r_j, T_j, rho_j, m_j, R_ph)
            or None (stalls: rho -> 0, the Gamma -> 1 separatrix)"""
            try:
                sx_ = solve_ivp(rhs_x, [r_e, 1.5*r_e], [np.log(Te), np.log(rhoe*(1 + eps)),
                                                        ye[1]], method='LSODA', rtol=1e-10,
                                atol=1e-12, dense_output=True, events=[ev_j2, ev_stall])
            finally:
                st.GM = gm_save
            if sx_.t_events[0].size == 0:
                return None
            y_ = sx_.y_events[0][0]
            rj_ = sx_.t_events[0][0]
            R_ = hopf_out(st, rj_, float(np.exp(y_[1])), tauj, G_CODE*y_[2], c['rtol'])[0]
            if R_ is None:
                return None
            return sx_, rj_, float(np.exp(y_[0])), float(np.exp(y_[1])), float(y_[2]), R_
        Rt = c['R_target']
        lo, hi = c.get('eps_lo', -0.5), c.get('eps_hi', 2.0)
        for e_ in (-0.5, -0.1, 0.0, 0.1, 0.5, 1.0, 2.0):
            rr_ = run(e_)
            log('    extension eps %+.2f -> %s' % (e_, 'stalls' if rr_ is None else
                                                  'r_j %.5e R_ph %.5e (target %.5e)'
                                                  % (rr_[1], rr_[5], Rt)))
        rh = run(hi)
        if rh is None or rh[5] > Rt:
            raise RuntimeError('mode mesa: extension shot: eps_hi %g does not give R_ph < target'
                               % hi)
        rl = run(lo)
        if rl is not None and rl[5] < Rt:
            raise RuntimeError('mode mesa: extension shot: eps_lo %g gives R_ph < target' % lo)
        for _ in range(60):
            mid = 0.5*(lo + hi)
            rm = run(mid)
            log('    shot eps %+.8e -> %s' % (mid, 'stalls' if rm is None else 'R_ph %.6e' % rm[5]))
            if rm is None or rm[5] > Rt:
                lo = mid
            else:
                hi, rh = mid, rm
            if hi - lo < 1e-9 or abs(rh[5]/Rt - 1) < 1e-7:
                break
        eps = hi
        sx, rj, Tj, rhoj, mj = rh[:5]
        log('  mode mesa: the MESA profile ends at P %.4e (r %.6e = %.3f Rsun, T %.4e K) below '
            'tau_j: OUTWARD extension (local MLT/MLT++, OUR kappa, G m(r)) from rho_e (1 + eps),'
            ' eps SHOT to R_ph = %.6e: eps = %+.6e (the density jump at the profile top); '
            'R_ph %.6e' % (np.exp(ye[0]), r_e, r_e/6.957e10, Te, Rt, eps, rh[5]))
        ext = dict(r_e=r_e, rj=rj, sol=sx, eps=eps)
    else:
        raise RuntimeError('mode mesa: neither tau_j nor the profile top reached')
    GMj = G_CODE*mj
    log('  mode mesa: tau-%g join at r_j %.8e (%.4f R_MESA), T %.2f K, rho %.4e, m %.6e'
        % (tauj, rj, rj/mesa['R_ph'] if mesa['R_ph'] else np.nan, Tj, rhoj, mj))
    # interior rows (MESA T(P) part, then the extension)
    rg = np.exp(np.arange(np.log(r0), np.log(rj), dl))
    rg = np.unique(np.concatenate([rg, [c['r_wall']]]))
    rg = rg[rg < rj]
    ra = rg[rg <= r_e]
    yi = s.sol(ra)
    Ti = np.exp(fT(yi[0]))
    guess[0] = float(np.exp(np.interp(np.log(r0), np.log(ms['r']), np.log(ms['rho']))))
    rhoi = np.array([rho_of(np.exp(a), b) for a, b in zip(yi[0], Ti)])
    mi = yi[1]
    nab = dT(yi[0])
    Pi = np.exp(yi[0])
    kap = st.kap(Ti, rhoi)
    F = LUM/(4*np.pi*ra**2)
    nrad = 3.0*kap*Pi*F/(4.0*A_RAD_*CL*(G_CODE*mi/ra**2)*Ti**4)
    fr = nab/nrad
    if ext is not None:
        rb = rg[rg > r_e]
        yb = ext['sol'].sol(rb)
        Tb, rhob, mb = np.exp(yb[0]), np.exp(yb[1]), yb[2]
        frb, nabb, nradb = np.empty(rb.size), np.empty(rb.size), np.empty(rb.size)
        gm_save = st.GM
        for k in range(rb.size):
            st.GM = G_CODE*mb[k]
            o_ = st.mlt_gradient(Tb[k], rhob[k], rb[k])
            frb[k], nabb[k], nradb[k] = o_[1], o_[0], o_[8]
        st.GM = gm_save
        Ti, rhoi, mi = np.r_[Ti, Tb], np.r_[rhoi, rhob], np.r_[mi, mb]
        fr, nab, nrad = np.r_[fr, frb], np.r_[nab, nabb], np.r_[nrad, nradb]
    F = LUM/(4*np.pi*rg**2)
    excess = fr > 1.0
    fm = np.clip(1.0 - fr, 0.0, 1.0)
    Fr = F*np.minimum(fr, 1.0)
    if surf:
        R, rho_ph, Tph = rj, rhoj, Tj
        drn = c['drn_atm_rph']*R
        r_h, T_h, rho_h = np.array([]), np.array([]), np.array([])
    else:
        # Hopf layer OUTWARD from (r_j, rho_j, tau_j) to tau 2/3
        R, rho_ph, s1 = hopf_out(st, rj, rhoj, tauj, GMj, c['rtol'])
        if R is None:
            raise RuntimeError('mode mesa: the Hopf layer never reached tau 2/3')
        Tph = (3.0*(LUM/(4*np.pi*R**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
        drn = c['drn_atm_rph']*R
        r_h = rj + drn*np.arange(1, int(np.ceil((R - rj)/drn)) + 1)
        r_h = np.concatenate([r_h[r_h < R*(1.0 - 1e-9)], [R]])
        yh = s1.sol(r_h)
        tau_h = np.maximum(yh[1], 1e-8)
        T_h = (3.0*(LUM/(4*np.pi*r_h**2))/(4.0*SIG)*(tau_h + q_hopf(tau_h)))**0.25
        rho_h = np.exp(yh[0])
    # column floor (rule needs rho_ph) and the isothermal atmosphere above R_ph
    rho_floor = c['rho_floor'] if c.get('rho_floor') else \
        10.0**np.round(np.log10(c['rho_floor_rel']*rho_ph))
    rhs_iso = _rhs_iso(st, Tph, c['atm_rad'], GMj)

    def top_event(rr, y):
        return y[0] - np.log(rho_floor)
    top_event.terminal = True
    su = solve_ivp(rhs_iso, [R, 3.0*R], [np.log(rho_ph), 2.0/3.0], events=[top_event],
                   method='LSODA', rtol=c['rtol'], atol=1e-13, dense_output=True)
    rfl = su.t_events[0][0] if su.t_events[0].size else 3.0*R
    r_u = np.arange(R + drn, rfl, drn)
    yu = su.sol(r_u)
    r_f = np.arange((r_u[-1] if r_u.size else R) + drn, c['r_col_top'] + drn*c['fl_step'],
                    drn*c['fl_step'])
    r = np.concatenate([rg, r_h, r_u, r_f])
    T = np.concatenate([Ti, T_h, np.full(r_u.size + r_f.size, Tph)])
    rho = np.concatenate([rhoi, rho_h, np.exp(yu[0]), np.full(r_f.size, rho_floor)])
    Fr_all = np.concatenate([Fr, LUM/(4*np.pi*np.concatenate([r_h, r_u, r_f])**2)])
    fm_all = np.concatenate([fm, np.zeros(r.size - fm.size)])
    eint = rho*10.0**eos.le(np.log10(T), np.log10(rho))
    E = _ecol(st, r, T, R, Tph)
    # mass of our column, the wall mass
    mw = float(np.interp(c['r_wall'], rg, mi))
    out = dict(r=r, rho=rho, T=T, eint=eint, Fr=Fr_all, E=E, fmlt=fm_all, rph=R, tph=Tph,
               rho_R=rho_ph, rj=rj, r_floor=rfl, scheme='mesa', n_deep=0,
               rho_floor=rho_floor, m_wall=mw, m_top=float(mj), GM_wall=G_CODE*mw,
               excess=dict(n=int(excess.sum()),
                           max=float((fr[excess] - 1).max()) if excess.any() else 0.0,
                           rlo=float(rg[excess].min()) if excess.any() else np.nan,
                           rhi=float(rg[excess].max()) if excess.any() else np.nan),
               interior=dict(r=rg, m=mi, nab=nab, nrad=nrad, frad=fr),
               ext=None if ext is None else dict(r_e=ext['r_e'], rj=ext['rj'], eps=ext['eps']))
    return out


def mesa_rin(ms, n_hp, cz_fmlt=0.01):
    """the r_in rule on MESA's structure: n_hp H_p below the base of MESA's Fe CZ (the
    contiguous F_MLT,M/F > cz_fmlt region containing the max of kappa_M at 1e5-1e6 K; the
    Fe peak itself if MESA is radiative there)"""
    r, T, kap, fm, Hp = ms['r'], ms['T'], ms['kap'], ms['fmlt'], ms['Hp']
    sel = (T > 1e5) & (T < 1e6)
    ip = np.where(sel)[0][np.argmax(kap[sel])]
    cz = fm > cz_fmlt
    a = ip
    if cz[ip]:
        while a > 0 and cz[a-1]:
            a -= 1
    rb = r[a]
    x, h = r[:a+1][::-1], Hp[:a+1][::-1]
    acc = 0.0
    for k in range(1, x.size):
        d = 0.5*(1/h[k] + 1/h[k-1])*(x[k-1] - x[k])
        if acc + d >= n_hp:
            rin = x[k-1] - (n_hp - acc)/(d/(x[k-1] - x[k]))
            return rin, dict(r_fe=r[ip], fe_convective=bool(cz[ip]), r_base=rb)
        acc += d
    return None, dict(r_fe=r[ip], fe_convective=bool(cz[ip]), r_base=rb)


def mesa_zones(ms, cz_fmlt=0.01):
    r, fm = ms['r'], ms['fmlt']
    zs, i = [], 0
    while i < r.size:
        if fm[i] > cz_fmlt:
            k = i
            while k + 1 < r.size and fm[k+1] > cz_fmlt:
                k += 1
            zs.append((r[i], r[k], float(fm[i:k+1].max())))
            i = k + 1
        else:
            i += 1
    return zs


def hopf_out(st, rj, rhoj, tauj, GM, rtol):
    """the Hopf layer (gas HSE with the radiation force of F) OUTWARD from (r_j, rho_j,
    tau_j) to tau = 2/3 -> (R_ph, rho_ph, solution) or (None, None, None)"""
    rhs1 = _rhs_hopf(st, GM)

    def ev_ph(rr, y):
        return y[1] - 2.0/3.0
    ev_ph.terminal = True
    s1 = solve_ivp(rhs1, [rj, 2.0*rj], [np.log(rhoj), tauj], events=[ev_ph], method='LSODA',
                   rtol=rtol, atol=1e-13, dense_output=True)
    if s1.t_events[0].size == 0:
        return None, None, None
    return s1.t_events[0][0], float(np.exp(s1.y_events[0][0][0])), s1


def two_sided(st, c, mesa, ms, s, fT, dT, rho_of, guess, log):
    """10-16: the MESA profile is truncated below tau_j (L. Ma's BSG file).  Two-sided fit
    instead of a density jump:
      outward: MESA T(P) + our HSE (G m(r)) from the column bottom to the profile top r_e;
      inward:  the mode-global surface column (Hopf + iso atmosphere, local MLT of the spec,
               MESA defaults = plain Cox MLT alpha 2) from R_ph down into [0.9 r_e, r_e],
               point mass G m(r_e);
    FREE PARAMETER: R_ph (within c['fit_dR'] = 4 % of the target; fit points in [0.75, 1] r_e).  For each R_ph the fit point
    r* is where P_in = P_out (P continuous); R_ph is chosen so that T is continuous there too
    (then rho is continuous through the EOS).  If no root: the R_ph with the smallest |dlnT|."""
    from scipy.optimize import brentq
    from .column import column_agcar_ge
    LUM = st.LUM
    r_e = s.t_events[1][0]
    m_e = float(s.y_events[1][0][1])
    gm_save, R_save = st.GM, st.RSTAR
    rlo = c.get('fit_rlo', 0.75)*r_e
    rr = np.linspace(rlo, r_e*(1 - 1e-9), 600)
    yo = s.sol(rr)
    lnPo, lnTo = yo[0], fT(yo[0])
    Rt = c['R_target']

    def track(R, dlnr):
        st.RSTAR, st.GM = R, G_CODE*m_e
        cc = dict(c, scheme='agcar_ge', legacy_deep=False, r_bot=rlo*0.995, r_wall=rlo*0.995,
                  dlnr=dlnr, drn_atm=c['drn_atm_rph']*R,
                  rho_floor=c.get('rho_floor') or 1e-16)
        try:
            return column_agcar_ge(st, cc)
        finally:
            st.GM = gm_save

    def resid(R, dlnr=1e-3, full=False):
        try:
            col = track(R, dlnr)
        except Exception as e:                           # noqa: BLE001
            log('    two-sided fit: the inward column for R_ph %.6e failed (%s)' % (R, str(e)[:80]))
            return (None, None, None) if full else np.nan
        m = (col['r'] >= rlo*0.995) & (col['r'] < col['rj'])
        P = st.thermo(col['T'][m], col['rho'][m])[0]
        lPi = np.interp(rr, col['r'][m], np.log(P))
        lTi = np.interp(rr, col['r'][m], np.log(col['T'][m]))
        dP = lPi - lnPo
        k = np.where(np.sign(dP[1:]) != np.sign(dP[:-1]))[0]
        if k.size == 0:
            return (None, None, col) if full else np.nan
        k = k[-1]
        w = dP[k]/(dP[k] - dP[k+1])
        rs = rr[k] + w*(rr[k+1] - rr[k])
        dTs = (1 - w)*(lTi[k] - lnTo[k]) + w*(lTi[k+1] - lnTo[k+1])
        return (rs, dTs, col) if full else dTs
    dR = c.get('fit_dR', 0.04)
    Rs = Rt*(1 + dR*np.linspace(-1, 1, 17))
    fs = [resid(R) for R in Rs]
    for R, f in zip(Rs, fs):
        log('    two-sided fit scan: R_ph %.6e (%.4f target): dlnT at the P-match point %s'
            % (R, R/Rt, 'no P match in [%g, 1] r_e' % (rlo/r_e) if f != f else '%+.4e' % f))
    ok = [i for i in range(len(fs) - 1) if fs[i] == fs[i] and fs[i+1] == fs[i+1] and
          np.sign(fs[i]) != np.sign(fs[i+1])]
    if ok:
        i = min(ok, key=lambda j: abs(Rs[j] - Rt))
        Rbest = brentq(lambda R: resid(R, 2e-4), Rs[i], Rs[i+1], xtol=1e-7*Rt)
        how = 'root'
    else:
        good = [i for i in range(len(fs)) if fs[i] == fs[i]]
        if not good:
            st.RSTAR = R_save
            raise RuntimeError('two-sided fit: no P match for R_ph within +-%g' % dR)
        i = min(good, key=lambda j: abs(fs[j]))
        Rbest, how = Rs[i], 'min |dlnT| (no root within +-%g)' % dR
    if not c.get('rho_floor'):          # the floor rule on the fitted photosphere
        c['rho_floor'] = 10.0**np.round(np.log10(c['rho_floor_rel']*track(Rbest, 1e-3)['rho_R']))
    rs, dTs, col = resid(Rbest, c['dlnr'], full=True)
    # the outward MESA part below r*, the inward column at and above r*
    dl = c['dlnr']
    r0 = c['r_bot']
    rg = np.exp(np.arange(np.log(r0), np.log(rs), dl))
    rg = np.unique(np.concatenate([rg, [c['r_wall']]]))
    rg = rg[rg < rs]
    yi = s.sol(rg)
    Ti = np.exp(fT(yi[0]))
    guess[0] = float(np.exp(np.interp(np.log(r0), np.log(ms['r']), np.log(ms['rho']))))
    rhoi = np.array([rho_of(np.exp(a), b) for a, b in zip(yi[0], Ti)])
    mi = yi[1]
    nab = dT(yi[0])
    Pi = np.exp(yi[0])
    kap = st.kap(Ti, rhoi)
    F = LUM/(4*np.pi*rg**2)
    nrad = 3.0*kap*Pi*F/(4.0*A_RAD_*CL*(G_CODE*mi/rg**2)*Ti**4)
    fr = nab/nrad
    excess = fr > 1.0
    fm = np.clip(1.0 - fr, 0.0, 1.0)
    Fr = F*np.minimum(fr, 1.0)
    eint = rhoi*10.0**st.eos.le(np.log10(Ti), np.log10(rhoi))
    up = col['r'] >= rs
    out = dict(col)
    for k_, v_ in (('r', rg), ('rho', rhoi), ('T', Ti), ('eint', eint), ('Fr', Fr),
                   ('E', A_RAD_*Ti**4), ('fmlt', fm)):
        out[k_] = np.concatenate([v_, col[k_][up]])
    # continuity residuals at r*
    Ps = np.exp(np.interp(rs, rr, lnPo))
    Tin = np.interp(rs, col['r'], col['T'])
    rin = np.exp(np.interp(rs, col['r'], np.log(col['rho'])))
    Tout = float(np.exp(fT(np.log(Ps))))
    guess[0] = rin
    rout = rho_of(Ps, Tout)
    Pin = st.thermo(Tin, rin)[0]
    mw = float(np.interp(c['r_wall'], rg, mi))
    st.RSTAR, st.GM = Rbest, G_CODE*m_e
    out.update(scheme='mesa+two_sided', rho_floor=c['rho_floor'],
               m_wall=mw, m_top=m_e, GM_wall=G_CODE*mw,
               excess=dict(n=int(excess.sum()),
                           max=float((fr[excess] - 1).max()) if excess.any() else 0.0,
                           rlo=float(rg[excess].min()) if excess.any() else np.nan,
                           rhi=float(rg[excess].max()) if excess.any() else np.nan),
               interior=dict(r=rg, m=mi, nab=nab, nrad=nrad, frad=fr), ext=None,
               fit=dict(R=Rbest, how=how, r_star=rs, r_e=r_e, dlnP=float(np.log(Pin/Ps)),
                        dlnT=float(np.log(Tin/Tout)), dlnrho=float(np.log(rin/rout))))
    log('  mode mesa TWO-SIDED FIT (%s): R_ph %.8e (%.5f x target %.6e); fit point r* %.6e '
        '(%.4f Rsun; MESA top reached at %.4f Rsun); continuity residuals dlnP %+.2e dlnT %+.2e '
        'dln rho %+.2e' % (how, Rbest, Rbest/Rt, Rt, rs, rs/6.957e10, r_e/6.957e10,
                           out['fit']['dlnP'], out['fit']['dlnT'], out['fit']['dlnrho']))
    return out
