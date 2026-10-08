"""Required radiative Eddington factor for a steady spherical wind through the observed
Betelgeuse density profile.

  v = Mdot/(4 pi r^2 rho_obs);  beyond the radius where v reaches v_inf the wind coasts at
  v_inf and rho = Mdot/(4 pi r^2 v_inf).
  v dv/dr = -GM/r^2 (1 - Gamma) - (1/rho) d(P_gas + rho v_t^2)/dr
  => Gamma_req = 1 + r^2/(GM) [v dv/dr + (1/rho) dP_gas/dr + v_t^2 dln rho/dr]
P_gas = rho k T/(mu m_H), mu from the FastChem CE table (ck_lib) at (T_obs, P).
"""
import os
import numpy as np
os.environ.setdefault('CK_DATA', '/orion/ptmp/jinma/rsg_wind_1008/dace/data/')
import ck_lib as ck  # noqa: E402
import prof as P  # noqa: E402

CE = ck.read_ce(ck.DATA + 'CE_tables/FastChem_ck_1x_int_hiT2_lowP.txt')


def mu_of(T, rho):
    mu = np.full_like(T, 1.3)
    for _ in range(30):
        Pg = rho*P.KB*T/(mu*P.MH)
        mn = np.array([ck.comp(CE, t, p*1e-6)[0] for t, p in zip(T, Pg)])
        if np.max(np.abs(mn - mu)/mu) < 1e-7:
            break
        mu = mn
    return mn


def xgrid(x0=1.2, x1=30., n=600):
    return np.exp(np.linspace(np.log(x0), np.log(x1), n))


def steady(st, prof='dent', mdot=2e-6, vinf=12e5, x=None, vts=(0., 5e5, 10e5),
           Tfix=None, mu_fix=None):
    x = xgrid() if x is None else x
    r = x*st['R']
    Md = mdot*P.MSUN/P.YR
    rho = P.rho(x, prof, st['d'])
    v = Md/(4*np.pi*r**2*rho)
    ic = np.where(v >= vinf)[0]
    xc = np.nan
    if len(ic):
        i = ic[0]
        xc = x[i]
        v = np.where(np.arange(len(x)) >= i, vinf, v)
        rho = Md/(4*np.pi*r**2*v)
    T = P.Tobs(x, prof) if Tfix is None else np.full_like(x, Tfix)
    mu = mu_of(T, rho) if mu_fix is None else np.full_like(x, mu_fix)
    Pg = rho*P.KB*T/(mu*P.MH)
    dvdr = np.gradient(v, r)
    dPdr = np.gradient(Pg, r)
    dlr = np.gradient(np.log(rho), r)
    gr = P.G*st['M']/r**2
    G0 = 1 + (v*dvdr + dPdr/rho)/gr
    out = dict(x=x, r=r, rho=rho, v=v, T=T, mu=mu, P=Pg, dvdr=dvdr, dlnrho_dr=dlr, g=gr,
               x_vinf=xc, mdot=mdot, prof=prof, Hrho=-1/dlr,
               term_acc=v*dvdr/gr, term_gas=dPdr/rho/gr)
    out['Greq'] = {vt: G0 + vt**2*dlr/gr for vt in vts}
    return out


def vt_equiv(w, Gav):
    """Isotropic turbulent/wave velocity whose pressure rho v_t^2 closes Gamma_req(v_t=0)
    - Gamma_avail; nan where no extra support is needed or where rho does not fall."""
    need = w['Greq'][0.] - Gav
    v2 = need*w['g']/(-w['dlnrho_dr'])
    return np.where((need > 0) & (w['dlnrho_dr'] < 0), np.sqrt(np.maximum(v2, 0)), np.nan)


def crossings(x, y, level=1.):
    s = np.sign(y - level)
    i = np.where(s[1:] != s[:-1])[0]
    return [float(np.exp(np.interp(0., [y[k] - level, y[k+1] - level][::(1 if y[k] < y[k+1]
                                                                       else -1)],
                                   np.log([x[k], x[k+1]])[::(1 if y[k] < y[k+1] else -1)])))
            for k in i]


VTS26 = (0., 5e5, 10e5, 14e5, 19.5e5)


def mdot26(st, vinf=9e5):
    """Mdot [Msun/yr] of the SEM wind term (eq. 1): 4 pi R*^2 n_* m_H v_inf / 0.998^1.5."""
    return 4*np.pi*st['R']**2*P.NST*P.NSCALE*P.MPH*vinf/0.998**P.DEL*P.YR/P.MSUN


def steady26(st, x=None, vinf=9e5, vts=VTS26, Tmode=None, mu_fix=None, Tfix=None):
    """Steady wind through the Dent+2026 SEM: v(r) from the wind term (P.v26), rho = total
    SEM rho (chromosphere + wind); v_c = continuity velocity Mdot/(4 pi r^2 rho_total)."""
    x = xgrid(1.15, 30., 600) if x is None else x
    r = x*st['R']
    rho = P.rho(x, 'dent26')
    v = P.v26(x, vinf)
    md = mdot26(st, vinf)
    vc = md*P.MSUN/P.YR/(4*np.pi*r**2*rho)
    T = P.T26(x, Tmode) if Tfix is None else np.where(x >= 1.2, Tfix, P.T26(x, Tmode))
    mu = mu_of(T, rho) if mu_fix is None else np.full_like(x, mu_fix)
    Pg = rho*P.KB*T/(mu*P.MH)
    dvdr = np.gradient(v, r)
    dPdr = np.gradient(Pg, r)
    dlr = np.gradient(np.log(rho), r)
    gr = P.G*st['M']/r**2
    G0 = 1 + (v*dvdr + dPdr/rho)/gr
    out = dict(x=x, r=r, rho=rho, v=v, vc=vc, T=T, mu=mu, P=Pg, dvdr=dvdr, dlnrho_dr=dlr,
               g=gr, mdot=md, prof='dent26', Hrho=-1/dlr, x_vinf=np.nan,
               term_acc=v*dvdr/gr, term_gas=dPdr/rho/gr,
               f_chrom=1 - P.nH26_wind(x)/P.nH26(x))
    out['Greq'] = {vt: G0 + vt**2*dlr/gr for vt in vts}
    return out
