"""COPY of grains/grainlib.py for the Betelgeuse observed-density analysis (betelgeuse/).
Changes: Betelgeuse stars (FT eq. 7 fitted to the observed profile), dens='obs' parcel density,
T_gas may be a function of r, Td table to 100 R*, run_steady() along a steady-wind v(r).
Original docstring: Grain growth of Mg2SiO4 in Fuller & Tsuna (2024) chromospheric parcels.
Units cgs. See report for every parameter and its source."""
import numpy as np
import prof as P
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator

D = '/orion/ptmp/jinma/rsg_wind_1008/grains/out/'
G, c, kB, mu_, Msun, Lsun, Rsun = 6.674e-8, 2.99792458e10, 1.380649e-16, 1.66054e-24, \
    1.989e33, 3.828e33, 6.957e10
YR = 3.156e7

# ---------- material (forsterite) ----------
RHO_D = 3.21                         # g/cm^3 (Hoefner 2008; Gail & Sedlmayr 1999: 3.21)
M_FO = 140.69 * mu_                  # g per Mg2SiO4 formula unit
V_MON = M_FO / RHO_D                 # cm^3
M_MG = 24.305 * mu_
A_SEED = 1e-7                        # 1 nm

STARS = {  # MESA columns, rsg_wind_1008/mesa/col_<m>_mft_clamp.npz (FT-way v_con)
    'golden16': dict(M=15.26, L=1.149e5, R=669.9, Teff=4106.3, rho_ph=1.0357e-9,
                     v_con=8.175e5, npz='m15lgl5.1'),
    'm20lgl5.5': dict(M=19.37, L=3.15e5, R=1312., Teff=3776.4, rho_ph=6.5316e-10,
                      v_con=8.298e5, npz='m20lgl5.5'),
}


def bet_star(sset, M, prof='dent'):
    st = P.stellar(sset, M)
    f = P.ft_fit(st, prof)
    return dict(M=M, L=st['L']/Lsun, R=st['R']/Rsun, Teff=st['Teff'], rho_ph=f['rho_ph'],
                v_con=f['vcon'], d=st['d'], prof=prof)


for _s, _m in (('A', 18.), ('A', 20.), ('B', 16.5), ('B', 19.)):
    STARS[f'bet{_s}{_m:g}'] = bet_star(_s, _m)


class Star:
    def __init__(self, name):
        s = STARS[name]
        self.d, self.prof = s.get('d'), s.get('prof')
        self.name = name
        self.M, self.L, self.R = s['M'] * Msun, s['L'] * Lsun, s['R'] * Rsun
        self.Teff, self.rho_ph, self.v_con = s['Teff'], s['rho_ph'], s['v_con']
        self.vesc = np.sqrt(2 * G * self.M / self.R)
        self.kE = 4 * np.pi * c * G * self.M / self.L          # cm^2/g gas
        self.mdot0 = 4 * np.pi * self.R**2 * self.rho_ph * self.v_con

    def W(self, r):
        x = np.minimum(self.R / r, 1.0)
        return 0.5 * (1 - np.sqrt(1 - x**2))

    def rho_ft(self, r):
        """FT24 eq. 7 time-averaged chromosphere (continued beyond R_d = 5 R)."""
        x = self.R / np.maximum(r, self.R)
        return self.rho_ph * x**2 * np.exp(-(self.vesc / self.v_con) * np.sqrt(1 - x))

    def rho_obs(self, r):
        return float(P.rho(np.array([r / self.R]), self.prof, self.d)[0])

    def T_obs(self, r):
        return float(P.Tobs(np.array([r / self.R]), self.prof)[0])

    def mdot_gt(self, v0):
        """Mass flux launched faster than v0: FT24 eq. 15 form, 4 pi R^2 rho_ph v_con
        exp(-v0/v_con) (their exponential velocity distribution; eq. 7 <-> eq. 15)."""
        return self.mdot0 * np.exp(-v0 / self.v_con)


class Optics:
    def __init__(self, var):
        o = np.load(D + 'optics.npz')
        iv = list(o['variants']).index(var)
        self.var = var
        self.la = np.log(o['a'])
        self.lT = np.log(o['T'])
        self.QaT = o['QaT'][iv]            # [a, T]
        self.QprT = o['QprT'][iv]
        self.o = o
        self.iv = iv

    def Qpr_star(self, a, Teff):
        q = np.array([np.interp(np.log(Teff), self.lT, self.QprT[i])
                      for i in range(len(self.la))])
        return np.interp(np.log(a), self.la, q)

    def Td_table(self, Teff, W):
        """radiative equilibrium: W <Qabs>_Teff Teff^4 = <Qabs>_Td Td^4 (no gas heating)."""
        Td = np.zeros((len(self.la), len(W)))
        for i in range(len(self.la)):
            qs = np.interp(np.log(Teff), self.lT, self.QaT[i])
            f = np.log(self.QaT[i]) + 4 * self.lT                 # ln(Q T^4), monotonic
            Td[i] = np.exp(np.interp(np.log(W * qs) + 4 * np.log(Teff), f, self.lT))
        return Td


class Chem:
    """ln a0(T,P): forsterite activity in an equilibrium gas (FastChem, no condensation)."""
    def __init__(self):
        ch = np.load(D + 'chem.npz')
        self.ch = ch
        self.f_lna = RegularGridInterpolator((ch['T'], ch['lP']), ch['lna'],
                                             bounds_error=False, fill_value=None)
        self.f_nt = RegularGridInterpolator((ch['T'], ch['lP']), ch['ntot'] / ch['nH'],
                                            bounds_error=False, fill_value=None)
        self.f_h2o = RegularGridInterpolator((ch['T'], ch['lP']), ch['x_H2O1'],
                                             bounds_error=False, fill_value=None)
        self.eMg, self.eSi = float(ch['eps_Mg']), float(ch['eps_Si'])
        # gas mass per H (Asplund 2009: He 0.0851, metals Z = 0.0134 by mass)
        self.mH = 1.008 * mu_ / 0.7381          # X = 0.7381 (Asplund 2009)

    def lnS(self, Td, nH, fcond):
        """ln of the Mg (key species) supersaturation for a grain at Td, using the
        vapour in equilibrium at Td with the same element densities (Mg depleted by fcond,
        SiO by 0.614 fcond, H2O by 3 eMg/2 fcond / x_H2O). S_key = a^(1/2) (nu_Mg = 2)."""
        P = nH * 0.586 * kB * Td / 1e6
        for _ in range(2):
            lP = np.log10(np.clip(P, 1e-16, 1e-2))
            P = nH * self.f_nt((Td, lP)) * kB * Td / 1e6
        lP = np.log10(np.clip(P, 1e-16, 1e-2))
        la = self.f_lna((np.clip(Td, 300, 3000), lP))
        xh2o = self.f_h2o((np.clip(Td, 300, 3000), lP))
        dep = (2 * np.log(max(1 - fcond, 1e-12)) + np.log(max(1 - 0.5 * self.eMg / self.eSi
                                                                * fcond, 1e-12))
               + 3 * np.log(max(1 - 1.5 * self.eMg * fcond / max(xh2o, 1e-30), 1e-12)))
        return 0.5 * (la + dep)


class Model:
    def __init__(self, star, var, chem=None):
        self.st = Star(star)
        self.op = Optics(var)
        self.ch = chem or Chem()
        self.lW = np.linspace(np.log(self.st.W(100 * self.st.R)), np.log(0.5), 300)
        self.Td = self.op.Td_table(self.st.Teff, np.exp(self.lW))    # [a, W]
        self.fTd = RegularGridInterpolator((self.op.la, self.lW), np.log(self.Td),
                                           bounds_error=False, fill_value=None)
        aa = np.exp(self.op.la)
        self.qpr = self.op.Qpr_star(aa, self.st.Teff)

    def Td_at(self, a, r):
        return float(np.exp(self.fTd((np.log(max(a, 1e-7)), np.log(self.st.W(r))))))

    def kpr_dust(self, a):
        """flux-mean radiation-pressure opacity per gram of dust, cm^2/g"""
        return 3 * np.interp(np.log(a), self.op.la, self.qpr) / (4 * a * RHO_D)

    def gamma_d(self, a, xseed):
        """Gamma of the dust per gram of GAS: xseed pi a^2 Qpr / (m_gas per H) / kE"""
        q = np.interp(np.log(a), self.op.la, self.qpr)
        return xseed * np.pi * a**2 * q / self.ch.mH / self.st.kE

    def fcond(self, a, xseed):
        return xseed * (4 / 3 * np.pi * (a**3 - A_SEED**3) * RHO_D / M_FO) \
            / (0.5 * self.ch.eMg)

    def dadt(self, a, r, nH, Tg, xseed, alpha):
        """da/dt = V_mon * alpha * (n_Mg vbar/4)/2 * (1 - 1/S_eff),
        1/S_eff = sqrt(Td/Tg)/S_key (Gail & Sedlmayr 1999 key-species form)."""
        f = min(self.fcond(a, xseed), 0.999999)
        Td = self.Td_at(a, r)
        lnS = self.ch.lnS(Td, nH, f)
        vbar = np.sqrt(8 * kB * Tg / (np.pi * M_MG))
        J = alpha * nH * self.ch.eMg * (1 - f) * vbar / 4 / 2
        return V_MON * J * (1 - np.sqrt(Td / Tg) * np.exp(-lnS)), Td, lnS

    def run(self, v0, alpha, xseed, Tg, gam_gas, dens='ft', tmax=3e10, rout=20.):
        st = self.st
        R = st.R
        mdot = st.mdot_gt(v0)

        def nH_of(r, v):
            if dens == 'ft':
                rho = st.rho_ft(r)
            elif dens == 'obs':
                rho = st.rho_obs(r)
            else:   # mass-conserving stream of all gas launched faster than v0
                rho = mdot / (4 * np.pi * r**2 * max(abs(v), st.v_con))
            return rho / self.ch.mH

        Tgf = Tg if callable(Tg) else (lambda r_: Tg)

        def rhs(t, y):
            r, v, a = y
            nH = nH_of(r, v)
            Tg_ = Tgf(r)
            Td = self.Td_at(a, r)
            # seeds only once the seed grain could survive (S(seed) > 1); else a stays
            if a <= A_SEED * 1.0001:
                lnS0 = self.ch.lnS(self.Td_at(A_SEED, r), nH, 0.0)
                if lnS0 <= 0:
                    da = 0.0
                else:
                    da = self.dadt(A_SEED, r, nH, Tg_, xseed, alpha)[0]
            else:
                da = self.dadt(a, r, nH, Tg_, xseed, alpha)[0]
            if a <= A_SEED and da < 0:
                da = 0.0
            gd = self.gamma_d(a, xseed) if a > A_SEED * 1.0001 else 0.0
            return [v, -G * st.M / r**2 * (1 - gam_gas - gd), da]

        def ev_fall(t, y):
            return y[0] - R * 0.999
        ev_fall.terminal, ev_fall.direction = True, -1

        def ev_out(t, y):
            return y[0] - rout * R
        ev_out.terminal = True
        sol = solve_ivp(rhs, (0, tmax), [R, v0, A_SEED], method='LSODA', rtol=1e-6,
                        atol=[1e8, 1e2, 1e-10], events=[ev_fall, ev_out], max_step=2e7,
                        dense_output=False)
        t, (r, v, a) = sol.t, sol.y
        a = np.maximum(a, A_SEED)
        if sol.status == 1 and len(sol.t_events[1]):
            out = 'wind'
        elif sol.status == 1:
            out = 'fall'
        else:
            out = 'bound(tmax)'
        gd = np.array([self.gamma_d(x, xseed) if x > A_SEED * 1.0001 else 0 for x in a])
        res = dict(t=t, r=r / R, v=v, a=a, gd=gd, outcome=out,
                   f=np.array([self.fcond(x, xseed) for x in a]))
        res['Td'] = np.array([self.Td_at(x, y) for x, y in zip(a, r)])
        res['nH'] = np.array([nH_of(x, y) for x, y in zip(r, v)])
        res['rmax'] = r.max() / R
        m3 = r > 3 * R
        dt = np.diff(t, prepend=t[0])
        res['t3'] = dt[m3].sum()
        res['amax'] = a.max()
        res['fmax'] = res['f'].max()
        res['gdmax'] = gd.max()
        if out == 'wind':
            ve = np.sqrt(2 * G * st.M / r[-1]) * np.sqrt(max(1 - gam_gas - gd[-1], 0))
            res['vinf'] = np.sqrt(max(v[-1]**2 - ve**2, 0))
            res['mdot'] = mdot
        # drift (Kwok 1975 interpolation): F_rad = pi a^2 Qpr L/(4 pi r^2 c) =
        # pi a^2 rho vd sqrt(vd^2 + (4/3 vbar_gas)^2)   [gas vbar with mu = 2.34 at Tg]
        q = np.interp(np.log(a), self.op.la, self.qpr)
        Frad = q * st.L / (4 * np.pi * r**2 * c)
        rho = res['nH'] * self.ch.mH
        Tga = np.array([Tgf(x) for x in r])
        vth = 4 / 3 * np.sqrt(8 * kB * Tga / (np.pi * 2.34 * mu_))
        res['vd'] = np.sqrt(0.5 * (-vth**2 + np.sqrt(vth**4 + 4 * (Frad / rho)**2)))
        res['tstop'] = 4 * a * RHO_D / (3 * rho * np.sqrt(res['vd']**2 + vth**2))
        return res


    def run_steady(self, vfun, alpha, xseed, Tg, x0=1.2, x1=40.):
        """Grain growth along a prescribed steady-wind trajectory dr/dt = v(r) through the
        observed density (kinematic: no feedback of the dust on v). vfun(r) [cm/s]."""
        st = self.st
        R = st.R
        Tgf = Tg if callable(Tg) else (lambda r_: Tg)

        def rhs(r, y):
            a, t = y
            nH = st.rho_obs(r) / self.ch.mH
            v = vfun(r)
            if a <= A_SEED * 1.0001:
                lnS0 = self.ch.lnS(self.Td_at(A_SEED, r), nH, 0.0)
                da = self.dadt(A_SEED, r, nH, Tgf(r), xseed, alpha)[0] if lnS0 > 0 else 0.
            else:
                da = self.dadt(a, r, nH, Tgf(r), xseed, alpha)[0]
            if a <= A_SEED and da < 0:
                da = 0.0
            return [da / v, 1. / v]
        rr = R * np.exp(np.linspace(np.log(x0), np.log(x1), 400))
        sol = solve_ivp(rhs, (rr[0], rr[-1]), [A_SEED, 0.], method='LSODA', t_eval=rr,
                        rtol=1e-6, atol=[1e-11, 1e3], max_step=0.02 * R)
        a = np.maximum(sol.y[0], A_SEED)
        gd = np.array([self.gamma_d(x, xseed) if x > A_SEED * 1.0001 else 0 for x in a])
        Td = np.array([self.Td_at(x, y) for x, y in zip(a, sol.t)])
        return dict(r=sol.t / R, a=a, gd=gd, Td=Td, t=sol.y[1], ok=sol.success,
                    f=np.array([self.fcond(x, xseed) for x in a]))


def rmax_dustfree(st, v0, gam):
    x = v0**2 / (st.vesc**2 * (1 - gam))
    return 1 / (1 - x) if x < 1 else np.inf
