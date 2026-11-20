# M35MMS (= T35L5.2, Schultz, Bildsten & Jiang 2022 ApJL 924 L11) IC for AthenaK he_star_m1.
# Construction (Jiang+2018: "take the typical luminosity, density and gravity ... in our 1D models and construct
# envelopes in hydrostatic and thermal equilibrium"; S22 sect 2.1: MESA model -> Athena++ IC): an envelope in
# hydrostatic AND radiative equilibrium, all of L carried by radiative diffusion ("radiative transfer only", no MLT),
# with the SIMULATION EOS (ideal gas mu, gamma 5/3, + a T^4/3) and the simulation Rosseland table, MATCHED TO THE
# MESA SURFACE: started at the outermost MESA row (r_ph 13.737 Rsun, tau 0.67, T, rho of MESA) and integrated inward
# (stable) to below 9.7 Rsun, and outward (thin, T ~ const) until the density floor.
# Anchoring at the MESA Fe peak instead is impossible: kappa_R there is 1.09 kappa_Edd, so the outward radiative-only
# solution from MESA's (rho, T) inverts the density and T collapses to 0 at ~13.3 Rsun (tested 11-15).
# MESA source: Will Schultz 35M_save_profile_MS.data (35 Msun, Xc 0.34, X_surf 0.70 Z 0.02).
# Outputs (cgs): ic/ic_m35mms.txt   he_ic_file, he_ic_cols 5 (+T, fmlt=0): r rho eint F_r E | T fmlt
#                ic/prof_m35mms.txt he_ic_linhao_file (inner ghosts): r T rho M kappa F_r
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
OUT = '/viper/ptmp2/jinma/m35mms_1115/ic/'
F = '/viper/ptmp2/jinma/schultz_ics_1113/x/initial_condition_help/IC_from_MESA/35M_save_profile_MS.data'
TAB = '/viper/ptmp2/jinma/run12_viper_1022/in/rosseland_mesakap_gs98_x0.70_z0.02.txt'
Rsun = 6.955e10                      # deck radii (Will's script Rsun); MESA rsun 6.957e10 used to read MESA r
RsunM = 6.957e10; MsunM = 1.98841e33; GM_G = 6.67430e-8   # MESA r21.12.1 header constants
mu = 0.6215                          # [IC] MESA mu column in 9.7-13.3 Rsun (0.6212-0.6225; fully ionised X0.70 Z0.02)
mu_u = 1.66053906660e-24; kB = 1.380649e-16; a = 7.5657332e-15; c = 2.99792458e10; gam = 5/3
GM = GM_G*35.0*MsunM                 # [PAPER] M_core 35 Msun (S23 Table 1), point mass; MESA star_mass 35.0
RLO, RHI = 9.0, 15.6                 # IC file range (Rsun): paper domain 9.7-15.3 + ghosts + coarse-smoke margins
RHOFLOOR = 1.0e-17                   # [NS] IC density floor above the photosphere
names = open(F).read().splitlines()[5].split()
d = np.genfromtxt(F, skip_header=6, names=names)
rM = 10**d['logR']*RsunM; TM = 10**d['logT']; dM = 10**d['logRho']; LM = d['L'][0]; kM = d['kappa']
ia = 0   # outermost MESA row (photosphere, tau 0.67)
ra, Ta, rhoa = rM[ia], TM[ia], dM[ia]
A = open(TAB).read().splitlines(); gi = [i for i, l in enumerate(A) if l.startswith('# grid')][0]
h = A[gi+1].lstrip('#').split(); nT, nD = int(h[0]), int(h[1]); lT0, dlT, lD0, dlD = map(float, h[2:6])
v = np.array([float(x) for x in [l for l in A if l.strip() and not l.startswith('#')][:nT*nD]]).reshape(nT, nD)
fk = RGI((lT0+dlT*np.arange(nT), lD0+dlD*np.arange(nD)), v)   # bilinear in log-log


def kap(T, rho):
    lT = np.clip(np.log10(np.atleast_1d(T)), lT0, lT0+dlT*(nT-1))
    lD = np.clip(np.log10(np.atleast_1d(rho)), lD0, lD0+dlD*(nD-1))
    return 10**fk(np.c_[lT, lD])


def rhs(r, y):   # y = (P_gas, T); dPtot/dr = -GM rho/r^2, dPr/dr = -kappa rho L/(4 pi r^2 c)
    pg, t = y; rh = max(pg, 1e-30)*mu*mu_u/(kB*t); k = kap(t, rh)[0]
    dpr = -k*rh*LM/(4*np.pi*r**2*c); dpt = -GM*rh/r**2
    return np.array([dpt-dpr, dpr/(4*a*t**3/3)])


def integrate(r0, y0, r1, n):
    rr = np.linspace(r0, r1, n+1); hh = rr[1]-rr[0]; R = [r0]; Y = [np.array(y0)]; y = np.array(y0)
    for i in range(n):
        r = rr[i]; k1 = rhs(r, y); k2 = rhs(r+hh/2, y+hh/2*k1); k3 = rhs(r+hh/2, y+hh/2*k2); k4 = rhs(r+hh, y+hh*k3)
        y = y+hh/6*(k1+2*k2+2*k3+k4)
        if y[0]*mu*mu_u/(kB*y[1]) < RHOFLOOR:
            break
        R.append(rr[i+1]); Y.append(y)
    Y = np.array(Y); return np.array(R), Y[:, 0], Y[:, 1]


pga = rhoa*kB*Ta/(mu*mu_u)
r1, p1, t1 = integrate(ra, [pga, Ta], RLO*Rsun, 100000)
r2, p2, t2 = integrate(ra, [pga, Ta], RHI*Rsun, 100000)
r = np.r_[r1[::-1], r2[1:]]; pg = np.r_[p1[::-1], p2[1:]]; T = np.r_[t1[::-1], t2[1:]]
rho = pg*mu*mu_u/(kB*T)
rfl = np.nan
if r[-1] < RHI*Rsun*0.9999:   # floor region: rho = floor, T held at its last value
    rfl = r[-1]/Rsun
    rf = np.linspace(r[-1], RHI*Rsun, 2001)[1:]
    r = np.r_[r, rf]; rho = np.r_[rho, np.full(rf.size, RHOFLOOR)]; T = np.r_[T, np.full(rf.size, T[-1])]
rg = np.exp(np.linspace(np.log(r[0]), np.log(r[-1]), 40001))
T = np.exp(np.interp(rg, r, np.log(T))); rho = np.exp(np.interp(rg, r, np.log(rho))); r = rg
k = kap(T, rho); Fr = LM/(4*np.pi*r**2); E = a*T**4; eint = rho*kB*T/(mu*mu_u)/(gam-1)
tau = np.r_[np.cumsum(((k*rho)[::-1][:-1])*(-np.diff(r[::-1])))[::-1], 0.0]
dm = 4*np.pi*r**2*rho; Mext = np.r_[np.cumsum((0.5*(dm[1:]+dm[:-1])*np.diff(r))[::-1])[::-1], 0.0]
M = GM/GM_G-Mext
sel = (r > 8.5*Rsun) & (r < 15*Rsun); jfe = np.argmax(np.where(sel, k, 0)); jph = np.argmin(abs(tau-1))
print('anchor (MESA top row): r %.4f Rsun T %.4e rho %.4e kappa_MESA %.4f kappa_tab %.4f; L %.5e (log %.4f, Lsun 3.828e33);'
      ' GM %.6e; kappa_Edd %.4f' % (ra/Rsun, Ta, rhoa, kM[ia], kap(Ta, rhoa)[0], LM, np.log10(LM/3.828e33), GM,
                                    4*np.pi*GM*c/LM))
print('IC: Fe peak (max kappa) r %.4f T %.4e kappa %.4f; tau=1 at r %.4f T %.4e; floor from r %.4f'
      % (r[jfe]/Rsun, T[jfe], k[jfe], r[jph]/Rsun, T[jph], rfl))
for rr in [9.7, 10.9, 11.1, 12.0, 12.7, 13.14, 13.5, 13.73, 14.0, 15.3]:
    j = np.argmin(abs(r-rr*Rsun)); jm = np.argmin(abs(rM-rr*Rsun))
    print('r %6.2f: IC T %.4e rho %.4e tau %.4g kap %.3f | MESA T %.4e rho %.4e tau %.4g kap %.3f | dlnT %+.3f dlnrho %+.3f'
          % (rr, T[j], rho[j], tau[j], k[j], TM[jm], dM[jm], d['tau'][jm], kM[jm], np.log(T[j]/TM[jm]), np.log(rho[j]/dM[jm])))
j0 = np.argmin(abs(r-9.7*Rsun)); j1 = np.argmin(abs(r-15.3*Rsun))
print('M_env on 9.7-15.3 Rsun = %.4e Msun; MESA m_ext(9.7) = %.4e' % ((Mext[j0]-Mext[j1])/MsunM,
      d['mass'][0]-d['mass'][np.argmin(abs(rM-9.7*RsunM))]))
dr = np.gradient(rho, r); inv = (dr > 0) & (r > 9.7*Rsun) & (rho > RHOFLOOR*1.01)
if inv.any():
    print('density inversion between %.4f and %.4f Rsun' % (r[inv][0]/Rsun, r[inv][-1]/Rsun))
# opacity check along the MESA profile: our table vs MESA kappa
sm = (rM > 9.7*RsunM) & (rM < 13.73*RsunM)
rat = kap(TM[sm], dM[sm])/kM[sm]
print('table/MESA kappa along the MESA profile 9.7-13.73: min %.4f max %.4f' % (rat.min(), rat.max()))
hdr = ('# he_star_m1 IC (he_ic_cols 5, 7 columns): M35MMS = T35L5.2 (Schultz+2022 ApJL 924 L11; Schultz+2023 Table 1)\n'
       '# radiative-only HSE + thermal equilibrium (Jiang+2018 construction) matched to the MESA surface row'
       ' r %.5f Rsun, T %.5e, rho %.5e\n'
       '# of Will Schultz 35M_save_profile_MS.data; L %.6e erg/s (MESA), GM %.6e (35 Msun, MESA constants), mu %.4f,'
       ' gamma 5/3, a %.7e\n'
       '# kappa_R mesakap gs98 X0.70 Z0.02; rho floor %.1e above the photosphere; F_r = L/(4 pi r^2); E = a T^4\n'
       '# r[cm] rho eint F_r E | T[K] fmlt') % (ra/Rsun, Ta, rhoa, LM, GM, mu, a, RHOFLOOR)
np.savetxt(OUT+'ic_m35mms.txt', np.c_[r, rho, eint, Fr, E, T, np.zeros_like(r)], fmt='%.12e', header=hdr, comments='')
hdr2 = ('# he_bc_linhao profile (r[cm] T[K] rho[g/cc] M[g] kappa F_r) for the inner ghosts: T, rho = the IC column (fixed),'
        ' F_r = L/(4 pi r^2),\n# L = %.6e (MESA L of 35M_save_profile_MS.data; S22 sect 2.1 fixed T, rho, F_r ghosts)\n'
        '# r T rho M kappa F_r') % LM
np.savetxt(OUT+'prof_m35mms.txt', np.c_[r, T, rho, M*MsunM, k, Fr], fmt='%.12e', header=hdr2, comments='')
print('F_in at 9.7 Rsun (Rsun 6.955e10) = %.10e' % (LM/(4*np.pi*(9.7*Rsun)**2)))
print('GM_for_deck %.10e' % GM)
