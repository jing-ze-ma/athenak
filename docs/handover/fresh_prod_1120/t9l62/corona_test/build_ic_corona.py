# build_ic_corona.py: T9L6.2 test IC (corona test, 11-20) = ic_t9l62_LIC.txt with the floor gas above the
# photosphere replaced by a static, hydrostatic hot corona.
# Below the join radius r_j the IC is unchanged (Y.-F. Jiang init_outburst.txt + the ghost extension). Above it:
#   rho_c(r) = RHO_TOP (r/r_top)^-2, r_top = deck x1max, RHO_TOP = 2 x the test dfloor (1e-16), so the corona stays
#   above the floor everywhere and the code floor never acts on it;
#   P_c from dP/dr = -rho_c GM/r^2 (1 - kappa_R(T, rho_c) L_IC/(4 pi c GM)), integrated inward from r_top with
#   T_top from the isothermal-free power-law solution kT/(mu m_u) = GM (1 - Gamma)/(3 r) (rho ~ r^-2, T ~ 1/r).
# r_j = the first radius above the photosphere where the IC gas pressure falls to P_c: P is continuous there, rho and
# T jump (a static contact). The radiation columns (F_r = L_IC/(4 pi r^2), E) are kept from the original IC, so the
# only change is the gas above r_j. eint is recomputed from T (the deck reads he_ic_eint_from_t = true anyway).
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
SRC = '/viper/ptmp2/jinma/t9l62_1113/ic/ic_t9l62_LIC.txt'
OUT = '/viper/ptmp2/jinma/t9l62_1113/ic/ic_t9l62_corona.txt'
TAB = '/viper/ptmp2/jinma/run12_viper_1022/in/rosseland_mesakap_gs98_x0.70_z0.02.txt'
mu = 0.674; mu_u = 1.66053906660e-24; kB = 1.380649e-16; c = 2.99792458e10; gam = 5/3
GM = 7.4033817458e27; LIC = 2.47875e39; Rsun = 6.955e10
R_TOP = 5.632159e13                 # deck x1max (809.8 Rsun)
DFLOOR = 1.0e-16                     # test deck dfloor
RHO_TOP = 2.0*DFLOOR
A = open(TAB).read().splitlines(); gi = [i for i, l in enumerate(A) if l.startswith('# grid')][0]
h = A[gi+1].lstrip('#').split(); nT, nD = int(h[0]), int(h[1]); lT0, dlT, lD0, dlD = map(float, h[2:6])
v = np.array([float(x) for x in [l for l in A if l.strip() and not l.startswith('#')][:nT*nD]]).reshape(nT, nD)
fk = RGI((lT0+dlT*np.arange(nT), lD0+dlD*np.arange(nD)), v)


def kap(T, rho):
    lT = np.clip(np.log10(np.atleast_1d(T)), lT0, lT0+dlT*(nT-1))
    lD = np.clip(np.log10(np.atleast_1d(rho)), lD0, lD0+dlD*(nD-1))
    return 10**fk(np.c_[lT, lD])


hdr0 = [l for l in open(SRC) if l.startswith('#')]
d = np.loadtxt(SRC); r, rho, eint, Fr, E, T, fm = d.T
rho_c = lambda x: RHO_TOP*(x/R_TOP)**-2
# corona pressure: RK4 inward on a fine log grid from beyond the last IC row
rmax = r.max()*1.0001
gam_top = kap(1e7, rho_c(rmax))[0]*LIC/(4*np.pi*c*GM)
P = rho_c(rmax)*GM*(1-gam_top)/(3*rmax)
lg = np.linspace(np.log(rmax), np.log(60*Rsun), 200001)
rr = np.exp(lg); PP = np.empty_like(rr); PP[0] = P


def dPdr(x, p):
    rc = rho_c(x); t = p*mu*mu_u/(kB*rc)
    return -rc*GM/x**2*(1-kap(t, rc)[0]*LIC/(4*np.pi*c*GM))


for n in range(1, rr.size):
    x0, x1 = rr[n-1], rr[n]; hh = x1-x0; p = PP[n-1]
    k1 = dPdr(x0, p); k2 = dPdr(x0+hh/2, p+hh/2*k1); k3 = dPdr(x0+hh/2, p+hh/2*k2); k4 = dPdr(x1, p+hh*k3)
    PP[n] = p+hh/6*(k1+2*k2+2*k3+k4)
Pc = np.exp(np.interp(np.log(r), lg[::-1], np.log(PP[::-1])))
Pg = rho*kB*T/(mu*mu_u)
iph = np.argmin(abs(r-93.9*Rsun))   # tau = 1 of the IC (T9L62_SETUP sect 2)
cand = np.nonzero((np.arange(r.size) > iph) & (Pg <= Pc))[0]
j = cand[0]
Tc = Pc*mu*mu_u/(kB*rho_c(r))
new_rho = np.where(np.arange(r.size) >= j, rho_c(r), rho)
new_T = np.where(np.arange(r.size) >= j, Tc, T)
new_e = new_rho*kB*new_T/(mu*mu_u)/(gam-1)
G_c = kap(new_T, new_rho)*LIC/(4*np.pi*c*GM)
sel = np.arange(r.size) >= j
print('join r_j %.4f Rsun: IC rho %.3e T %.3e P %.3e -> corona rho %.3e T %.3e (density jump %.1f)'
      % (r[j]/Rsun, rho[j], T[j], Pg[j], new_rho[j], new_T[j], rho[j]/new_rho[j]))
print('corona: T %.3e K at r_j, %.3e K at x1max; rho %.3e .. %.3e; Gamma_R %.3f .. %.3f; min rho/dfloor %.2f'
      % (new_T[j], np.interp(R_TOP, r, new_T), new_rho[j], np.interp(R_TOP, r, new_rho), G_c[sel].min(),
         G_c[sel].max(), new_rho[sel & (r <= R_TOP*1.02)].min()/DFLOOR))
M_c = np.trapz(4*np.pi*r[sel & (r <= R_TOP)]**2*new_rho[sel & (r <= R_TOP)], r[sel & (r <= R_TOP)])
print('corona mass (r_j..x1max, full sphere) %.3e g = %.3e Msun' % (M_c, M_c/1.989e33))
# HSE check of the corona as written
gP = np.gradient(new_rho*kB*new_T/(mu*mu_u), r)
res = (gP+new_rho*GM/r**2*(1-G_c))/(new_rho*GM/r**2)
s2 = sel & (r < R_TOP) & (np.arange(r.size) > j+2)
print('corona HSE residual (fraction of rho g): median %.1e max %.1e' % (np.median(abs(res[s2])), abs(res[s2]).max()))
cs = np.sqrt(gam*new_rho*kB*new_T/(mu*mu_u)/new_rho)
print('corona sound speed %.2e cm/s at r_j, %.2e at x1max' % (cs[j], np.interp(R_TOP, r, cs)))
hdr = (''.join(hdr0[:-1]).rstrip('\n') + '\n'
       '# CORONA TEST 11-20 (build_ic_corona.py): above r_j = %.6e cm (%.4f Rsun) the floor gas is replaced by a static'
       ' hydrostatic corona rho = %.1e (r/%.6e)^-2, P from HSE with gravity GM (1 - kappa_R L_IC/(4 pi c GM));'
       ' F_r and E columns unchanged; T column = gas T\n' % (r[j], r[j]/Rsun, RHO_TOP, R_TOP) + hdr0[-1].rstrip('\n'))
np.savetxt(OUT, np.c_[r, new_rho, new_e, Fr, E, new_T, fm], fmt='%.12e', header=hdr, comments='')
print('wrote', OUT)
