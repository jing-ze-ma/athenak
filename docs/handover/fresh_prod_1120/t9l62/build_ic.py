# T9L6.2 IC for AthenaK he_star_m1 from Yan-Fei Jiang's init_outburst.txt (r/Rsun, T/1.66724e5 K, rho/3.6e-9 g/cc)
# Outputs (cgs):  ic/ic_t9l62_<arm>.txt   he_ic_file, he_ic_cols 5 (+T, fmlt=0): r rho eint F_r E | T fmlt
#                 ic/prof_t9l62_<arm>.txt he_ic_linhao_file format for he_bc_linhao ghosts: r T rho M kappa F_r
# The file below r_file_bottom (24.2286 Rsun) is extended by integrating HSE + radiative diffusion downward
# with the IC's own GM, L_IC, mu and the Rosseland table (only ghost cells of arm A use it).
import numpy as np, sys
OUT='/viper/ptmp2/jinma/t9l62_1113/ic/'
SRC='/viper/ptmp2/jinma/schultz_ics_1113/x/initial_condition_help/IC_from_Yan-Fei/init_outburst.txt'
TAB='/viper/ptmp2/jinma/run12_viper_1022/in/rosseland_mesakap_gs98_x0.70_z0.02.txt'
T0=1.66724e5; rho0=3.6e-9; Rsun=6.955e10          # Will's id_initial_condition.py / Jiang+2018 Fig 2 units
mu=0.674; mu_u=1.66053906660e-24; kB=1.380649e-16  # Will's script mu, Rgas = 8.31447e7/mu (k/m_u)
a=7.5657332e-15; c=2.99792458e10; G=6.67408e-8; Msun=1.989e33; Lsun=3.828e33
gam=5.0/3.0
LBASE=10**6.20*Lsun                                 # Schultz+2020 Table 1 log L = 6.20; Jiang+2018 'fixed luminosity from the bottom boundary'
def tab(fn):
    A=open(fn).read().splitlines(); gi=[i for i,l in enumerate(A) if l.startswith('# grid')][0]
    h=A[gi+1].lstrip('#').split(); nT,nD=int(h[0]),int(h[1]); lT0,dlT,lD0,dlD=map(float,h[2:6])
    L=[l for l in A if l.strip() and not l.startswith('#')]
    v=np.array([float(x) for x in L[:nT*nD]]).reshape(nT,nD)
    from scipy.interpolate import RegularGridInterpolator as R
    f=R((lT0+dlT*np.arange(nT),lD0+dlD*np.arange(nD)),v)
    return lambda T,rho: 10**f(np.c_[np.atleast_1d(np.log10(T)),np.atleast_1d(np.log10(rho))])
kap=tab(TAB)
d=np.loadtxt(SRC); r=d[:,0]*Rsun; T=d[:,1]*T0; rho=d[:,2]*rho0
k=kap(T,rho)
Pg=rho*kB*T/(mu*mu_u); Pr=a*T**4/3
dPg=np.gradient(Pg,r); dPr=np.gradient(Pr,r)
GMr=-(dPg+dPr)/rho*r**2; Lr=-4*np.pi*r**2*c/(k*rho)*dPr
sel=(r>27*Rsun)&(r<92*Rsun)&(np.gradient(rho,r)<0)   # smooth radiative part, skip the bottom rows and the flat floor
GM=np.median(GMr[sel]); LIC=np.median(Lr[sel])
print('GM fit = %.6e (= %.4f G Msun, Msun 1.989e33), L_IC = %.5e (log L/Lsun = %.4f), spread GM %.2e L %.2e'%(GM,GM/G/Msun,LIC,np.log10(LIC/Lsun),np.std(GMr[sel])/GM,np.std(Lr[sel])/LIC))
# downward extension 24.2286 -> 18.0 Rsun (RK4 in r on P_g, T: dPtot/dr = -GM rho/r^2, dPr/dr = -kappa rho L/(4 pi r^2 c))
def rhs(rr,y):
    pg,t=y; rh=pg*mu*mu_u/(kB*t); kk=kap(t,rh)[0]
    dpr=-kk*rh*LIC/(4*np.pi*rr**2*c); dpt=-GM*rh/rr**2
    dt=dpr/(4*a*t**3/3); return np.array([dpt-dpr,dt])
y=np.array([Pg[0],T[0]]); rr=r[0]; h=-(r[1]-r[0])*50
ext=[]
while rr>18.0*Rsun:
    k1=rhs(rr,y);k2=rhs(rr+h/2,y+h/2*k1);k3=rhs(rr+h/2,y+h/2*k2);k4=rhs(rr+h,y+h*k3)
    y=y+h/6*(k1+2*k2+2*k3+k4); rr=rr+h; ext.append((rr,y[1],y[0]*mu*mu_u/(kB*y[1])))
ext=np.array(ext[::-1])
re=np.r_[ext[:,0],r]; Te=np.r_[ext[:,1],T]; de=np.r_[ext[:,2],rho]
print('extension: %d rows down to %.4f Rsun, T %.4e rho %.4e'%(len(ext),re[0]/Rsun,Te[0],de[0]))
ke=kap(Te,de)
# IC radiation flux: the diffusion flux of the IC itself (L_IC) in the optically thick part; above the file's
# tau = 1 radius it is held at L_IC/(4 pi r^2) (constant-L free streaming); F_r = L_IC/(4 pi r^2) everywhere then.
Fr=LBASE/(4*np.pi*re**2)   # IC F_r = the imposed base luminosity (pgen requires F_col(r_in) = implicit_flux_x1min); thick interior relaxes F to diffusion at once
eint=de*kB*Te/(mu*mu_u)/(gam-1.0)
E=a*Te**4
# enclosed envelope mass for the profile's M column (diagnostic only; not read by the BC)
dm=4*np.pi*re**2*de; M=GM/G+np.r_[0,np.cumsum(0.5*(dm[1:]+dm[:-1])*np.diff(re))]
hdr=('# he_star_m1 IC (he_ic_cols 5, 7 columns): T9L6.2 (Jiang+2018 Nature 561 498; Schultz+2020 ApJ 902 67 Table 1) from Y.-F. Jiang init_outburst.txt\n'
     '# units of the source: r/Rsun (6.955e10), T/1.66724e5 K, rho/3.6e-9 g/cc; ideal gas gamma 5/3 mu %.3f (k/m_u); a %.7e\n'
     '# below 24.2286 Rsun: RK4 HSE + diffusion with GM %.6e, L_IC %.5e (log %.4f), kappa_R mesakap gs98 X0.70 Z0.02\n'
     '# F_r = L_BASE/(4 pi r^2), L_BASE = 10^6.20 Lsun (the imposed base luminosity; the T column itself carries L_IC); E = a T^4\n# r[cm] rho eint F_r E | T[K] fmlt\n')%(mu,a,GM,LIC,np.log10(LIC/Lsun))
np.savetxt(OUT+'ic_t9l62.txt',np.c_[re,de,eint,Fr,E,Te,np.zeros_like(re)],fmt='%.12e',header=hdr.rstrip('\n'),comments='')
np.savetxt(OUT+'ic_t9l62_LIC.txt',np.c_[re,de,eint,LIC/(4*np.pi*re**2),E,Te,np.zeros_like(re)],fmt='%.12e',header=hdr.rstrip('\n').replace('L_BASE/(4 pi r^2), L_BASE = 10^6.20 Lsun (the imposed base luminosity; the T column itself carries L_IC)','L_IC/(4 pi r^2) [VARIANT: base L = the IC column luminosity]'),comments='')
hdr2=('# he_bc_linhao profile (r[cm] T[K] rho[g/cc] M[g] kappa F_r[erg/cm2/s]) for the inner GHOSTS only: T, rho = the IC column (fixed),\n'
      '# F_r = L_BASE/(4 pi r^2), L_BASE = 10^6.20 Lsun (Lsun 3.828e33): Schultz+2020 Table 1 / Jiang+2018 fixed bottom luminosity. IC column L_IC = 10^%.4f\n# r T rho M kappa F_r')%np.log10(LIC/Lsun)
np.savetxt(OUT+'prof_t9l62_L620.txt',np.c_[re,Te,de,M,ke,LBASE/(4*np.pi*re**2)],fmt='%.12e',header=hdr2,comments='')
np.savetxt(OUT+'prof_t9l62_LIC.txt',np.c_[re,Te,de,M,ke,LIC/(4*np.pi*re**2)],fmt='%.12e',header=hdr2.replace('L_BASE/(4 pi r^2), L_BASE = 10^6.20 Lsun','L_IC/(4 pi r^2) [VARIANT: base flux = the IC column luminosity]'),comments='')
for rb in [24.2286,35.0]:
    print('r_in %.4f Rsun: F_in(L 6.20) = %.10e  F_in(L_IC) = %.10e'%(rb,LBASE/(4*np.pi*(rb*Rsun)**2),LIC/(4*np.pi*(rb*Rsun)**2)))
print('GM_for_deck %.10e'%GM)
