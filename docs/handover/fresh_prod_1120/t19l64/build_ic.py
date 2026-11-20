# build_ic.py: T19L6.4 IC for AthenaK he_star_m1 from Y.-F. Jiang init_80.txt (r/Rsun, T/1.6668e5 K, rho/5e-9 g/cc).
# Clone of t9l62_1113/jobs/build_ic.py (USER 11-13 rule for T9 "follow the ic file": L_base = the IC column's L).
# Outputs (cgs): ic/ic_t19l64.txt  he_ic_file, he_ic_cols 5 (+T, fmlt=0): r rho eint F_r E | T fmlt
# Below the file bottom (14.8292 Rsun) the column is extended by RK4 HSE + diffusion (ghosts/margin only).
import numpy as np
OUT='/viper/ptmp2/jinma/t19l64_1119/ic/'
SRC='/viper/ptmp2/jinma/schultz_ics_1113/x/initial_condition_help/IC_from_Yan-Fei/init_80.txt'
TAB='/viper/ptmp2/jinma/run12_viper_1022/in/rosseland_mesakap_gs98_x0.70_z0.02.txt'
T0=1.6668e5; rho0=5e-9; Rsun=6.955e10
mu=0.674; mu_u=1.66053906660e-24; kB=1.380649e-16
a=7.5657332e-15; c=2.99792458e10; G=6.67408e-8; Msun=1.989e33; Lsun=3.828e33; gam=5.0/3.0
def tab(fn):
    A=open(fn).read().splitlines(); gi=[i for i,l in enumerate(A) if l.startswith('# grid')][0]
    h=A[gi+1].lstrip('#').split(); nT,nD=int(h[0]),int(h[1]); lT0,dlT,lD0,dlD=map(float,h[2:6])
    L=[l for l in A if l.strip() and not l.startswith('#')]
    v=np.array([float(x) for x in L[:nT*nD]]).reshape(nT,nD)
    from scipy.interpolate import RegularGridInterpolator as R
    f=R((lT0+dlT*np.arange(nT),lD0+dlD*np.arange(nD)),v)
    return lambda T,rho: 10**f(np.c_[np.atleast_1d(np.log10(T)),np.atleast_1d(np.log10(rho))])
kap=tab(TAB)
d=np.loadtxt(SRC); d=d[np.r_[True,np.diff(d[:,0])>0]]
r=d[:,0]*Rsun; T=d[:,1]*T0; rho=d[:,2]*rho0
k=kap(T,rho); Pg=rho*kB*T/(mu*mu_u); Pr=a*T**4/3
GMr=-(np.gradient(Pg,r)+np.gradient(Pr,r))/rho*r**2
Lr=-4*np.pi*r**2*c/(k*rho)*np.gradient(Pr,r)
for lo,hi in [(14.9,16),(16,18),(18,20),(20,24),(24,28),(28,30)]:
    s=(r>lo*Rsun)&(r<hi*Rsun)
    print('%5.1f-%5.1f GM %.4f Msun  logL %.4f'%(lo,hi,np.median(GMr[s])/G/Msun,np.log10(np.median(Lr[s])/Lsun)))
sel=(r>16*Rsun)&(r<28*Rsun)&(np.gradient(rho,r)<0)
GM=np.median(GMr[sel]); LIC=np.median(Lr[sel])
print('FIT 16-28 Rsun: GM %.10e (%.4f Msun) spread %.2e | L_IC %.6e (log %.4f) spread %.2e'%(GM,GM/G/Msun,
      np.std(GMr[sel])/GM,LIC,np.log10(LIC/Lsun),np.std(Lr[sel])/LIC))
def rhs(rr,y):
    pg,t=y; rh=pg*mu*mu_u/(kB*t); kk=kap(t,rh)[0]
    dpr=-kk*rh*LIC/(4*np.pi*rr**2*c); dpt=-GM*rh/rr**2
    return np.array([dpt-dpr,dpr/(4*a*t**3/3)])
y=np.array([Pg[0],T[0]]); rr=r[0]; h=-(r[1]-r[0])*50; ext=[]
while rr>12.0*Rsun:
    k1=rhs(rr,y);k2=rhs(rr+h/2,y+h/2*k1);k3=rhs(rr+h/2,y+h/2*k2);k4=rhs(rr+h,y+h*k3)
    y=y+h/6*(k1+2*k2+2*k3+k4); rr=rr+h; ext.append((rr,y[1],y[0]*mu*mu_u/(kB*y[1])))
ext=np.array(ext[::-1])
re=np.r_[ext[:,0],r]; Te=np.r_[ext[:,1],T]; de=np.r_[ext[:,2],rho]
print('extension: %d rows down to %.4f Rsun, T %.4e rho %.4e'%(len(ext),re[0]/Rsun,Te[0],de[0]))
eint=de*kB*Te/(mu*mu_u)/(gam-1.0); E=a*Te**4; Fr=LIC/(4*np.pi*re**2)
hdr=('# he_star_m1 IC (he_ic_cols 5, 7 columns): T19L6.4 (Jiang+2018 Nature 561 498; Schultz+2020 Table 1) from Y.-F. Jiang init_80.txt\n'
     '# source units: r/Rsun (6.955e10), T/1.6668e5 K, rho/5e-9 g/cc; ideal gas gamma 5/3 mu %.3f; a %.7e\n'
     '# below 14.8292 Rsun: RK4 HSE + diffusion with GM %.6e (%.3f Msun), L_IC %.5e (log %.4f), kappa_R mesakap gs98 X0.70 Z0.02\n'
     '# F_r = L_IC/(4 pi r^2) (follow the IC file: S20 Table 1 log L 6.40 NOT used); E = a T^4\n'
     '# r[cm] rho eint F_r E | T[K] fmlt')%(mu,a,GM,GM/G/Msun,LIC,np.log10(LIC/Lsun))
np.savetxt(OUT+'ic_t19l64_LIC.txt',np.c_[re,de,eint,Fr,E,Te,np.zeros_like(re)],fmt='%.12e',header=hdr,comments='')
print('GM_for_deck %.10e  LIC %.10e  F_in(14.8292) %.10e'%(GM,LIC,LIC/(4*np.pi*(re[len(ext)])**2)))
