# Dent et al. 2026 (arXiv:2608.19339) App. B eq.(1)-(2) SEM; R*=812 Rsun (Rosseland, d=172 pc)
import numpy as np
Rsun=6.957e10; Rs=812*Rsun; mH=2.27e-24
nchrom=1.34e12; Rmin=1.144; H=0.06; nst=2.96e9; g=0.45; d=1.5
def nH(x): 
    x=np.asarray(x,float); c=np.where(x>=Rmin,nchrom*np.exp(-(x-Rmin)/H),np.nan)
    return c+nst*x**-2*(0.998-x**-g)**-d
def xe(x): return np.where(x<2, 2.5e-5+(5e-4-2.5e-5)*(x-Rmin)/(2-Rmin),5e-4)
# T read off Fig 5 (approx +-50 K)
TT={1.144:2520,1.2:1780,1.3:2870,1.5:2940,1.7:3200,2:3500,2.5:3500,3:3180,3.5:2620,4:2000,4.5:1530,5:1410}
base={1.2:(8.3e11,2570),1.5:(9.3e10,3040),2:(7.6e9,3800),3:(1.6e9,3290),4:(6.8e8,2660),5:(3.7e8,2260)}
Rb=1014*Rsun
print(" x   r[cm]     nH        rho      ne      T(fig)  | Dent24 nH(same x) ratio | Dent24 interp at same r_cm ratio")
xb=np.array(sorted(base)); nb=np.array([base[k][0] for k in xb])
for x in [1.144,1.2,1.3,1.5,1.7,2,2.5,3,3.5,4,4.5,5,7,10,20,30]:
    n=float(nH(x)); s=f"{x:5.2f} {x*Rs:.3e} {n:.3e} {n*mH:.3e} {n*float(xe(x)):.2e} {TT.get(x,'-')!s:>5}"
    if x in base: s+=f" | {base[x][0]:.2e} {n/base[x][0]:.2f}"
    xb_eq=x*Rs/Rb
    if xb[0]<=xb_eq<=xb[-1]:
        nbi=10**np.interp(np.log10(xb_eq),np.log10(xb),np.log10(nb)); s+=f" | r/R*(1014)={xb_eq:.2f} {nbi:.2e} {n/nbi:.2f}"
    print(s)
x=np.array([2,3,5,10,30.]); n=nH(x); print("local slope dlnn/dlnr:",[round(v,2) for v in np.gradient(np.log(nH(np.linspace(2,30,500))),np.log(np.linspace(2,30,500)))[[0,36,107,286,499]]])
# Mdot check
v=9e5; print("Mdot(9 km/s, 1.4 mu) Msun/yr at 30R*:", 4*np.pi*(30*Rs)**2*nH(30)*1.4*1.6726e-24*v/1.989e33*3.156e7)
