import numpy as np, sys
kB=1.380649e-16; mH=1.6605390666e-24; a=7.5657332e-15; c=2.99792458e10
def table(fn):
    hdr=[l for l in open(fn) if l.startswith('#')]
    i=[k for k,l in enumerate(hdr) if 'nT nD' in l][0]
    nT,nD,lT0,dlT,lD0,dlD=[float(x) for x in hdr[i+1][1:].split()]
    v=np.loadtxt(fn,comments='#').reshape(int(nT),int(nD))
    def k(T,rho):
        x=(np.log10(T)-lT0)/dlT; y=(np.log10(rho)-lD0)/dlD
        x=np.clip(x,0,nT-1.001); y=np.clip(y,0,nD-1.001)
        i=x.astype(int); j=y.astype(int); fx=x-i; fy=y-j
        return 10**((1-fx)*(1-fy)*v[i,j]+fx*(1-fy)*v[i+1,j]+(1-fx)*fy*v[i,j+1]+fx*fy*v[i+1,j+1])
    return k
name,ic,tab,mu,GM,L,x1min,x1max,efromT=sys.argv[1:10]
mu,GM,L,x1min,x1max=map(float,(mu,GM,L,x1min,x1max)); efromT=efromT=='1'
d=np.loadtxt(ic,comments='#'); r,rho,eint,F,E,T=d[:,:6].T
kap=table(tab)
m=(r>=x1min)&(r<=x1max)
Te=(2/3)*eint*mu*mH/(rho*kB)
Pg=rho*kB*T/(mu*mH) if efromT else (2/3)*eint
k=kap(T,rho)
dPg=np.gradient(Pg,r); dE=np.gradient(E,r)
g=GM/r**2
res=(dPg+rho*g-rho*k*F/c)/(rho*g)
Fd=-c*dE/(3*k*rho)
Lr=4*np.pi*r**2*F
tau=np.concatenate([np.cumsum((k*rho*np.gradient(r))[::-1])[::-1]])
# tau from x1max down
mm=r<=x1max; tt=np.zeros_like(r); 
s=(k*rho*np.gradient(r))*mm; tt=np.cumsum(s[::-1])[::-1]
iph=np.where(tt>=2/3)[0].max()
Gam=k*L/(4*np.pi*c*GM)
print(f"== {name}: rows {len(r)}, in domain {m.sum()}, r {r[m].min()/6.957e10:.3f}-{r[m].max()/6.957e10:.3f} Rsun (deck {x1min/6.957e10:.3f}-{x1max/6.957e10:.3f})")
print(f"  T(eint,mu)/T_col: med {np.median(Te[m]/T[m]):.5f} range {np.min(Te[m]/T[m]):.5f}-{np.max(Te[m]/T[m]):.5f}")
print(f"  E/aT^4: {np.median(E[m]/(a*T[m]**4)):.5f}  L(r)/L: {np.min(Lr[m]/L):.4f}-{np.max(Lr[m]/L):.4f}")
dense=m&(rho>1e3*rho[m].min())&(tt>3)
for lab,sel in [('below photosphere tau>3',dense)]:
    print(f"  HSE residual (dPg/dr + rho g - rho k F/c)/(rho g) [{lab}]: med {np.median(res[sel]):+.2e}, |max| {np.max(np.abs(res[sel])):.2e} at r {r[sel][np.argmax(np.abs(res[sel]))]/6.957e10:.2f}")
    q=Fd[sel]/F[sel]; print(f"  diffusion F_diff/F: med {np.median(q):.4f} range {q.min():.3f}-{q.max():.3f}")
    print(f"  Gamma_Edd(kR) range {Gam[sel].min():.3f}-{Gam[sel].max():.3f}")
print(f"  photosphere (tau_R=2/3 from x1max): r {r[iph]/6.957e10:.2f} Rsun, T {T[iph]:.0f} K, rho {rho[iph]:.2e}; Teff from L,r_ph {(L/(4*np.pi*r[iph]**2*5.6704e-5))**0.25:.0f} K")
print(f"  rho at x1max {rho[m][-1]:.2e}, min rho in domain {rho[m].min():.2e}")
