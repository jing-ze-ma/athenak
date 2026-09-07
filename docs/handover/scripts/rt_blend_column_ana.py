import numpy as np, sys
sig=5.670374419e-5
def load(d):
    a=np.loadtxt(d+'/col.txt'); hdr=open(d+'/col.txt').readline()
    return a
def freedman(T,p,met=0.0):
    T1=np.clip(T,75,4000); p1=np.clip(p,1,3e8); lT=np.log10(T1); lp=np.log10(p1)
    c1,c2,c3,c4,c5,c6,c7=10.602,2.882,6.09e-15,2.954,-2.526,0.843,-5.490
    lo=T1<800
    c8=np.where(lo,-14.051,82.241);c9=np.where(lo,3.055,-55.456);c10=np.where(lo,0.024,8.754)
    c11=np.where(lo,1.877,0.7048);c12=np.where(lo,-0.445,-0.0414);c13=0.8321
    lkl=c1*np.arctan(lT-c2)-c3/(lp+c4)*np.exp((lT-c5)**2)+c6*met+c7
    lkh=c8+c9*lT+c10*lT**2+lp*(c11+c12*lT)+c13*met*(0.5+np.arctan((lT-2.5)/0.2)/np.pi)
    return 10**lkl+10**lkh
Fint=4.7352e6
for name in sys.argv[1:]:
    a=load(S+'/'+name) if False else np.loadtxt(name+'/col.txt')
    i,r,p,T,Flw,Qsw,g1,gad=a[:,0],a[:,1],a[:,2],a[:,3],a[:,4],a[:,5],a[:,6],a[:,7]
    tau=a[:,8] if a.shape[1]>8 else np.zeros_like(r); w=a[:,9] if a.shape[1]>9 else np.zeros_like(r)
    hdr=open(name+'/col.txt').readlines()[1]; icut=int(hdr.split('icut = ')[1].split()[0])
    print('=== %s  icut=%d'%(name,icut))
    # face quantities: r is the face radius; p,T are LEVEL (face) values in the ck dump
    n=len(r)
    # diffusive flux at interior faces from face T gradient (cell T = average of level T)
    # estimate kappa_R*rho at the faces via Freedman on level (p,T) and rho from p/(R T)? use tau differences instead
    dtau=-(np.diff(tau))  # tau(face i) - tau(face i+1) = kappa rho dr of cell i  (>0)
    dr=np.diff(r)
    krho_cell=dtau/dr
    # cell temperatures and face gradient
    Tc=0.5*(T[1:]+T[:-1]); rc=0.5*(r[1:]+r[:-1])
    Fd=np.zeros(n); 
    for f in range(1,n-1):
        kr=0.5*(krho_cell[f-1]+krho_cell[f])
        if kr<=0: continue
        Tf=T[f]; dTdr=(Tc[f]-Tc[f-1])/(rc[f]-rc[f-1])
        F=-16*sig*Tf**3/(3*kr)*dTdr
        Fd[f]=F/np.sqrt(1+(F/(sig*Tf**4))**2)
    Ftot=(1-w)*Flw+w*Fd
    print(' i  p[bar]      T[K]    tau_R    w      F_lw/Fint  F_diff/Fint  F_tot/Fint')
    sel=[k for k in range(n) if (name!='old' and 0<w[k]<1) or k in (n-1,n//2,icut-1,icut,icut+1,0,1,2,3)]
    for k in sorted(set(sel)):
        print('%3d %10.3e %8.1f %9.3g %6.3f %9.4f %10.4f %10.4f'%(i[k],p[k],T[k],tau[k],w[k],Flw[k]/Fint,Fd[k]/Fint,Ftot[k]/Fint))
    # source per cell from the dumped fluxes: -(Ftot(i+1)-Ftot(i))/dr + Qsw   [Flw dumped is unweighted]
    src=-(Ftot[1:]-Ftot[:-1])/dr
    if name!='old':
        ov=np.where((w[:-1]>0)&(w[:-1]<1))[0]
        print(' overlap faces:',ov.min(),ov.max(),' p range %.3g-%.3g bar'%(p[ov.max()],p[ov.min()]))
        print(' F_tot/Fint in overlap: min %.4f max %.4f'%((Ftot[ov]/Fint).min(),(Ftot[ov]/Fint).max()))
        lo=max(ov.min()-3,0); hi=min(ov.max()+3,n-2)
        print(' cell source [erg/s/cm3] around the overlap (cells %d..%d):'%(lo,hi)); print(np.array2string(src[lo:hi+1],precision=3))
    else:
        print(' cell source around icut: ',np.array2string(src[max(icut-4,0):icut+4],precision=3))
# explicit diffusion dt per cell in the new arm
a=np.loadtxt('new/col.txt'); r,p,T,tau,w=a[:,1],a[:,2],a[:,3],a[:,8],a[:,9]
dr=np.diff(r); krho=-np.diff(tau)/dr; Tc=0.5*(T[1:]+T[:-1]); pc=np.sqrt(p[1:]*p[:-1])*1e6
mu=2.3; kB=1.380649e-16; mH=1.6726e-24; rho=pc*mu*mH/(kB*Tc); cv=1.5*kB/(mu*mH)  # rough
D=16*sig*Tc**3/(3*krho*rho*cv); dtd=0.5*dr**2/D
wc=0.5*(w[1:]+w[:-1]); m=wc>0
print('\nexplicit diffusion dt [s] where w>0: min %.3g at cell %d (p=%.3g bar, T=%.0f, w=%.2f); '%(dtd[m].min(), np.argmin(np.where(m,dtd,1e99)), pc[np.argmin(np.where(m,dtd,1e99))]/1e6, Tc[np.argmin(np.where(m,dtd,1e99))], wc[np.argmin(np.where(m,dtd,1e99))]))
for k in np.where(m)[0][::6]: print('  cell %3d p=%9.3g T=%6.0f w=%.2f  dr=%.2e  D=%.3g  dt_diff=%.3g'%(k,pc[k]/1e6,Tc[k],wc[k],dr[k],D[k],dtd[k]))
