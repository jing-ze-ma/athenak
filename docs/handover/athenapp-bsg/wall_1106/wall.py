import numpy as np, glob, sys
Rs=6.957e10; Rg=8.314462618e7; a_r=7.5657e-15; c=2.99792458e10
T0=1e5; rho0=1e-9; mu=.62
vu=np.sqrt(Rg*T0/mu); tu=Rs/vu; Pu=rho0*vu**2; Eru=a_r*T0**4; Fu=Eru*c; GM=2.84560241e3*Rs**3/tu**2
Eu=Pu*Rs**3; dOm=np.pi/3; gam=5/3
def G1(Pg,Pr):
    b=Pg/(Pg+Pr); return b+(4-3*b)**2*(gam-1)/(b+12*(gam-1)*(1-b))
def N2(r,rho,Pg,Er):
    Pr=Er/3; P=Pg+Pr; g=GM/r**2; x=np.log(P); y=np.log(rho)
    dx=np.gradient(x,r,edge_order=1); dy=np.gradient(y,r,edge_order=1)
    return g*(dx/G1(Pg,Pr)-dy)
days=[0,.5,1,2,3,5,7,10]
out=open('wall_profiles.txt','w')
out.write('# Athena++ e47fabb7 + bsg_star.cpp, 1-D, cgs. r[Rsun] rho Pgas Er(lab) T=mu Pg/(R rho) vr Fr1(lab) Er0(comoving) Fr01(comoving) N2[s^-2]; first 24 active cells. dumps at the first cycle past each 0.5 d\n')
for a in ['A_af1','A_tc4']:
    print('==',a)
    h=np.loadtxt(f'{a}/out/bsg.hst'); fs=sorted(glob.glob(f'{a}/out/bsg.block0.out2.*.tab'))
    Ls=None
    for dd in days:
        if int(round(dd*2))>=len(fs): continue
        f=fs[int(round(dd*2))]; t=float(open(f).readline().split('time=')[1].split()[0])*tu
        if abs(t/86400-dd)>0.05: print(f'{a}: no valid dump at {dd} d (nearest t={t/86400:.2f} d)'); continue
        d=np.loadtxt(f); r=d[:,1]*Rs
        print(f'--- {a} t={t/86400:.3f} d')
        rho,Pg,v,Er,F1=d[:,2]*rho0,d[:,3]*Pu,d[:,4]*vu,d[:,7]*Eru,d[:,8]*Fu
        Er0,F01=d[:,12]*Eru,d[:,13]*Fu
        T=d[:,3]/d[:,2]*T0
        if dd==0: rho00,T00=rho.copy(),T.copy()
        n=N2(r,rho,Pg,Er)
        # coarsen to Caltech dr 0.3125 (pairs of cells: simple means)
        rc=(r[:-1:2]+r[1::2])/2
        def cg(q): return .5*(q[:-1:2]+q[1::2])
        nc=N2(rc,cg(rho),cg(Pg),cg(Er))
        # face N2 between cells i,i+1: no odd-even aliasing
        Pr=Er/3;P=Pg+Pr;dr=r[1]-r[0];gf=GM/((r[:-1]+r[1:])/2)**2
        G=G1(.5*(Pg[:-1]+Pg[1:]),.5*(Pr[:-1]+Pr[1:]))
        nf=gf*(np.diff(np.log(P))/dr/G-np.diff(np.log(rho))/dr)
        hse=(np.diff(P)/dr+.5*(rho[:-1]+rho[1:])*gf)/(.5*(rho[:-1]+rho[1:])*gf)
        print('   rho/IC : '+' '.join(f'{rho[i]/rho00[i]:7.3f}' for i in range(12)))
        print('   T/IC   : '+' '.join(f'{T[i]/T00[i]:7.3f}' for i in range(12)))
        print('   Nf e-10: '+' '.join(f'{nf[i]*1e10:7.1f}' for i in range(12)))
        print('   a/g fac: '+' '.join(f'{hse[i]:7.2f}' for i in range(12)))
        print(f'   n(N2f<0)/24={np.sum(nf[:24]<0)}  rms N2f cells 0-23={np.sqrt(np.mean(nf[:24]**2))*1e10:.1f}e-10; 24-96 ={np.sqrt(np.mean(nf[24:96]**2))*1e10:.1f}e-10')
        dV=4*np.pi*r**2*(r[1]-r[0])
        if Ls is None: Ls=np.median(4*np.pi*r[5:150]**2*F1[5:150])
        out.write(f'\n# {a} t={t:.1f} s = {t/86400:.4f} d\n')
        for i in range(24):
            out.write(f'{i:3d} '+' '.join(f'{x:.17e}' for x in [r[i]/Rs,rho[i],Pg[i],Er[i],T[i],v[i],F1[i],Er0[i],F01[i],n[i]])+'\n')
        print(f't={t/86400:5.2f} d  rmin F/L*: '+' '.join(f'{4*np.pi*r[i]**2*F1[i]/Ls:+7.1f}' for i in [0,1,3,7]))
        print('   T/1e6  : '+' '.join(f'{T[i]/1e6:7.3f}' for i in range(12)))
        print('   vr km/s: '+' '.join(f'{v[i]/1e5:7.2f}' for i in range(12)))
        print('   N2 e-10: '+' '.join(f'{n[i]*1e10:7.1f}' for i in range(12)))
        print('   N2 dr.3125 at r='+' '.join(f'{rc[i]/Rs:.2f}' for i in range(6)))
        print('   N2 e-10: '+' '.join(f'{nc[i]*1e10:7.1f}' for i in range(6)))
        # Global energy (wedge -> sphere)
        Etot=np.sum((Pg/(gam-1)+.5*rho*v**2+Er-GM/r*rho)*np.diff(np.r_[r-(r[1]-r[0])/2,r[-1]+(r[1]-r[0])/2]**3)/3*4*np.pi)
        M=np.sum(rho*np.diff(np.r_[r-(r[1]-r[0])/2,r[-1]+(r[1]-r[0])/2]**3)/3*4*np.pi)
        if dd==0: E00=Etot
        print(f'   dEtot={Etot-E00:+.4e} erg full sphere; mean rate {(Etot-E00)/max(t,1)/Ls:+.1f} L*; Mass={M/1.989e33:.6f} Msun')
out.close()
