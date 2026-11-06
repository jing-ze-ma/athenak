import numpy as np,glob,re,sys,os
Rs=6.957e10;Rg=8.314462618e7;gam=5/3;mu=.62;yr=3.156e7;c=2.99792458e10
r0,T0,gm=6.15245e-5,1.07887e6,260.795;vu=np.sqrt(Rg*T0/mu);tu=Rs/vu;Pu=r0*vu**2;Eru=7.5657e-15*T0**4;GM=gm*Rs**3/tu**2;Eu=Pu*Rs**3;LS=5.581e38
def G1(Pg,Pr):
    b=Pg/(Pg+Pr);return b+(4-3*b)**2*(gam-1)/(b+12*(gam-1)*(1-b))
def load(d,k):
    fs=sorted(glob.glob(f'{d}/out/bsg.block*.out2.{k:05d}.tab'),key=lambda f:int(re.search(r'block(\d+)',f).group(1)))
    t=float(open(fs[0]).readline().split('time=')[1].split()[0])*tu/86400
    q=np.vstack([np.loadtxt(f) for f in fs]);q=q[np.argsort(q[:,1])];return t,q
for d in sys.argv[1:]:
    h=np.loadtxt(d+'/out/bsg.hst');t=h[:,0]*tu/86400;k=(h[:,6]+h[:,7]+h[:,8])*Eu
    np.savetxt(d+'_ke.txt',np.column_stack([t,k]),fmt='%.6e',header=f'{d}: t[d] KE[erg, pi/3 sr], hst cadence (0.1 code = 578 s)')
    print('=====',d,'t_end %.2f dt_end %.2f s nan %s; KE max %.2e at %.2f d'%(t[-1],h[-1,1]*tu,np.isnan(h).any(),k.max(),t[k.argmax()]))
    for a in range(0,int(t[-1]+.5),2):
        s=(t>=a)&(t<a+2)
        if s.any():print(f' KE {a:2d}-{a+2:2d} median {np.median(k[s]):.1e} min {k[s].min():.1e} max {k[s].max():.1e}')
    s=t>2;kk=k[s];mins=np.where((kk[1:-1]<kk[:-2])&(kk[1:-1]<kk[2:])&(kk[1:-1]<0.5*np.median(kk)))[0];print(' KE dips [d]',np.round(t[s][mins][:30],1))
    n=np.array(re.findall(r'niter: (\d+)',open(d+'/run.log').read()),int);print(' Jacobi solves',n.size,'mean %.2f p99 %.0f max %d capped %d'%(n.mean(),np.percentile(n,99),n.max(),(n>=100).sum()),'; FATAL',open(d+'/run.log').read().count('FATAL'))
    nd=len(glob.glob(f'{d}/out/bsg.block0.out2.*.tab'))
    out=open(d+'_wall_profiles.txt','w');out.write(f'# {d} (1-D, his radial grid via SMR, his deck units/tables/GM, cfl_rad .025, nzeta 5), cgs: i r[Rsun] rho Pgas Er(lab) T vr Fr1(lab) N2_face(i,i+1)[s^-2]\n')
    E0=None
    print(' t  T/IC c0-3 | rho/IC c0-2 | N2f f0-2 | rms N2f 0-23 | dEgrav dU [1e44,pi/3] | dE[1e45] | Lc0/L* | Mdot/Menv max r<25,25-35,35-45 | vmax r<45 | R(1e-10)')
    for k in range(nd):
        tt,q=load(d,k)
        if np.isnan(q).any():print(tt,'NaN');break
        r=q[:,1]*Rs;rho=q[:,2]*r0;Pg=q[:,3]*Pu;v=q[:,4]*vu;Er=q[:,7]*Eru;F=q[:,8]*Eru*c;T=q[:,3]/q[:,2]*T0
        if k==0:T00,R00=T.copy(),rho.copy()
        Pr=Er/3;P=Pg+Pr;rc=.5*(r[:-1]+r[1:]);g=GM/rc**2;G=G1(.5*(Pg[:-1]+Pg[1:]),.5*(Pr[:-1]+Pr[1:]));nf=g*(np.diff(np.log(P))/np.diff(r)/G-np.diff(np.log(rho))/np.diff(r))
        ff=np.r_[r[0]-(r[1]-r[0])/2,rc,r[-1]+(r[-1]-r[-2])/2];V=np.diff(ff**3)/3*(np.pi/3)
        Eg=np.sum(-GM/r*rho*V);U=np.sum((Pg/(gam-1)+Er)*V);Et=np.sum((Pg/(gam-1)+.5*rho*v**2+Er-GM/r*rho)*V)*12
        if E0 is None:E0=(Eg,U,Et)
        Me=np.sum(rho*V)*12;Md=4*np.pi*r**2*rho*v/Me*yr;R=r/Rs;mx=lambda a,b:np.abs(Md[(R>a)&(R<b)]).max()
        out.write(f'\n# t={tt:.4f} d\n')
        for j in range(24):out.write(f'{j:3d} '+' '.join(f'{x:.17e}' for x in [R[j],rho[j],Pg[j],Er[j],T[j],v[j],F[j],nf[j]])+'\n')
        if abs(tt-round(tt))<.1 and round(tt) in (0,1,2,3,5,7,10,14,16,18,20,22,25,28,30) or abs(tt-.5)<.1:
            print(f'{tt:5.1f} '+' '.join(f'{T[j]/T00[j]:.3f}' for j in range(4))+' | '+' '.join(f'{rho[j]/R00[j]:.2f}' for j in range(3))+' | '+' '.join(f'{nf[j]*1e10:+.0f}' for j in range(3))+f' | {np.sqrt(np.mean(nf[:24]**2))*1e10:.0f} | {(Eg-E0[0])/1e44:+.2f} {(U-E0[1])/1e44:+.2f} | {(Et-E0[2])/1e45:+.2f} | {4*np.pi*r[0]**2*F[0]/LS:+.0f} | {mx(20.5,25):.1f} {mx(25,35):.1f} {mx(35,45):.2f} | {np.abs(v[R<45]).max()/1e5:.1f} | {R[np.where(rho>1e-10)[0].max()]:.1f}')
    print(' grid check: nx',len(r),'dr near 20:',round(R[1]-R[0],5),'near 45:',round(R[np.argmin(abs(R-45))+1]-R[np.argmin(abs(R-45))],5),'near 70:',round(R[np.argmin(abs(R-70))+1]-R[np.argmin(abs(R-70))],5))
