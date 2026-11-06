import numpy as np,glob,sys
Rs=6.957e10;Rg=8.314462618e7;a_r=7.5657e-15;gam=5/3;mu=.62
def G1(Pg,Pr):
    b=Pg/(Pg+Pr);return b+(4-3*b)**2*(gam-1)/(b+12*(gam-1)*(1-b))
arms=[('af1_cflrad  (old units/padded tables, cfl_rad .025)','af1_cflrad',1e-9,1e5,2.84560241e3),
('his1d       (his units/tables, nzeta5, cfl_rad .025)','his1d',6.15245e-5,1.07887e6,260.795),
('his1d_cfl1  (his units/tables, nzeta5, cfl_rad 1)','his1d_cfl1',6.15245e-5,1.07887e6,260.795)]
for lab,d,r0,T0,gm in arms:
    vu=np.sqrt(Rg*T0/mu);tu=Rs/vu;Pu=r0*vu**2;Eru=a_r*T0**4;GM=gm*Rs**3/tu**2
    fs=sorted(glob.glob(f'{d}/out/bsg.block0.out2.*.tab'));print('==',lab,'dumps',len(fs))
    log=open(f'{d}/run.log').read(); print('   jacobi stops at 101:',log.count('niter: 101'),' total reports',log.count('Iteration stops'), 'fatal',log.count('FATAL'))
    h=np.loadtxt(f'{d}/out/bsg.hst'); print(f'   hst reaches {h[-1,0]*tu/86400:.2f} d, dt {h[-1,1]*tu:.1f} s, nan {np.isnan(h).any()}')
    tf=[];
    for f in fs: tf.append(float(open(f).readline().split('time=')[1].split()[0])*tu/86400)
    tf=np.array(tf);E0=None
    for dd in [0,.5,1,2,3,5,7,10]:
        k=np.argmin(abs(tf-dd))
        if abs(tf[k]-dd)>.06:continue
        q=np.loadtxt(fs[k])
        if np.isnan(q).any():print(f'   {dd} d NaN');continue
        r=q[:,1]*Rs;rho=q[:,2]*r0;Pg=q[:,3]*Pu;v=q[:,4]*vu;Er=q[:,7]*Eru;T=q[:,3]/q[:,2]*T0
        if dd==0:T00,R00=T.copy(),rho.copy()
        dr=r[1]-r[0];Pr=Er/3;P=Pg+Pr;g=GM/((r[:-1]+r[1:])/2)**2;Ga=G1(.5*(Pg[:-1]+Pg[1:]),.5*(Pr[:-1]+Pr[1:]))
        n=g*(np.diff(np.log(P))/dr/Ga-np.diff(np.log(rho))/dr)
        e=np.diff(np.r_[r-dr/2,r[-1]+dr/2]**3)/3*4*np.pi
        Et=np.sum((Pg/(gam-1)+.5*rho*v**2+Er-GM/r*rho)*e)
        if E0 is None:E0=Et
        print(f'   {tf[k]:5.2f} d T/IC '+' '.join(f'{T[i]/T00[i]:.3f}' for i in range(4))+' rho/IC '+' '.join(f'{rho[i]/R00[i]:.2f}' for i in range(3))+' N2f '+' '.join(f'{n[i]*1e10:+.0f}' for i in range(3))+f' rms(0-23) {np.sqrt(np.mean(n[:24]**2))*1e10:.0f} KE(pi/3) {np.sum(.5*rho*v**2*e)/12:.2e} dE {(Et-E0)/1e45:+.2f}e45 Lc0/L* {4*np.pi*r[0]**2*q[0,8]*Eru*2.998e10/5.58e38:+.0f}')
