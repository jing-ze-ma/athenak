import numpy as np,glob
Rs=6.957e10;Rg=8.314462618e7;a_r=7.5657e-15;T0=1e5;rho0=1e-9;mu=.62
vu=np.sqrt(Rg*T0/mu);tu=Rs/vu;Pu=rho0*vu**2;Eru=a_r*T0**4;GM=2.84560241e3*Rs**3/tu**2;gam=5/3;LS=5.581e38
def G1(Pg,Pr):
    b=Pg/(Pg+Pr);return b+(4-3*b)**2*(gam-1)/(b+12*(gam-1)*(1-b))
def load(f):
    t=float(open(f).readline().split('time=')[1].split()[0])*tu;d=np.loadtxt(f);return t,d
def N2f(r,rho,Pg,Er):
    Pr=Er/3;P=Pg+Pr;dr=r[1]-r[0];g=GM/((r[:-1]+r[1:])/2)**2;G=G1(.5*(Pg[:-1]+Pg[1:]),.5*(Pr[:-1]+Pr[1:]))
    return g*(np.diff(np.log(P))/dr/G-np.diff(np.log(rho))/dr)
for a in ['A_af1','A_tc4']:
  for lab,D in [('SEEDED',f'/viper/ptmp2/jinma/athpp_ke_1106/wall_1106/{a}/out'),('NOSEED',f'{a}/noseed/out')]:
    fs=sorted(glob.glob(f'{D}/bsg.block0.out2.*.tab'));print('==',a,lab,'dumps',len(fs))
    E0=None
    for dd in [0,.5,1,2,3,5,7,10]:
        k=int(round(dd*2))
        if k>=len(fs):continue
        t,d=load(fs[k])
        if abs(t/86400-dd)>.05 or np.isnan(d).any():print(f'  {dd} d: no valid dump');continue
        r=d[:,1]*Rs;rho=d[:,2]*rho0;Pg=d[:,3]*Pu;v=d[:,4]*vu;Er=d[:,7]*Eru;T=d[:,3]/d[:,2]*T0
        if dd==0:r0,T0a,rho0a=r,T.copy(),rho.copy()
        n=N2f(r,rho,Pg,Er);dr=r[1]-r[0];e=np.diff(np.r_[r-dr/2,r[-1]+dr/2]**3)/3*4*np.pi
        ke=np.sum(.5*rho*v**2*e)/12
        Et=np.sum((Pg/(gam-1)+.5*rho*v**2+Er-GM/r*rho)*e)
        if E0 is None:E0=Et
        print(f'  {t/86400:5.2f} d T/IC '+' '.join(f'{T[i]/T0a[i]:.3f}' for i in range(5))+' rho/IC '+' '.join(f'{rho[i]/rho0a[i]:.2f}' for i in range(3))+' N2f '+' '.join(f'{n[i]*1e10:+.0f}' for i in range(3))+f' rmsN2f(0-23) {np.sqrt(np.mean(n[:24]**2))*1e10:.1f} KE(wedge) {ke:.2e} dE {(Et-E0)/1e45:+.2f}e45 maxv(r<45) {np.max(abs(v[r/Rs<45]))/1e5:.1f}')
