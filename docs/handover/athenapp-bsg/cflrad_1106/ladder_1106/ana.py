import numpy as np,glob,os
Rs=6.957e10;Rg=8.314462618e7;mu=.62;gam=5/3;a_r=7.5657e-15;c=2.99792458e10
def u(rho,T):
    vu=np.sqrt(Rg*T/mu);return dict(rho=rho,T=T,vu=vu,tu=Rs/vu,Pu=rho*vu**2,Eru=a_r*T**4,Eu=rho*vu**2*Rs**3)
H=u(6.15245e-5,1.07887e6);A=u(1e-9,1e5)
arms=[('H_base',H,'his1d baseline'),('U_units',A,'units only (physics same)'),('G_gm20',H,'GM = 20 Msun'),('Z_nzeta4',H,'nzeta 4'),('F_floors',H,'floors = old'),('T_tables',H,'padded tables'),('S_seed',H,'seed thresh = old'),('ALL_his_units',H,'G+Z+F+T+S'),('ALLnoG',H,'Z+F+T+S, his GM'),('../af1_cflrad',A,'af1_cflrad')]
def KEat(d,U,w):
    h=np.loadtxt(d+'/out/bsg.hst');t=h[:,0]*U['tu']/86400;i=np.argmin(abs(t-w));return (h[i,6]+h[i,7]+h[i,8])*U['Eu'],t[i],t[-1]
for d,U,lab in arms:
    if not os.path.exists(d+'/out/bsg.hst'):continue
    k=[KEat(d,U,w) for w in [.1,.25,.5,1]]
    fs=sorted(glob.glob(d+'/out/bsg.block0.out2.*.tab'))
    ts=np.array([float(open(f).readline().split('time=')[1].split()[0])*U['tu']/86400 for f in fs])
    GM={'G_gm20':2845.60241/10.7887,'ALL_his_units':2845.60241/10.7887}.get(d,260.795 if U is H else (2845.60241 if 'af1' in d else 260.795*10.7887))*Rs**3/U['tu']**2
    def st(w):
        q=np.loadtxt(fs[np.argmin(abs(ts-w))]);r=q[:,1]*Rs;rho=q[:,2]*U['rho'];P=q[:,3]*U['Pu'];v=q[:,4]*U['vu'];Er=q[:,7]*U['Eru'];F=q[:,8]*U['Eru']*c
        dr=r[1]-r[0];V=np.diff(np.r_[r-dr/2,r[-1]+dr/2]**3)/3*4*np.pi
        E=np.sum((P/(gam-1)+.5*rho*v**2+Er-GM/r*rho)*V);T=q[:,3]/q[:,2]*U['T']
        return dict(T=T,E=E,L=4*np.pi*r[0]**2*F[0]/5.58e38,v=v,r=r/Rs,rho=rho,P=P,Er=Er,V=V)
    s0=st(0);s5=st(.5);s1=st(1)
    m=np.argmax(abs(s5['v'])*(s5['r']<45))
    ek=.5*s5['rho']*s5['v']**2*s5['V']/12;c1=np.cumsum(ek)/ek.sum();r50=s5['r'][np.searchsorted(c1,.5)]
    print(f"{lab:26s} tend {k[0][2]:.2f} KE .1/.25/.5/1 d {k[0][0]:.2e} {k[1][0]:.2e} {k[2][0]:.2e} {k[3][0]:.2e} | 1d T/IC {s1['T'][0]/s0['T'][0]:.3f} {s1['T'][1]/s0['T'][1]:.3f} Lc0 {s1['L']:+.0f} dE {(s1['E']-s0['E'])/1e45:+.2f}e45 | .5d max|v| {s5['v'][m]/1e5:+.2f} km/s at {s5['r'][m]:.1f} KE r50 {r50:.1f}")
# HSE residual of the IC for each GM (r 22-40 Rsun, smoothed by a 10-cell fit to remove the seed): (dP/dr)/(rho g)+1
q=np.loadtxt('H_base/out/bsg.block0.out2.00000.tab');r=q[:,1]*Rs;rho=q[:,2]*H['rho'];P=q[:,3]*H['Pu']+q[:,7]*H['Eru']/3
d=np.loadtxt('/viper/ptmp2/jinma/bsg_1001/ma2026_files/BSG/MESA_profile_combined_cgs.txt')
print('profile mass at r=20 Rsun [Msun]:',np.interp(20*Rs,d[:,0],d[:,3])/1.98847e33,' at 80:',np.interp(80*Rs,d[:,0],d[:,3])/1.98847e33)
for lab,gm in [('his GM 260.795',260.795*Rs**3/H['tu']**2),('old GM (20 Msun)',2845.60241/10.7887*Rs**3/H['tu']**2)]:
    w=(r/Rs>21)&(r/Rs<40);c0=np.polyfit(r[w],np.log(P[w]),6);dp=np.polyval(np.polyder(c0),r[w])*P[w]
    fr=np.polyval(np.polyfit(r[w],np.log(rho[w]),6),r[w]);res=dp/(np.exp(fr)*gm/r[w]**2)+1
    print(f'{lab} GM={gm:.4e} = {gm/6.6743e-8/1.98847e33:.3f} Msun; IC (dP_tot/dr)/(rho g)+1, r 21-40: median {np.median(res):+.4f} at 21: {res[0]:+.4f} at 25 {res[np.argmin(abs(r[w]/Rs-25))]:+.4f} at 40 {res[-1]:+.4f}')
