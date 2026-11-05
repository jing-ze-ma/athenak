import numpy as np, os
O='/raven/ptmp/jinma/athpp_audit_1015/bsg1d/'
RS=6.957e10; ru=1e-9; Tu=1e5; lu=RS; mu=0.62
R=8.314462618e7/mu; a=7.565733250033928e-15; c=2.99792458e10; G=6.67408e-8
cs=np.sqrt(R*Tu); tu=lu/cs
prat=a*Tu**3/(ru*R); crat=c/cs
d=np.loadtxt('/raven/ptmp/jinma/bsglinhao_1014/files/MESA_profile_combined_cgs.txt')
out=np.column_stack([d[:,0]/lu, d[:,1]/Tu, d[:,2]/ru, d[:,3]/(ru*lu**3), d[:,4]*ru*lu, d[:,5]/(c*a*Tu**4)])
np.savetxt(O+'MESA_profile_combined.txt', out, fmt='%.12e')
def tab(fn):
    for ln in open(fn):
        if ln.startswith('# ') and len(ln.split())==7 and ln.split()[1].isdigit():
            nT,nD,lT0,dlT,lD0,dlD=[float(v) for v in ln.split()[1:]]
    v=np.loadtxt(fn,comments='#').reshape(int(nT),int(nD)); return v,lT0,dlT,lD0,dlD
D='/raven/ptmp/jinma/bsg_raven_1002/ma2026_tables/'
vr,lT0,dlT,lD0,dlD=tab(D+'rosseland_ma2026_x0.7_z0.008_padded.txt')
vp,lT0p,dlTp,_,_=tab(D+'planck_ma2026_x0.7_z0.008_padded.txt')
assert abs(lT0-lT0p)<1e-9 and vp.shape==vr.shape
sel=np.arange(0,vr.shape[0],10)
lT=lT0+dlT*sel; lD=lD0+dlD*np.arange(vr.shape[1])
np.savetxt(O+'log_T.txt',lT,fmt='%.9f'); np.savetxt(O+'log_rho.txt',lD,fmt='%.9f')
np.savetxt(O+'log_kappa_ross.txt',vr[sel],fmt='%.7f'); np.savetxt(O+'log_kappa_planck.txt',vp[sel],fmt='%.7f')
np.savetxt(O+'log_kappa_electron.txt',np.full(vr[sel].shape,np.log10(0.34)),fmt='%.7f')
GM=G*20*1.98892e33*tu**2/lu**3
print('tu %.6e s prat %.6f crat %.3f GM_code %.8e  min kappa_R in table %.3f'%(tu,prat,crat,GM,10**vr.min()))
print('day in code %.6f'%(86400/tu))
