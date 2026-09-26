# usage: hst.py <run dir> [nrows] [gamma]: sph_wedge history summary
import sys, glob, numpy as np
d=sys.argv[1]; nr=int(sys.argv[2]) if len(sys.argv)>2 else 12; gam=float(sys.argv[3]) if len(sys.argv)>3 else 5/3
f=glob.glob(d+'/*.user.hst')[0]
hdr=[l for l in open(f) if l.startswith('#')][-1]
names=[s.split('=')[1] for s in hdr[1:].split()]
a=np.loadtxt(f); a=np.atleast_2d(a)
c={n:a[:,i] for i,n in enumerate(names)}
t=c['time']; L0=c['L_in']
mach=np.sqrt(2*c['KE_int']/(gam*c['PV_int'])); machr=np.sqrt(2*c['KEr_int']/(gam*c['PV_int']))
vr=c['Mr_int']/c['M_int']
idx=np.unique(np.linspace(0,len(t)-1,nr).astype(int))
print("%10s %7s %10s %10s %10s %10s %10s %10s %10s %11s %11s"%("time","row","Lbot/Lin","Lmid/Lin","Lint/Lin","Ltop/Lin","Mach_rms","Mach_r","<vr>/cs","dErad/E0","degas/e0"))
cs=np.sqrt(gam*c['PV_int']/c['M_int'])
for i in idx:
    print("%10.4e %7d %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %11.3e %11.3e"%(t[i],i,c['L_bot'][i]/L0[i],c['L_mid'][i]/L0[i],c['L_int'][i]/L0[i],c['L_top'][i]/L0[i],mach[i],machr[i],vr[i]/cs[i],c['E_rad'][i]/c['E_rad'][0]-1,c['e_gas'][i]/c['e_gas'][0]-1))
