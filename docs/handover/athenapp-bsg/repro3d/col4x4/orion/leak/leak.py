#!/usr/bin/env python3
"""TASK-2026-11-06t: leak vs dt. Band-averaged ln(Er_1e-6/Er_1e-9), column at several cfl_rad, and 3-D wedge."""
import numpy as np,h5py,glob,re
ru,Tu,lu,mu=6.15245e-5,1.07887e6,6.957e10,.62
Rg=8.314462618e7/mu;vu=(Rg*Tu)**.5;tu=lu/vu;Eu=ru*vu*vu*lu**3;a=7.565733250033928e-15;Eru=a*Tu**4;GM=260.795*lu**3/tu**2;D=tu/86400;gam=5/3
Pu=ru*vu*vu
B=[(20,24),(24,35),(35,45)]
def ix(f): return {x.decode():i for i,x in enumerate(f.attrs['VariableNames'])}
def ld(g):
    with h5py.File(g,'r') as f:
        I=ix(f);h=f['hydro'];x1=f['x1v'][:];x1f=f['x1f'][:];x2f=f['x2f'][:];x3f=f['x3f'][:];T=f.attrs['Time']*D
        S=np.zeros(3);V=np.zeros(3);U=Eg=KE=KEr=KEl=0.
        for b in range(x1.shape[0]):
            dO=(np.cos(x2f[b,:-1])-np.cos(x2f[b,1:]))[None,:,None]*np.diff(x3f[b])[:,None,None]
            dV=dO*((x1f[b,1:]**3-x1f[b,:-1]**3)/3)[None,None,:];rc=x1[b][None,None,:]
            Er=h[I['Er'],b];rho=h[I['rho'],b];P=h[I['press'],b];v1=h[I['vel1'],b];v2=h[I['vel2'],b];v3=h[I['vel3'],b]
            for n,(l,u) in enumerate(B):
                m=(x1[b]>=l)&(x1[b]<u);S[n]+=np.sum(Er[:,:,m]*dV[:,:,m]);V[n]+=np.sum(dV[:,:,m])
            U+=np.sum(P/(gam-1)*dV)*Eu+np.sum(Er*dV)*Eru*lu**3;Eg+=np.sum(-rho*dV/rc)*GM*ru*lu**2
            KEr+=np.sum(.5*rho*v1**2*dV)*Eu;KEl+=np.sum(.5*rho*(v2**2+v3**2)*dV)*Eu
    return T,S/V,U,Eg,KEr,KEl
def series(d,pat):
    fs=sorted(glob.glob(d+pat));return [ld(g) for g in fs]
def stats(d):
    log=open(d+'/run.log',errors='ignore').read();n=np.array(list(map(int,re.findall(r'niter: (\d+)',log))))
    c=re.search(r'cpu time used\s*=\s*([\d.e+-]+)',log);cy=re.findall(r'cycle=(\d+)',log)
    return n.mean(),n.max(),int(cy[-1]),float(c.group(1)) if c else float('nan')
def comp(lab,d6,d9,pat,ts,o):
    s6=series(d6,pat);s9=series(d9,pat);h6=np.loadtxt(d6+'/output/bsg.hst');h9=np.loadtxt(d9+'/output/bsg.hst')
    t6=np.array([x[0] for x in s6]);t9=np.array([x[0] for x in s9])
    for T in ts:
        i=np.argmin(abs(t6-T));j=np.argmin(abs(t9-T))
        if abs(t6[i]-T)>.015 or abs(t9[j]-T)>.015: continue
        dts=np.interp(T,h6[:,0]*D,h6[:,1])*tu;Ns=T*86400/dts
        x=np.log(s6[i][1]/s9[j][1])
        print('%s t=%.3f dt=%.2f s | per step %s | per day %s | dE6 %.3e dE9 %.3e | KE6 %.3e KE9 %.3e | KElat6 %.3e KElat9 %.3e'%(lab,s6[i][0],dts,' '.join('%+.2e'%(y/Ns) for y in x),' '.join('%+.2e'%(y/T) for y in x),
          s6[i][2]+s6[i][3]+s6[i][4]+s6[i][5]-(s6[0][2]+s6[0][3]),s9[j][2]+s9[j][3]+s9[j][4]+s9[j][5]-(s9[0][2]+s9[0][3]),s6[i][4]+s6[i][5],s9[j][4]+s9[j][5],s6[i][5],s9[j][5]),file=o)
    print(lab,'stats 1e-6 mean,max,cyc,cpu',stats(d6),' 1e-9',stats(d9),file=o)
R='/orion/ptmp/jinma/athpp_repro3d/'
with open('leak_vs_dt.txt','w') as o:
    print('# r bands 20-24 24-35 35-45; volume-weighted Er in band; ln(Er_1e-6/Er_1e-9); steps N=t/dt(t). dE=U+Egrav+KE change, wedge cgs',file=o)
    ts=(0.05,0.1,0.25,0.5,1.0)
    comp('col cfl0.025',R+'col4x4/prof3d',R+'col4x4/el1e-9','/output/bsg.out2.*.athdf',ts,o)
    for cf in ('0.011','0.05','0.1','0.25'):
        comp('col cfl'+cf,R+f'col4x4/cfl/c{cf}_e1e-6',R+f'col4x4/cfl/c{cf}_e1e-9','/output/bsg.out2.*.athdf',ts,o)
    comp('w3d',R+'w3d/e1e-6',R+'w3d/e1e-9','/output/bsg.out1.*.athdf',(0.025,),o)
