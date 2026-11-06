#!/usr/bin/env python3
"""4x4 column: KE history, shell-mean wall profiles, dump stats, comparison with his1d_g."""
import numpy as np, glob, re, subprocess
RUN='/orion/ptmp/jinma/athpp_repro3d/col4x4/run/'
ru,Tu,lu,mu=6.15245e-5,1.07887e6,6.957e10,.62
Rg=8.314462618e7/mu; vu=(Rg*Tu)**.5; tu=lu/vu; Pu=ru*vu*vu; Eu=Pu*lu**3; a=7.565733250033928e-15; c=2.99792458e10
Eru=a*Tu**4; Fu=c*Eru; GM=260.795*lu**3/tu**2; gam=5/3; Rs=lu
def G1(Pg,Pr):
    b=Pg/(Pg+Pr); return b+(4-3*b)**2*(gam-1)/(b+12*(gam-1)*(1-b))
h=np.loadtxt(RUN+'output/bsg.hst'); t=h[:,0]*tu/86400
with open('col4x4_ke.txt','w') as f:
    f.write('# Orion 4x4 Athena++ column (job 227389): t[d] KE KE_r KE_lat [erg, pi/3 sr], hst cadence (0.1 code = 578 s)\n')
    for i in range(len(t)): f.write('%.6e %.6e %.6e %.6e\n'%(t[i],(h[i,6]+h[i,7]+h[i,8])*Eu,h[i,6]*Eu,(h[i,7]+h[i,8])*Eu))
log=open(RUN+'run.log',errors='ignore').read()
cyc=[(float(m[0]),float(m[1])) for m in re.findall(r'cycle=\d+ time=([\d.e+-]+) dt=([\d.e+-]+)',log)]
ct=np.array([x[0] for x in cyc]); cd=np.array([x[1] for x in cyc])
# niter per interval: walk log once
nit=[];pos=[]
for m in re.finditer(r'cycle=\d+ time=([\d.e+-]+)|niter: (\d+)',log):
    if m.group(1): cur=float(m.group(1))
    else: nit.append(int(m.group(2))); pos.append(cur)
nit=np.array(nit);pos=np.array(pos)
fs=sorted(glob.glob(RUN+'output/bsg.block0.out2.*.tab'))
W=open('col4x4_wall_profiles.txt','w'); W.write('# Orion 4x4 column, shell mean over 16 lateral cells (solid-angle weighted), cgs: i r[Rsun] rho Pgas Er(lab) T vr Fr1(lab) Fr01(comoving) N2_face(i,i+1)[s^-2]\n')
S=open('col4x4_dumpstats.txt','w'); S.write('# t[d] dt[s] jac_mean jac_max jac_capped rms_vth[cm/s] rms_vph[cm/s] maxabs_vth[cm/s] maxabs_vph[cm/s] rms_lat(drho/rho)_first24 rms_lat(dT/T)_first24  (lateral stats over cells i<24 for rms_lat; v over all r<60)\n')
prof={}
prev=0.
for n,_ in enumerate(fs):
    tag='%05d'%n; rows=[]
    for g in glob.glob(RUN+f'output/bsg.block*.out2.{tag}.tab'):
        rows.append(np.loadtxt(g))
    if not rows: continue
    d=np.vstack(rows); tc=float(open(RUN+f'output/bsg.block0.out2.{tag}.tab').readline().split('time=')[1].split()[0]); T=tc*tu/86400
    r=d[:,1];th=d[:,3];ph=d[:,5];rs=np.unique(r);ths=np.unique(th);phs=np.unique(ph)
    dth=ths[1]-ths[0]; w=np.sin(th)*dth   # dOmega
    ir=np.searchsorted(rs,r);nr=len(rs)
    sm=lambda q:np.bincount(ir,w*q,nr)/np.bincount(ir,w,nr)
    rho,Pg,v1,v2,v3,Er,F1,Er0,F01=d[:,6]*ru,d[:,7]*Pu,d[:,8]*vu,d[:,9]*vu,d[:,10]*vu,d[:,11]*Eru,d[:,12]*Fu,d[:,16]*Eru,d[:,17]*Fu
    Tg=d[:,7]/d[:,6]*Tu
    mr,mP,mE,mT,mv,mF,mF0=[sm(q) for q in (rho,Pg,Er,Tg,v1,F1,F01)]
    Pr=mE/3;P=mP+Pr;rc=rs*lu;fr=.5*(rc[:-1]+rc[1:]);gf=GM/fr**2;dr=np.diff(rc)
    N=gf*(np.diff(np.log(P))/dr/G1(.5*(mP[:-1]+mP[1:]),.5*(Pr[:-1]+Pr[1:]))-np.diff(np.log(mr))/dr)
    W.write('\n# t=%.4f d\n'%T)
    for i in range(24): W.write('%3d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n'%(i,rs[i],mr[i],mP[i],mE[i],mT[i],mv[i],mF[i],mF0[i],N[i]))
    prof[round(T*2)/2]=(rs,mv,mr,mP,mE,mT,mF,N)
    i24=rs[23]; m24=r<=i24; m60=r<=60.
    lat=lambda q,qm:(np.sqrt(np.sum(w*((q-qm[ir])/qm[ir])**2*m24)/np.sum(w*m24)))
    mask=(pos>prev)&(pos<=tc); j=nit[mask] if mask.any() else np.array([])
    k=np.argmin(abs(ct-tc)) if len(ct) else 0
    S.write('%.4f %.3f %.3f %d %d %.4e %.4e %.4e %.4e %.4e %.4e\n'%(T,cd[k]*tu,j.mean() if len(j) else float('nan'),j.max() if len(j) else -1,int((j>=100).sum()) if len(j) else -1,
        (np.sum(w*v2**2*m60)/np.sum(w*m60))**.5,(np.sum(w*v3**2*m60)/np.sum(w*m60))**.5,abs(v2[m60]).max(),abs(v3[m60]).max(),lat(rho,mr),lat(Tg,mT)))
    prev=tc
W.close();S.close()
# his1d_g comparison
ref=subprocess.run(['git','-C','/orion/u/jinma/ATHENAK/athenak','show','origin/bsg-files-1009:docs/handover/athenapp-bsg/cflrad_1106/his1d_g_ke.txt'],capture_output=True,text=True).stdout
rk=np.array([list(map(float,l.split())) for l in ref.splitlines() if l and l[0]!='#'])
wp=subprocess.run(['git','-C','/orion/u/jinma/ATHENAK/athenak','show','origin/bsg-files-1009:docs/handover/athenapp-bsg/cflrad_1106/his1d_g_wall_profiles.txt'],capture_output=True,text=True).stdout
blk={};cur=None
for l in wp.splitlines():
    if l.startswith('# t='): cur=float(l.split('=')[1].split()[0]);blk[cur]=[]
    elif l and l[0]!='#' and cur is not None: blk[cur].append(list(map(float,l.split())))
blk={k:np.array(v) for k,v in blk.items()}
ke=(h[:,6]+h[:,7]+h[:,8])*Eu
print('t_last %.3f d, dumps %d'%(t[-1],len(fs)))
for lo,hi in ((0,1),(1,2),(2,4),(4,6),(6,8),(8,10),(10,14),(14,18),(18,22),(22,25)):
    m=(t>=lo)&(t<hi);n=(rk[:,0]>=lo)&(rk[:,0]<hi)
    if m.any(): print(f'{lo}-{hi} d: col median {np.median(ke[m]):.2e} max {ke[m].max():.2e}  his1d_g median {np.median(rk[n,1]):.2e} max {rk[n,1].max():.2e}  n={m.sum()}')
print('col peak KE %.3e at %.3f d; his1d_g peak %.3e at %.3f d'%(ke.max(),t[ke.argmax()],rk[:,1].max(),rk[rk[:,1].argmax(),0]))
for T in (0,1,2,3,4):
    if T in prof:
        rs,mv,*_=prof[T]; print('t=%g d vr(r=22.66): col %.1f m/s; '%(T,np.interp(22.66,rs,mv)/100),end='')
        k=min(blk,key=lambda x:abs(x-T));b=blk[k];print('his1d_g %.1f m/s (table covers %.2f..%.2f)'%(np.interp(22.66,b[:,1],b[:,6])/100,b[0,1],b[-1,1]))
        rs2,mv2,mr2,mP2,mE2,mT2,mF2,N2=prof[T];
        print('   i=0: r %.3f rho %.5e/%.5e T %.4e/%.4e Fr1 %.4e/%.4e N2 %.3e/%.3e (col/his1d_g)'%(rs2[0],mr2[0],b[0,2],mT2[0],b[0,5],mF2[0],b[0,7],N2[0],b[0,8]))
