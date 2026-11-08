#!/usr/bin/env python3
"""an.py <run dir> [dump numbers]: per-dump metrics of one R136 arm -> <run dir>/an.json. Uses r136noise_1106/an1.py's readers and dtmap."""
import sys,json,glob,re,os
RD=sys.argv[1]
src=open('/resnick/groups/carnegie_poc/jingze/r136noise_1106/an1.py').read().split('res={}')[0]
sys.argv=[sys.argv[0]]+sys.argv[2:]
exec(src.replace("RUN=R+'/runB/prod'","RUN=%r"%RD))
ns=[int(a) for a in sys.argv[1:]] or sorted(int(f.split('.')[-2]) for f in glob.glob(RD+'/bin/r136.hydro_w.*.bin'))
TTOP=6.31e6*0.98   # EOS table top log T 6.8
out={}
fp='/resnick/groups/carnegie_poc/jingze/r136open_1106/an_'+os.path.basename(RD.rstrip('/'))+'.json'
if os.path.exists(fp): out=json.load(open(fp))
GM=None
for n in ns:
    if str(n) in out: continue
    t,w,m,D=load(n)
    rho=w['dens'].astype(np.float64);v1=w['velx'].astype(np.float64);v2=w['vely'].astype(np.float64);v3=w['velz'].astype(np.float64)
    T,Pg=eos.temp_pgas(w['dens'],w['eint']);T=T.astype(np.float64);Pg=Pg.astype(np.float64)
    E=m['m1_e'].astype(np.float64);kr=kap(w['dens'],T.astype(np.float32)).astype(np.float64)
    vt=np.sqrt(v2**2+v3**2);dt_,cs_,cg=dtmap(rho,v1,v2,v3,Pg,E,kr,D);tm=np.minimum(np.minimum(dt_[0],dt_[1]),dt_[2])
    dtau=kr*rho*dr[None,None,:];tau=np.flip(np.cumsum(np.flip(dtau,axis=2),axis=2),axis=2)
    Rph=r[np.clip((tau>=1).sum(axis=2)-1,0,len(r)-1)]
    q={'t_d':t/DAY,'hasD':D is not None,'Mwedge':float((rho*vol).sum())}
    q['Rph_med'],q['Rph_p1'],q['Rph_p99']=[float(x)/RSUN for x in np.percentile(Rph,[50,1,99])]
    q['Rph_mean']=float(np.sum(Rph*dom)/np.sum(dom)/RSUN)
    abv=r[None,None,:]>Rph[:,:,None]
    q['M_above_Rph']=float((rho*vol)[abv].sum()); q['KE_above_Rph']=float((.5*rho*(v1**2+vt**2)*vol)[abv].sum())
    for a,b in [(35,50),(50,90),(90,100),(35,100),(19.45,24),(24,27),(27,35)]:
        s=(r>=a*RSUN)&(r<b*RSUN);ro=rho[...,s];v=v1[...,s];vv=vol[...,s];vtt=vt[...,s];TT=T[...,s];Mb=(ro*vv).sum();c=cg[...,s]
        k='%g_%g'%(a,b)
        q[k]=dict(M=float(Mb),rho_med=float(np.median(ro)),rho_min=float(ro.min()),rho_p1=float(np.percentile(ro,1)),vrmax=float(np.abs(v).max()/1e5),vrmin=float(v.min()/1e5),vrpos=float(v.max()/1e5),vtmax=float(vtt.max()/1e5),
          vr_rms=float(np.sqrt(np.mean(v**2))/1e5),vt_rms=float(np.sqrt(np.mean(vtt**2))/1e5),vr_rmsw=float(np.sqrt((ro*v**2*vv).sum()/Mb)/1e5),vt_rmsw=float(np.sqrt((ro*vtt**2*vv).sum()/Mb)/1e5),
          f100=float(np.mean(np.abs(v)>1e7)),f400=float(np.mean(np.abs(v)>4e7)),fin100=float(np.mean(v<-1e7)),fflo=float(np.mean(ro<=1.05*DF)),n1e15=int((ro<1e-15).sum()),fT=float(np.mean(TT<=1.05*TF)),T_med=float(np.median(TT)),T_max=float(TT.max()),
          T_mw=float((TT*ro*vv).sum()/Mb),Mach_p99=float(np.percentile(np.abs(v)/c,99)),dtmin=float(tm[...,s].min()),Mdot=float((ro*v*vv/dr[s][None,None,:]).sum()*12/s.sum()),Mdot_in=float((ro*np.minimum(v,0)*vv/dr[s][None,None,:]).sum()*12/s.sum()))
    # top shells: fate of gas at the boundary
    q['top']=dict(rho_med=float(np.median(rho[...,-1])),rho_min=float(rho[...,-1].min()),rho_max=float(rho[...,-1].max()),vr_med=float(np.median(v1[...,-1])/1e5),vr_min=float(v1[...,-1].min()/1e5),vr_max=float(v1[...,-1].max()/1e5),fin=float(np.mean(v1[...,-1]<0)),T_med=float(np.median(T[...,-1])))
    # rho profile (shell median, shell mean), T profile for plot
    q['prof']=dict(r=[float(x)/RSUN for x in r[::8]],rho_med=[float(x) for x in np.median(rho,axis=(0,1))[::8]],rho_mean=[float(x) for x in rho.mean(axis=(0,1))[::8]],T_mean=[float(x) for x in T.mean(axis=(0,1))[::8]],vr_rms=[float(x)/1e5 for x in np.sqrt((v1**2).mean(axis=(0,1)))[::8]],vt_rms=[float(x)/1e5 for x in np.sqrt((vt**2).mean(axis=(0,1)))[::8]],tmr=[float(x) for x in tm.min(axis=(0,1))[::8]])
    # dt setters
    dir_=np.argmin(np.stack(dt_),axis=0);top=np.argsort(tm,axis=None)[:6];kt,jt,it=np.unravel_index(top,tm.shape);L=[]
    for K,J,I in zip(kt,jt,it):
        L.append(dict(dt=float(tm[K,J,I]),dir=int(dir_[K,J,I])+1,r=float(r[I]/RSUN),k=int(K),j=int(J),i=int(I),rho=float(rho[K,J,I]),T=float(T[K,J,I]),vr=float(v1[K,J,I]/1e5),vt=float(vt[K,J,I]/1e5),cg=float(cg[K,J,I]/1e5),tau=float(tau[K,J,I])))
    q['dtcells']=L; q['dtdirmin']=[float(x.min()) for x in dt_]
    u=r[None,None,:]>27*RSUN
    q['ntop_T']=int((T>=TTOP).sum()); q['ntop_T_27']=int(((T>=TTOP)&u).sum()); q['nhot3e5']=int(((T>3e5)&(rho<1e-11)&u).sum()); q['Tmax27']=float(T[...,r>27*RSUN].max())
    for lim,nm in ((300e5,'vt300'),(1000e5,'vt1000'),(3000e5,'vt3000'),(0.99e9,'vtceil')): q['n'+nm]=int(((vt>lim)&u).sum())
    q['nvrceil']=int(((np.abs(v1)>0.99e9)&u).sum()); q['nvrtop']=[int(((np.abs(v1)>l)&u).sum()) for l in (5e7,1e8,3e8)]
    pk=rho<1e-15; q['npock']=int(pk.sum()); q['pock_r']=[float(r[np.where(pk)[2]].min()/RSUN),float(np.median(r[np.where(pk)[2]])/RSUN),float(r[np.where(pk)[2]].max()/RSUN)] if pk.any() else []
    q['floor_all']=int((rho<=1.05*DF).sum())
    # dt percentiles on tau>1 vs tau<1, by region
    q['dt_tau_gt1']=float(tm[tau>1].min()); q['dt_tau_lt1']=float(tm[tau<=1].min()); q['dt_r35']=float(tm[...,r>35*RSUN].min())
    out[str(n)]=q; json.dump(out,open(fp,'w'))
    print(n,'t=%.3f Rph %.2f dtmin %.2f Mabv %.3e'%(q['t_d'],q['Rph_mean'],tm.min(),q['M_above_Rph']),flush=True)
    del w,m,D,rho,v1,v2,v3,T,Pg,E,kr,vt,tm,tau,dt_,cs_,cg
