import sys,os,json,glob,re
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/r136a1_1104/page/pipe')
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/r136a1_1104/page')
os.environ.setdefault('MPLBACKEND','Agg')
import numpy as np
from scipy import ndimage as ndi
import bin_convert_repro as bc
R='/resnick/groups/carnegie_poc/jingze/r136a1_1104'
RUN=R+'/runB/prod'
RSUN,DAY,C,GN=6.957e10,86400.,2.99792458e10,6.674e-8
GM0=3.1017784815e28; DF,TF,PF=1e-17,16500.,1e-8; SPD=2e-13
class Table:
    def __init__(self,fn):
        hdr=[l for l in open(fn) if l.startswith('#')]
        g=[l for l in hdr if len(l.split())==7 and l.split()[1].isdigit()][0].split()[1:]
        self.nT,self.nD=int(g[0]),int(g[1]); self.t0,self.dt,self.d0,self.dd=map(float,g[2:])
        self.v=np.loadtxt(fn,comments='#',dtype=np.float32).reshape(self.nT,self.nD)
    def __call__(self,rho,T):
        x=np.clip((np.log10(T)-self.t0)/self.dt,0,self.nT-1.001);y=np.clip((np.log10(rho)-self.d0)/self.dd,0,self.nD-1.001)
        i=x.astype(np.int32);j=y.astype(np.int32);fx=(x-i).astype(np.float32);fy=(y-j).astype(np.float32);v=self.v
        return 10.0**((1-fx)*(1-fy)*v[i,j]+fx*(1-fy)*v[i+1,j]+(1-fx)*fy*v[i,j+1]+fx*fy*v[i+1,j+1])
from r136_page import EOS, assemble, EOSF
eos=EOS(EOSF)
kap=Table(R+'/opac/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt')
rows=[l.split() for l in open(RUN+'/r136.x1grid.txt') if l.strip() and not l.startswith('#')]
i_=np.array([int(x[0]) for x in rows]);x1f=np.array([float(x[1]) for x in rows]);x1f=x1f[(i_>=3)&(i_<=834)]
r=.5*(x1f[1:]+x1f[:-1]);dr=np.diff(x1f);dvr=(x1f[1:]**3-x1f[:-1]**3)/3
x2f=np.linspace(np.pi/3,2*np.pi/3,129);x3f=np.linspace(0,np.pi/3,129);thv=.5*(x2f[1:]+x2f[:-1])
dom=np.diff(x3f)[:,None]*(np.cos(x2f[:-1])-np.cos(x2f[1:]))[None,:]
vol=dom[:,:,None]*dvr[None,None,:]
dth=x2f[1]-x2f[0];dph=x3f[1]-x3f[0]
BANDS=[(35,50),(50,70),(70,90),(90,100),(35,100)]
def load(n):
    dh=bc.read_binary(RUN+'/bin/r136.hydro_w.%05d.bin'%n);t=dh['time'];w=assemble(dh);del dh
    dm=bc.read_binary(RUN+'/bin/r136.m1.%05d.bin'%n);m=assemble(dm);del dm
    D=None
    for f in glob.glob(RUN+'/bin/r136.m1_vet.*.bin'):
        if int(f.split('.')[-2])*2==n:
            dv=bc.read_binary(f)
            if abs(dv['time']-t)<1: D=assemble(dv)
            del dv
    return t,w,m,D
def dtmap(rho,v1,v2,v3,Pg,E,kr,D):
    cg2=(5./3.)*Pg.astype(np.float64)/rho
    h=[dr[None,None,:],(r*dth)[None,None,:],r[None,None,:]*np.sin(thv)[None,:,None]*dph]
    kt=kr.astype(np.float64)*rho
    out=[];cs=[]
    for a,(vv,hh) in enumerate(zip((v1,v2,v3),h)):
        dd=D['m1_d%d%d'%(a+1,a+1)].astype(np.float64) if D is not None else 1./3.
        c2=cg2+(4./3.)*np.maximum(dd,0)*E/rho*(1-np.exp(-kt*hh))
        c=np.sqrt(c2); cs.append(c); out.append(0.3*hh/(np.abs(vv)+c))
    return out,cs,np.sqrt(cg2)
res={}
for n in [int(a) for a in sys.argv[1:]]:
    t,w,m,D=load(n)
    rho=w['dens'].astype(np.float64);v1=w['velx'].astype(np.float64);v2=w['vely'].astype(np.float64);v3=w['velz'].astype(np.float64)
    T,Pg=eos.temp_pgas(w['dens'],w['eint']);T=T.astype(np.float64);Pg=Pg.astype(np.float64)
    E=m['m1_e'].astype(np.float64);Fl=m['m1_f1'].astype(np.float64)
    kr=kap(w['dens'],T.astype(np.float32)).astype(np.float64)
    Prr=(D['m1_d11'].astype(np.float64)*E) if D is not None else E/3
    F0=Fl-v1*(E+Prr)
    dmass=(rho*vol).sum(axis=(0,1));M=GM0/GN+6*(np.cumsum(dmass)-.5*dmass)*2  # FAC=4pi/(pi/3)=12
    g=GN*M/r**2
    gam=kr*F0/(C*g[None,None,:])
    vt=np.sqrt(v2**2+v3**2);ke=.5*rho*(v1**2+v2**2+v3**2)*vol;ker=.5*rho*v1**2*vol
    cg=np.sqrt((5./3.)*Pg/rho)
    Mtot=(rho*vol).sum();KEtot=ke.sum()
    dt_,cs_,_=dtmap(rho,v1,v2,v3,Pg,E,kr,D)
    tm=np.minimum(np.minimum(dt_[0],dt_[1]),dt_[2])
    dtau=kr*rho*dr[None,None,:]
    tau=np.flip(np.cumsum(np.flip(dtau,axis=2),axis=2),axis=2)
    R_={'t_d':t/DAY,'Mtot':Mtot,'KEtot':KEtot,'hasD':D is not None}
    print('\n=== dump %d t=%.3f d  Mtot(wedge)=%.3e KEtot=%.3e vet=%s'%(n,t/DAY,Mtot,KEtot,D is not None))
    # per band
    for a,b in BANDS:
        s=(r>=a*RSUN)&(r<b*RSUN)
        ro=rho[...,s];v=v1[...,s];vtt=vt[...,s];TT=T[...,s];PP=Pg[...,s];cc=cg[...,s];G=gam[...,s];kk=kr[...,s];vv=vol[...,s];
        lo=ro<SPD
        f=lambda q:(100*np.mean(q))
        i=np.unravel_index(np.argmax(np.abs(v)),v.shape);it=np.unravel_index(np.argmax(vtt),v.shape)
        Mb=(ro*vv).sum();KEb=ke[...,s].sum();KErb=ker[...,s].sum();rw=np.sqrt((ro*v**2*vv).sum()/Mb);rwt=np.sqrt((ro*vtt**2*vv).sum()/Mb)
        hi=~lo
        print(f"r{a}-{b}: rho med {np.median(ro):.2e} p1 {np.percentile(ro,1):.2e} p99 {np.percentile(ro,99):.2e} max {ro.max():.2e}| %floor rho<=1.05df {f(ro<=1.05*DF):.2f} <1e-15 N={int((ro<1e-15).sum())} <1e-14 {f(ro<1e-14):.1f}% T<=1.05Tf {f(TT<=1.05*TF):.2f} P<=1.05pf {f(PP<=1.05*PF):.3f} %rho<2e-13 {f(lo):.1f}")
        print(f"   T med {np.median(TT):.3g} max {TT.max():.3g} | max|vr| {abs(v[i])/1e5:.0f} km/s (sign {np.sign(v[i]):.0f}, rho {ro[i]:.1e}, T {TT[i]:.2g}, r {r[s][i[2]]/RSUN:.1f}) vr min {v.min()/1e5:.0f} max {v.max()/1e5:.0f}; %|vr|>400 {f(np.abs(v)>4e7):.2f}, >100 {f(np.abs(v)>1e7):.2f}; max vt {vtt[it]/1e5:.0f} (rho {ro[it]:.1e}); Mach_r max {np.max(np.abs(v)/cc):.0f} p99 {np.percentile(np.abs(v)/cc,99):.1f}; cg med {np.median(cc)/1e5:.0f}")
        print(f"   rho-wt rms vr {rw/1e5:.0f} vt {rwt/1e5:.0f}; Mass frac {Mb/Mtot:.2e} (M={Mb*12/1.989e33:.2e} Msun 4pi) KE/KEtot {KEb/KEtot:.3e} KEr/KE {KErb/max(KEb,1e-99):.2f} | lowrho: |vr| med {np.median(np.abs(v[lo]))/1e5 if lo.any() else 0:.2f} p99 {np.percentile(np.abs(v[lo]),99)/1e5 if lo.any() else 0:.1f} vt med {np.median(vtt[lo])/1e5 if lo.any() else 0:.1f} p99 {np.percentile(vtt[lo],99)/1e5 if lo.any() else 0:.0f}; hirho: |vr| med {np.median(np.abs(v[hi]))/1e5 if hi.any() else 0:.1f} p99 {np.percentile(np.abs(v[hi]),99)/1e5 if hi.any() else 0:.0f} massfrac_hi {(ro[hi]*vv[hi]).sum()/Mb:.3f}")
        Mdot=(ro*v*vv/dr[s][None,None,:]).sum()*12/s.sum() # shell-mean 4pi r^2 rho v
        print(f"   Mdot_band {Mdot/1.989e33*3.156e7:.2e} Msun/yr net; down-part {((ro*np.minimum(v,0)*vv/dr[s][None,None,:]).sum()*12/s.sum())/1.989e33*3.156e7:.2e} up {((ro*np.maximum(v,0)*vv/dr[s][None,None,:]).sum()*12/s.sum())/1.989e33*3.156e7:.2e}; kappa med {np.median(kk):.3f} max {kk.max():.2f}; Gamma(F0) med {np.median(G):.3f} p1 {np.percentile(G,1):.2f} p99 {np.percentile(G,99):.2f}; Gamma hi-rho med {np.median(G[hi]) if hi.any() else 0:.2f}; dtmin {tm[...,s].min():.1f}s")
        # tau & dt
    # threshold test r 35-90
    s=(r>=35*RSUN)&(r<90*RSUN)
    edges=[1e-17,1e-16,1e-15,1e-14,5e-14,1e-13,1.5e-13,1.9e-13,2.1e-13,3e-13,5e-13,1e-12,1e-11,1e-9]
    print(' rho-bin: N, med|vr|, p99|vr|, med vt, med T (35-90 Rsun)')
    ro=rho[...,s]
    for lo_,hi_ in zip(edges[:-1],edges[1:]):
        q=(ro>=lo_)&(ro<hi_)
        if q.sum()==0: print(f'  {lo_:.1e}-{hi_:.1e}: 0');continue
        print(f'  {lo_:.1e}-{hi_:.1e}: {int(q.sum()):8d} {np.median(np.abs(v1[...,s][q]))/1e5:8.2f} {np.percentile(np.abs(v1[...,s][q]),99)/1e5:8.1f} {np.median(vt[...,s][q])/1e5:7.1f} {np.median(T[...,s][q]):9.3g}')
    # pockets
    pk=rho<1e-15
    if pk.any():
        kk_,jj,ii=np.where(pk)
        print(f' pockets rho<1e-15 N={pk.sum()} r range {r[ii].min()/RSUN:.1f}-{r[ii].max()/RSUN:.1f} med {np.median(r[ii])/RSUN:.1f}; T med {np.median(T[pk]):.3g}, E/ (aT^4) med {np.median(E[pk]/(7.5657e-15*T[pk]**4)):.3g}; |vr| med {np.median(np.abs(v1[pk]))/1e5:.1f} vt med {np.median(vt[pk])/1e5:.1f} eg/E med {np.median(w["eint"][pk]/E[pk]):.3e}')
    # hot cell at low rho
    hot=(r[None,None,:]>27*RSUN)&(T>3e5)&(rho<1e-11)
    print(' hot cells T>3e5 & rho<1e-11 r>27: N=',int(hot.sum()),(' rmed %.1f rho med %.2e T max %.3g'%(np.median(np.where(hot)[2]*0+r[np.where(hot)[2]])/RSUN,np.median(rho[hot]),T[hot].max()) if hot.any() else ''))
    u=r[None,None,:]>27*RSUN
    for lim in (300e5,1000e5,3000e5,0.99e9):
        q=u&(vt>lim)
        print(f' vt>{lim/1e5:.0f} km/s: N={int(q.sum())}'+(f' r {r[np.where(q)[2]].min()/RSUN:.1f}-{r[np.where(q)[2]].max()/RSUN:.1f} rho med {np.median(rho[q]):.2e} max {rho[q].max():.1e} T med {np.median(T[q]):.3g} max {T[q].max():.3g}' if q.any() else ''))
    ti=np.argmax(np.where(u,vt,0));K,J,I=np.unravel_index(ti,vt.shape)
    print(f' max vt {vt[K,J,I]/1e5:.0f} at r={r[I]/RSUN:.2f} j={J} k={K} rho={rho[K,J,I]:.2e} T={T[K,J,I]:.3g} v2={v2[K,J,I]/1e5:.0f} v3={v3[K,J,I]/1e5:.0f} vr={v1[K,J,I]/1e5:.0f}; dt2 there {dt_[1][K,J,I]:.2f} dt3 {dt_[2][K,J,I]:.2f}')
    tlat=np.minimum(dt_[1],dt_[2]);K,J,I=np.unravel_index(np.argmin(np.where(rho<1e-13,tm,1e9)),tm.shape)
    print(f' min dt in rho<1e-13: {tm[K,J,I]:.2f} r={r[I]/RSUN:.2f} rho={rho[K,J,I]:.2e} T={T[K,J,I]:.3g} vr={v1[K,J,I]/1e5:.0f} v2={v2[K,J,I]/1e5:.0f} v3={v3[K,J,I]/1e5:.0f} cg={cg[K,J,I]/1e5:.0f}')
    # dt
    tmin=[x.min() for x in dt_]
    print(' CFL dt dir mins',' '.join('%.2f'%x for x in tmin))
    top=np.argsort(tm,axis=None)[:200]
    kt,jt,it=np.unravel_index(top,tm.shape)
    print(' top200 cells r quantiles',np.round(np.percentile(r[it]/RSUN,[0,25,50,75,100]),2))
    for q in range(4):
        K,J,I=kt[q],jt[q],it[q];di=int(np.argmin([dt_[a][K,J,I] for a in range(3)]))
        print(f'  dt {tm[K,J,I]:.2f} dir{di+1} r={r[I]/RSUN:.2f} j={J} k={K} i={I} dr={dr[I]/RSUN:.4f} rho={rho[K,J,I]:.2e} T={T[K,J,I]:.3g} E={E[K,J,I]:.3g} vr={v1[K,J,I]/1e5:.0f} vt={vt[K,J,I]/1e5:.0f} cg={cg[K,J,I]/1e5:.0f} c{di+1}={cs_[di][K,J,I]/1e5:.0f} tau={tau[K,J,I]:.2f}')
    print(' min dt by r band (with limiting speeds):')
    for a,b in [(19.45,24),(24,26),(26,28),(28,30),(30,32),(32,35),(35,40),(40,50),(50,70),(70,100)]:
        s=(r>=a*RSUN)&(r<b*RSUN);sub=tm[...,s];K,J,I=np.unravel_index(np.argmin(sub),sub.shape);I2=np.where(s)[0][I]
        print(f'   {a}-{b}: dt {sub.min():6.2f} p1 {np.percentile(sub,1):6.2f} at r={r[I2]/RSUN:.2f} dr={dr[I2]/RSUN:.3f} rho={rho[K,J,I2]:.1e} T={T[K,J,I2]:.3g} |vr|={abs(v1[K,J,I2])/1e5:.0f} c1={cs_[0][K,J,I2]/1e5:.0f} cg={cg[K,J,I2]/1e5:.0f}')
    # tau
    i35=np.argmin(abs(r-35*RSUN));i50=np.argmin(abs(r-50*RSUN))
    print(f' tau(35) med {np.median(tau[...,i35]):.3f} max {tau[...,i35].max():.2f} p99 {np.percentile(tau[...,i35],99):.2f}; tau(50) med {np.median(tau[...,i50]):.3f} max {tau[...,i50].max():.2f}; columns tau35>1: {100*np.mean(tau[...,i35]>1):.1f}%')
    rph=[];
    ii=(tau>=1).sum(axis=2)
    Rphc=r[np.clip(ii-1,0,len(r)-1)]
    print(f' R_ph columns: p1 {np.percentile(Rphc,1)/RSUN:.1f} med {np.median(Rphc)/RSUN:.1f} p99 {np.percentile(Rphc,99)/RSUN:.1f} max {Rphc.max()/RSUN:.1f}')
    # dt excluding rho<1e-12 tenuous gas
    print(' dt min if only tau>1 cells: %.2f ; only tau<1: %.2f ; tau<1 & r>35: %.2f'%(tm[tau>1].min(),tm[tau<=1].min(),tm[(tau<=1)&(r[None,None,:]>35*RSUN)].min()))
    # blockiness: fast mask comps on shells 35-100
    s=np.where((r>=35*RSUN)&(r<100*RSUN))[0]
    fast=(v1[...,s]<-1e7)
    lab,nl=ndi.label(fast,structure=np.ones((3,3,3)))
    print(' fast(vr<-100) N',int(fast.sum()),'comps',nl)
    if nl:
        objs=ndi.find_objects(lab);sz=np.bincount(lab.ravel())[1:]
        wk=np.array([o[0].stop-o[0].start for o in objs]);wj=np.array([o[1].stop-o[1].start for o in objs]);wi=np.array([o[2].stop-o[2].start for o in objs])
        big=np.argsort(-sz)[:5]
        print('  comp size pct',np.percentile(sz,[50,90,99]),'max',sz.max(),'frac in single-cell comps %.3f'%(sz[sz==1].sum()/sz.sum()))
        print('  biggest 5 (Ncells, dk,dj,di, r range, dmass frac of all fast):',[(int(sz[b]),int(wk[b]),int(wj[b]),int(wi[b]),round(r[s][objs[b][2].start]/RSUN,1),round(r[s][objs[b][2].stop-1]/RSUN,1)) for b in big])
        print('  median widths dk %d dj %d di %d'%(np.median(wk),np.median(wj),np.median(wi)))
        # overlap with hi-rho
        print('  fast cells with rho>=2e-13 : %.3f ; hi-rho cells that are fast: %.3f'%(np.mean(rho[...,s][fast]>=SPD),np.mean(fast[rho[...,s]>=SPD])))
        # block-edge alignment: mask jumps across j
        dj=np.abs(np.diff(fast.astype(np.int8),axis=1)).sum(axis=(0,2));dk=np.abs(np.diff(fast.astype(np.int8),axis=0)).sum(axis=(1,2))
        print('  mask edges: at j=63|64 %d vs median per j %d ; k=63|64 %d vs median %d'%(dj[63],np.median(dj),dk[63],np.median(dk)))
    # checkerboard: lateral neighbor contrast (rho) vs radial, by band
    for a,b in [(35,50),(50,70),(70,100)]:
        ss=(r>=a*RSUN)&(r<b*RSUN);lr=np.log10(rho[...,ss])
        print(f'  lat/rad contrast r{a}-{b}: rms dlog rho dj {np.sqrt(np.mean(np.diff(lr,axis=1)**2)):.3f} dk {np.sqrt(np.mean(np.diff(lr,axis=0)**2)):.3f} di {np.sqrt(np.mean(np.diff(lr,axis=2)**2)):.3f}; rms dvr/kms dj {np.sqrt(np.mean(np.diff(v1[...,ss],axis=1)**2))/1e5:.1f} di {np.sqrt(np.mean(np.diff(v1[...,ss],axis=2)**2))/1e5:.1f} lag-1 corr vr j {np.corrcoef(v1[...,ss][:,:-1,:].ravel(),v1[...,ss][:,1:,:].ravel())[0,1]:.3f}')
    np.savez_compressed('/resnick/groups/carnegie_poc/jingze/r136noise_1106/red_%05d.npz'%n,t=t,r=r,dr=dr,
      tmr=tm.min(axis=(0,1)),eq_rho=w['dens'][:,64,:],eq_vr=w['velx'][:,64,:],eq_vt=vt[:,64,:].astype(np.float32),
      hist=np.histogram2d(np.log10(rho[...,(r>=35*RSUN)&(r<90*RSUN)]).ravel(),np.log10(np.maximum(np.abs(v1[...,(r>=35*RSUN)&(r<90*RSUN)]),1.)/1e5).ravel(),bins=[np.linspace(-17.5,-8,96),np.linspace(-3,3.5,66)])[0])
    sys.stdout.flush()
    del w,m,D,rho,v1,v2,v3,T,Pg,E,Fl,kr,tm,tau,gam,dt_,cs_,F0
