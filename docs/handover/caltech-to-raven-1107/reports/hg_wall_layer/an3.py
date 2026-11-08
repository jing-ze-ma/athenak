from lib import *
out=open('results3.txt','w')
def pr(*a):
    print(*a,flush=True); print(*a,file=out,flush=True)
PER=53.; TE=573.; ts=np.load('timescales.npz'); tth=ts['tth']; Rph=float(ts['Rph'])
z=np.load('hst_usr.npz'); n=list(z['names']); a=z['a']; t=a[:,0]/86400; Lin=a[-1,n.index('L_in')]; L0=a[(t>=361)&(t<TE),n.index('L_top')].mean()
W7=[(TE-(k+1)*PER,TE-k*PER) for k in range(7)][::-1]; xc7=np.array([.5*(u+v) for u,v in W7])
pr('windows',W7)
RR=[9,10,12,15,20,25,30,40,50,60,70,80,90,95]; ii=[int(np.argmin(abs(r/RS-x))) for x in RR]
wm=lambda q:np.array([q[(tp>=u)&(tp<v)].mean(0) for u,v in W7])
W={k:wm(Q[k]) for k in ('rho','Tgas','E','eint','vrrms','vtrms','LE','Lr','Lrcom','Lconv','Lkin','Lcv','Lbulk','Pm','U','Etotd','Lfl','mvr')}
np.savez('winmeans.npz',xc=xc7,**W)
def table(lab,M):
    pr('\n%s at fixed r: %%/100d: last 4 windows (361-573), all 7 (202-573); noise=%% residual std of a window mean; driftpertth=%% per t_th(r) (7w); t10[d]=10%%/|slope| (7w), t10lo uses |s|+2SE; last column: drift over 371 d / noise'%lab)
    pr('  r/Rs    t_th  s4    SE4    s7    SE7   noise  s/SE7 %/t_th    t10   t10lo  span/noise')
    for R,i in zip(RR,ii):
        y=np.log(M[:,i]); s4,se4,_,_=fit(y[-4:],xc7[-4:]); s7,se7,nz,_=fit(y,xc7)
        pr('%6.1f %6.1f %6.2f %5.2f %6.2f %5.2f %6.2f %5.1f %7.2f %7.0f %6.0f %6.1f'%(R,tth[i],1e4*s4,1e4*se4,1e4*s7,1e4*se7,100*nz,abs(s7)/se7,100*s7*tth[i],0.1/abs(s7),0.1/(abs(s7)+2*se7),abs(s7)*(xc7[-1]-xc7[0])/nz))
for lab,k in (('rho','rho'),('T_gas','Tgas'),('E_rad','E'),('e_gas','eint'),('vr_rms','vrrms'),('vt_rms','vtrms')): table(lab,W[k])
# per-window values at radii
pr('\nwindow means (rows = windows) rho at radii [1e-8 scale], T [1e4 K]')
for k,sc in (('rho',1),('Tgas',1)):
    pr(k,'ratio to last window, in %:'); pr(' '.join('%7.1f'%x for x in RR))
    for w in range(7): pr(' '.join('%7.2f'%(100*(W[k][w,i]/W[k][-1,i]-1)) for i in ii),' <- %.0f d'%xc7[w])
# fixed mass coordinate from the wall; r(m), T(m), rho(m) per dump
m=np.cumsum(Q['rho']*dV,axis=1)
Fr=[.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95]
Mref=m[(tp>=361)][:,r<=Rph][:,-1].mean(); pr('\nFIXED MASS: fractions of mass below Rph (%.3e g)'%Mref)
RM=np.zeros((len(tp),len(Fr))); TM=RM.copy(); DM=RM.copy()
for j in range(len(tp)):
    RM[j]=np.interp(np.array(Fr)*Mref,m[j],r); TM[j]=np.interp(RM[j],r,Q['Tgas'][j]); DM[j]=np.interp(RM[j],r,Q['rho'][j])
pr(' m/M   r/Rs  sR4 SE  sR7 SE noise | sT7 SE noise t10T | srho7 SE noise t10rho  (%/100d)')
for k,f in enumerate(Fr):
    o=[]
    for Y in (RM,TM,DM):
        yy=np.array([np.log(Y[(tp>=u)&(tp<v),k]).mean() for u,v in W7]); o.append((fit(yy[-4:],xc7[-4:]),fit(yy,xc7)))
    (a4,b7),(_,c7),(_,d7)=o
    pr('%.2f %6.1f  %+.2f %.2f  %+.2f %.2f %.2f | %+.2f %.2f %.2f %6.0f | %+.2f %.2f %.2f %6.0f'%(f,RM[:,k][tp>361].mean()/RS,1e4*a4[0],1e4*a4[1],1e4*b7[0],1e4*b7[1],100*b7[2],1e4*c7[0],1e4*c7[1],100*c7[2],.1/abs(c7[0]),1e4*d7[0],1e4*d7[1],100*d7[2],.1/abs(d7[0])))
np.savez('massc.npz',tp=tp,RM=RM,TM=TM,DM=DM,Fr=np.array(Fr))
# zone budget
E_dens=Q['Etotd']; Edges=[8.0,12,20,30,40,60,80,95]
eidx=[np.searchsorted(rf,x*RS) for x in Edges]; eidx[0]=0
pr('\nZONE BUDGET. Zones between radii; dE/dt = conservative energy incl KE + potential; U = eint+E only. LE(a), LE(b) = total energy flux at the zone faces (period-mean); units 1e36 erg/s')
pr(' zone[Rs]     U[1e45]  tth_U[d]  dU/dt4 SE  dE/dt4 SE  dE/dt7 SE  LEa-LEb(4w) SE_mean  ')
Zs=[]
for a_,b_ in zip(eidx[:-1],eidx[1:]):
    Ez=(E_dens[:,a_:b_]*dV[a_:b_]).sum(1); Uz=(Q['U'][:,a_:b_]*dV[a_:b_]).sum(1); Mz=(Q['rho'][:,a_:b_]*dV[a_:b_]).sum(1)
    ye=np.array([Ez[(tp>=u)&(tp<v)].mean() for u,v in W7]); yu=np.array([Uz[(tp>=u)&(tp<v)].mean() for u,v in W7]); ym=np.array([Mz[(tp>=u)&(tp<v)].mean() for u,v in W7])
    s4=fit(ye[-4:],xc7[-4:]*86400); s7=fit(ye,xc7*86400); u4=fit(yu[-4:],xc7[-4:]*86400); u7=fit(yu,xc7*86400); m7=fit(ym,xc7*86400)
    la=W['LE'][-4:,a_].mean(0) if a_>0 else Lin; lb=W['LE'][-4:,b_].mean(0); sdb=W['LE'][-4:,b_].std()/2
    pr('%4.0f-%-4.0f %9.3f %8.0f %+7.1f %5.1f %+7.1f %5.1f %+7.1f %5.1f %+8.1f %5.1f   dM/M per 100d %+.2f%%'%(r[a_]/RS,r[b_]/RS,yu.mean()/1e45,yu.mean()/L0/86400,u4[0]/1e36,u4[1]/1e36,s4[0]/1e36,s4[1]/1e36,s7[0]/1e36,s7[1]/1e36,(la-lb)/1e36,sdb/1e36,100*m7[0]*8.64e6/ym.mean()))
    Zs.append((r[a_]/RS,r[b_]/RS,s4[0],s4[1],la-lb))
# flux table
pr('\nFLUX (1e38) in period-mean windows: radius; L_E, L_rad(com), L_cv(enth+KE), for windows 4 last and their std (period to period)')
pr('  r/Rs   L_E: 361-414 414-467 467-520 520-573 |  std | L_rcom |L_cv | Lr(lab) Lfl | Pm[1e22 g/s] | Lin-Ecum(<r) slope')
ecum=np.cumsum((E_dens*dV),axis=1); ec4=np.array([ecum[(tp>=u)&(tp<v)].mean(0) for u,v in W7[-4:]]); xs=xc7[-4:]*86400
dEcum=np.array([np.polyfit(xs,ec4[:,i],1)[0] for i in range(len(r))])
for i,R in zip(ii,RR):
    pr('%6.1f '%R,' '.join('%6.3f'%(W['LE'][w,i]/1e38) for w in range(3,7)),'| %.3f | %.3f |%.3f | %.3f %.3f | %+.2f | %.3f'%(W['LE'][-4:,i].std()/1e38,W['Lrcom'][-4:,i].mean()/1e38,W['Lcv'][-4:,i].mean()/1e38,W['Lr'][-4:,i].mean()/1e38,W['Lfl'][-4:,i].mean()/1e38,W['Pm'][-4:,i].mean()/1e22,(Lin-dEcum[i])/1e38))
np.savez('flux.npz',r=r,LE=W['LE'],Lrcom=W['Lrcom'],Lcv=W['Lcv'],Lr=W['Lr'],Lfl=W['Lfl'],dEcum=dEcum,Lin=Lin,L0=L0,xc=xc7)
# wall layer
pr('\nWALL: i r/Rs T rho Lr Lrcom Lfl LE (361-573, 1e38), and T, rho s7 %/100d')
for i in list(range(0,8))+[12,16,24]:
    pr(i,'%.2f %.3e %.3e %.3f %.3f %.3f %.3f  sT %+.3f srho %+.3f mvr %.2e'%(r[i]/RS,W['Tgas'][-4:,i].mean(),W['rho'][-4:,i].mean(),W['Lr'][-4:,i].mean()/1e38,W['Lrcom'][-4:,i].mean()/1e38,W['Lfl'][-4:,i].mean()/1e38,W['LE'][-4:,i].mean()/1e38,1e4*fit(np.log(W['Tgas'][:,i]),xc7)[0],1e4*fit(np.log(W['rho'][:,i]),xc7)[0],W['mvr'][-4:,i].mean()))
# pulsation amplitude: sinusoid fit at the period in sliding windows
pr('\nPULSATION amplitude (sinusoid fit P free 45-60 d, 2-period sliding window)')
ti=tp; 
for nm,y,tt in (('Rph',Q['Rph_med'],tp),('L_top',a[:,n.index('L_top')],t)):
    pr(nm)
    for c in np.arange(205,540,53):
        s=(tt>=c-53)&(tt<c+53); best=None
        for Pp in np.linspace(45,60,31):
            A=np.c_[np.cos(2*np.pi*tt[s]/Pp),np.sin(2*np.pi*tt[s]/Pp),np.ones(s.sum()),tt[s]-c]; co,*_=np.linalg.lstsq(A,y[s],rcond=None); e=((y[s]-A@co)**2).sum()
            if best is None or e<best[0]: best=(e,Pp,np.hypot(*co[:2])/co[2])
        pr('  c %.0f P %.1f amp %.2f%% rms-rest %.2f%%'%(c,best[1],100*best[2],100*np.sqrt(best[0]/s.sum())/y[s].mean()))
