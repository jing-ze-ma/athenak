from lib import *
out=open('results2.txt','w')
def pr(*a):
    print(*a,flush=True); print(*a,file=out,flush=True)
z=np.load('hst_usr.npz'); n=list(z['names']); a=z['a']; t=a[:,0]/86400; gh=lambda k:a[:,n.index(k)]
PER=float(os.environ.get('PER','53')); TE=573.0
Lin=gh('L_in')[-1]; pr('L_in %.4e PER %.1f'%(Lin,PER))
# ---------- time scales
sr=(tp>=TE-3*PER)
Rph=Q['Rph_med'][sr].mean(); L0=gh('L_top')[(t>=TE-3*PER)&(t<TE)].mean()
U=Q['U'][sr].mean(0); pr('Rph %.2f Rsun (min %.1f max %.1f since 300 d) Teff %.0f L_top %.4e'%(Rph/RS,Q['Rph_med'][tp>300].min()/RS,Q['Rph_med'][tp>300].max()/RS,Q['Teff_med'][sr].mean(),L0))
vr=Q['vrrms'][sr].mean(0); vt=Q['vtrms'][sr].mean(0); rho=Q['rho'][sr].mean(0); T=Q['Tgas'][sr].mean(0)
Lrc=Q['Lrcom'][sr].mean(0);Lcv=Q['Lcv'][sr].mean(0);LE=Q['LE'][sr].mean(0)
below=r<=Rph
Hrho=-1/np.gradient(np.log(rho),r)
# P from eint & E: gas gamma ~ 5/3 assumed for the estimate only
Pt=Q['E'][sr].mean(0)/3+(2/3)*Q['eint'][sr].mean(0); HP=-1/np.gradient(np.log(Pt),r)
tth=np.array([(U[(r>=x)&below]*dV[(r>=x)&below]).sum()/L0 for x in r])/86400
# radial sound crossing + convective crossing
cs=np.sqrt(5/3*Pt/rho)
tcr=np.array([ (np.gradient(r)[(r>=x)&below]/vr[(r>=x)&below]).sum() for x in r])/86400
tol=HP/vr/86400
np.savez('timescales.npz',r=r,tth=tth,tcr=tcr,tol=tol,HP=HP,vr=vr,vt=vt,cs=cs,rho=rho,T=T,Lrc=Lrc,Lcv=Lcv,LE=LE,U=U,dV=dV,Rph=Rph)
pr('R/Rsun r/Rph m_frac rho T vr[km/s] vt Mach HP/Rs tth[d] tH=HP/vr[d] Lrcom/L Lcv/L')
m=np.cumsum(rho*dV); mt_=m[below][-1]
RR=[8.2,9,10,12,15,20,25,30,40,50,60,70,80,90,95,100]
ii=[int(np.argmin(abs(r/RS-x))) for x in RR]
for i in ii: pr('%6.1f %5.3f %6.4f %8.2e %8.2e %6.2f %6.2f %5.2f %6.2f %7.1f %7.2f %6.3f %6.3f'%(r[i]/RS,r[i]/Rph,m[i]/mt_,rho[i],T[i],vr[i]/1e5,vt[i]/1e5,vr[i]/cs[i],HP[i]/RS,tth[i],tol[i],Lrc[i]/L0,Lcv[i]/L0))
cz=(Lcv/L0>0.1)&below; icz=np.where(cz)[0]; pr('CZ (Lcv/L>0.1) from %.1f to %.1f Rsun, n=%d'%(r[icz[0]]/RS,r[icz[-1]]/RS,len(icz)))
tcz=(np.gradient(r)[cz]/vr[cz]).sum()/86400; pr('convective crossing int dr/vr over CZ = %.1f d, turnover (round trip) %.1f d'%(tcz,2*tcz))
pr('elapsed since t=200 d (end of collapse transient) %.0f d = %.1f turnovers; since t=325 %.0f d'%(TE-200,(TE-200)/(2*tcz),TE-325))
pr('t_th whole star below Rph %.1f d; from 12 Rs %.1f; from 20 Rs %.1f; from 40 Rs %.1f; from 60 %.1f; from 80 %.1f'%(tth[0],tth[ii[3]],tth[ii[5]],tth[ii[8]],tth[ii[10]],tth[ii[12]]))
# ---------- hst trends on period windows
def windows(t0,tend,offs=0):
    e=tend-offs; nw=int((e-t0)//PER); return [(e-(k+1)*PER,e-k*PER) for k in range(nw)][::-1]
pr('\nHST + DUMP TRENDS on %.0f-d window means (OLS on window means; SE from residuals). range [T0, 573]'%PER)
ser={'L_top':gh('L_top'),'L_mid':gh('L_mid'),'L_int':gh('L_int'),'L_bot':gh('L_bot'),'Ltop/Lin':gh('L_top')/gh('L_in'),'Lmid/Lin':gh('L_mid')/gh('L_in'),'Lmid/Ltop':gh('L_mid')/gh('L_top'),'Etot':gh('Etot'),'E_rad':gh('E_rad'),'e_gas':gh('e_gas'),'KE':gh('KE_int'),'KEr':gh('KEr_int'),'M_tot':gh('M_tot'),'M_int':gh('M_int')}
dser={'Rph':Q['Rph_med'],'Teff':Q['Teff_med']}
for T0 in (328.,202.):
    W=windows(T0,TE); xc=np.array([0.5*(u+v) for u,v in W])
    pr('--- T0 %.0f, %d windows'%(T0,len(W)))
    pr('%9s %11s %9s %9s %9s %8s %8s %9s'%('qty','mean','%/100d','SE%/100d','%noise','|s|/SE','t10[d]','t10lo[d]'))
    for nm,y in list(ser.items())+list(dser.items()):
        tt=tp if nm in dser else t
        ym=np.array([y[(tt>=u)&(tt<v)].mean() for u,v in W]); s,se,noise,me=fit(ym,xc)
        fs=100*s*100/me; fse=100*se*100/me
        t10=0.1/abs(s/me); t10lo=0.1/(abs(s/me)+2*se/me)
        pr('%9s %11.4e %9.3f %9.3f %9.3f %8.1f %8.0f %9.0f'%(nm,me,fs,fse,100*noise/me,abs(s)/se if se>0 else 0,t10,t10lo))
    # phase-offset sensitivity
    for nm in ('L_top','Etot','Rph'):
        ss=[]
        for of in (0,9,18,27,36):
            Wo=windows(T0,TE,of);xo=np.array([0.5*(u+v) for u,v in Wo]); y=ser.get(nm,dser.get(nm)); tt=tp if nm in dser else t
            ym=np.array([y[(tt>=u)&(tt<v)].mean() for u,v in Wo]); s=fit(ym,xo)[0]; ss.append(100*100*s/ym.mean())
        pr('  phase-offset spread %s %%/100d: %s'%(nm,np.round(ss,3)))
# energy balance
pr('\nGLOBAL ENERGY BALANCE (hst): DeltaEtot vs int(L_in-L_top+P_esrc-P_sponge)dt, over whole-period spans')
F=np.load('hst_usr_full.npy'); tf=F[:,0]/86400; gf=lambda k:F[:,n.index(k)]
for a0,b0 in ((328,573),(377,573),(475,573),(328,377),(377,426),(426,475),(475,524),(524,573)):
    s=(tf>=a0)&(tf<=b0); tt=F[s,0]; dE=gf('Etot')[s][-1]-gf('Etot')[s][0]; I=np.trapezoid(gf('L_in')[s]-gf('L_top')[s]+gf('P_esrc')[s]-gf('P_sponge')[s],tt); dtt=tt[-1]-tt[0]
    pr(' %3d-%3d dE/dt %+.3e  L_in-L_top-.. %+.3e  res %+.3e (%.2f%% of L) mean L_top %.4e; L_in-L_top/L %.3f%%; Mdot*Phi_top %.2e'%(a0,b0,dE/dtt,I/dtt,(dE-I)/dtt,100*(dE-I)/dtt/Lin,gf('L_top')[s].mean(),100*(Lin-gf('L_top')[s].mean())/Lin,gf('Mdot_top')[s].mean()*Phi[-1]))
pr('Etot t10 irrelevant: its offset is the potential. gain rate/U: t(10%% U) = %.0f d'%(0.1*(gh('E_rad')[-1]+gh('e_gas')[-1])/(Lin-gh('L_top')[t>328].mean())/86400))
pr('M_tot pct change 328->573: %.2e %%; int Mdot_top dt %.2e g'%(100*(gh('M_tot')[-1]/gh('M_tot')[t>328][0]-1),np.trapezoid(gf('Mdot_top')[tf>328],F[tf>328,0])))
