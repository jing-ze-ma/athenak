from lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.25,'grid.linewidth':.6,'legend.frameon':False,'legend.fontsize':8.5,'lines.linewidth':1.3})
CAT=['#0072B2','#D55E00','#009E73','#7A5195','#555555']
SEQ=plt.cm.Blues(np.linspace(.3,1,7))
PER=53.;TE=573.; xc=np.load('winmeans.npz')['xc']; W=np.load('winmeans.npz'); ts=np.load('timescales.npz'); Lg=np.load('lagr.npz'); tth=ts['tth']; Rph=float(ts['Rph'])
z=np.load('hst_usr.npz'); n=list(z['names']); a=z['a']; t=a[:,0]/86400; g=lambda k:a[:,n.index(k)]; Lin=g('L_in')[-1]
W7=[(TE-(k+1)*PER,TE-k*PER) for k in range(7)][::-1]
wm=lambda y,T:np.array([y[(T>=u)&(T<v)].mean() for u,v in W7])
# ---- Fig 1 time series
fig,ax=plt.subplots(5,1,figsize=(10,12),sharex=True); s=t>180; sp=tp>180
def panel(x,T,y,lab,yl,c=CAT[0]):
    x.plot(T,y,color=c,lw=.8,alpha=.65,label=lab); x.plot(xc,wm(y,T),'o-',color=c,ms=6,lw=1.6,label=lab+', 53-d means'); x.set_ylabel(yl)
panel(ax[0],t[s],g('L_top')[s]/1e38,'L_top','L [1e38 erg/s]'); ax[0].axhline(Lin/1e38,color=CAT[4],ls='--',lw=1,label='L_in'); ax[0].plot(t[s],g('L_bot')[s]/1e38,color=CAT[1],lw=1,label='L_bot (lab frame, wall)'); ax[0].set_ylim(1.15,1.65)
ax[0].annotate('strip 315 d (r>140 Rsun only)',(318,1.17),fontsize=8.5,color='#333'); ax[0].axvline(315,color=CAT[4],lw=.7,ls=':')
panel(ax[1],t[s],(g('E_rad')[s]+g('e_gas')[s])/1e45,'E_rad+e_gas','1e45 erg')
x=t>328; ax[2].plot(t[x],g('Etot')[x]/1e45,color=CAT[0],label='Etot after the strip'); I=np.cumsum((g('L_in')[x]-g('L_top')[x])*np.gradient(a[x,0])); ax[2].plot(t[x],(g('Etot')[x][0]+I)/1e45,color=CAT[1],ls='--',label='Etot(328) + int(L_in-L_top)dt'); ax[2].set_ylabel('1e45 erg')
panel(ax[3],t[s],g('KE_int')[s]/1e44,'KE','1e44 erg'); ax[3].plot(t[s],g('KEr_int')[s]/1e44,color=CAT[1],lw=.8,label='KE radial')
panel(ax[4],tp[sp],Q['Rph_med'][sp]/RS,'R_ph median','Rsun'); ax[4].set_xlabel('t [d]')
for x_ in ax: x_.legend(ncol=4,loc='upper right')
fig.suptitle('HG prodG: global series since 180 d, de-duplicated hst',y=.995); plt.tight_layout(); plt.savefig('fig1_series.png',dpi=95); plt.close()
# ---- Fig 2 shell profiles relative to last window
rs=r/RS; m=(rs<Rph/RS*.97)
fig,ax=plt.subplots(2,1,figsize=(10,8),sharex=True)
for k,(key,lab) in enumerate((('rho','rho'),('Tgas','T_gas'))):
    ref=W[key][-1]
    for w in range(6): ax[k].plot(rs[m],100*(W[key][w][m]/ref[m]-1),color=SEQ[w],label='%.0f-%.0f d'%W7[w])
    ax[k].axhline(0,color='#333',lw=.8); ax[k].set_ylabel('%s: % relative to 520-573 d'%lab if False else lab+': percent relative to 520-573 d'); ax[k].legend(ncol=3,loc='upper left'); ax[k].set_ylim(-10,20 if k==0 else 12); ax[k].set_xscale('log')
ax[1].set_xlabel('r [Rsun], below the photosphere'); ax[1].set_xticks([8,10,15,20,30,40,60,80,100]); ax[1].set_xticklabels(['8','10','15','20','30','40','60','80','100'])
fig.suptitle('Successive 53-d shell means: the envelope contraction is dying out'); plt.tight_layout(); plt.savefig('fig2_profiles.png',dpi=95); plt.close()
# ---- Fig 3 fluxes
fig,ax=plt.subplots(1,2,figsize=(12,5),sharey=False)
L=1.38e38; Lr=Lg['Lrcom'][-4:].mean(0); Lc=Lg['Lcv'][-4:].mean(0); B=Lg['Lst'][-4:].mean(0); A=Lg['Lcheck']
mm=(rs<95)
ax[0].plot(rs[mm],Lr[mm]/1e38,color=CAT[0],label='L_rad (comoving)'); ax[0].plot(rs[mm],Lc[mm]/1e38,color=CAT[1],label='L_conv (enthalpy + KE)'); ax[0].plot(rs[mm],B[mm]/1e38,color=CAT[2],label='sum from dumps'); ax[0].plot(rs[mm],A[mm]/1e38,color=CAT[3],ls='--',label='budget flux: L_in - dE/dt - mean flow')
ax[0].axhline(Lin/1e38,color=CAT[4],ls=':',label='L_in'); ax[0].set_xscale('log');ax[0].set_ylim(-.6,2.0);ax[0].set_xticks([8,10,15,20,30,40,60,80]);ax[0].set_xticklabels(['8','10','15','20','30','40','60','80']);ax[0].set_xlabel('r [Rsun]');ax[0].set_ylabel('1e38 erg/s');ax[0].set_title('Flux through the star, 361-573 d');ax[0].legend()
S=Lg['Stor']; ax[1].errorbar(S[:,0],100*(Lin-S[:,2])/Lin,yerr=100*S[:,3]/Lin,fmt='o-',color=CAT[0],ms=6,capsize=3,label='361-573 d'); ax[1].errorbar(S[:,0],100*(Lin-S[:,4])/Lin,yerr=100*S[:,5]/Lin,fmt='s-',color=CAT[1],ms=6,capsize=3,label='202-573 d'); ax[1].axhline(0,color='#333',lw=.8);ax[1].set_xlabel('m / M(below R_ph), from the wall');ax[1].set_ylabel('heat stored below m, percent of L_in');ax[1].set_title('Lagrangian storage (L_in - L(m)) / L_in');ax[1].legend()
plt.tight_layout();plt.savefig('fig3_flux.png',dpi=95);plt.close()
# ---- Fig 4 time scales + t10
fig,ax=plt.subplots(figsize=(10,5.5)); mm=(rs<Rph/RS*.95)
ax.loglog(rs[mm],tth[mm],color=CAT[0],label='thermal time above r: int U dV / L')
tcr=ts['tcr']; ax.loglog(rs[mm],ts['tol'][mm],color=CAT[1],label='H_P / v_r,rms')
ax.axhline(66,color=CAT[1],ls=':',label='CZ round trip 66 d'); ax.axhline(373,color=CAT[4],ls='--',label='elapsed since 200 d: 373 d'); ax.axhline(85,color=CAT[2],ls='-.',label='fitted contraction e-fold ~85 d')
ax.set_xlabel('r [Rsun]');ax.set_ylabel('days');ax.set_xticks([8,10,15,20,30,40,60,80]);ax.set_xticklabels(['8','10','15','20','30','40','60','80']);ax.legend();ax.set_title('Time scales');plt.tight_layout();plt.savefig('fig4_timescales.png',dpi=95);plt.close()
# ---- extra numbers
def ex(x,qi,A,tau): return qi+A*np.exp(-(x-xc[0])/tau)
for nm,Y,T in (('Rph',Q['Rph_med'],tp),('Teff',Q['Teff_med'],tp),('E_rad',g('E_rad'),t),('e_gas',g('e_gas'),t),('L_top',g('L_top'),t)):
    y=wm(Y,T); y=100*np.log(y/y[-1]); p,cv=curve_fit(ex,xc,y,p0=(0,y[0],100),bounds=([-50,-100,10],[50,100,5000]),maxfev=5000); e=y-ex(xc,*p); print(nm,'exp tau %.0f SE %.0f A %+.2f%% rem %+.2f%% noise %.2f%%; windows'%(p[2],np.sqrt(cv[2,2]),p[1],y[-1]-p[0],e.std(ddof=3)),np.round(wm(Y,T)/wm(Y,T)[-1],4))
print('B per window (last 4) at r:')
for R in (17,20,25,30,40,50,60,70): i=np.argmin(abs(rs-R)); print(R,np.round(Lg['Lst'][-7:,i]/1e38,3))
print('L_top windows',np.round(wm(g('L_top'),t)/1e38,3),'L_mid',np.round(wm(g('L_mid'),t)/1e38,3))
