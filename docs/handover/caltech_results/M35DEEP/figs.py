import sys, glob, re, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/pages_c1_1124/src'); import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, idort_style
idort_style.use(font=12); plt.rcParams.update({'axes.linewidth':1.2,'lines.linewidth':1.4,'figure.dpi':100,'savefig.bbox':None,'savefig.format':'png','xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':9})
RP=9.554e11; z=np.load('ring.npz'); t=z['t']; r=z['r']; R=z['R']/1e5; ir=lambda x:int(np.argmin(abs(r-x)))
fig,ax=plt.subplots(figsize=(7,4.6),constrained_layout=True)
for x,c in ((0.72,'C0'),(0.74,'k'),(0.76,'r')): ax.semilogy(t/86400,R[:,ir(x)],c,label='r = %.2f R_ph,IC'%x)
ax.axvline(1.8967e5/86400,color='0.5',ls=':',lw=1); ax.set_ylim(1e-6,5); ax.set_xlabel('t [d]'); ax.set_ylabel(r'ring rms $v_r$ [km/s], $19\leq|k|\leq 27$'); ax.legend(); ax.grid(lw=0.4)
ax.text(0.02,0.04,'dotted: binary switch (c16340)',transform=ax.transAxes,fontsize=9)
fig.savefig('fig_ring_growth.png',dpi=100)
d=np.load(sorted(glob.glob('out/q_*.npz'))[-1]); P=np.load('prof.npz'); fig,ax=plt.subplots(figsize=(7,4.6),constrained_layout=True)
for n,c in ((300,'C0'),(360,'C1'),(400,'C2'),(430,'r')):
    f='out/q_%05d.npz'%n; d=np.load(f); KY,KX=np.meshgrid(d['kx'],d['ky']); kk=np.hypot(KX,KY); ring=(kk>=18.5)&(kk<27.5); rr=d['r'][int(d['i0']):]/RP; ax.plot(rr,np.sqrt((abs(d['C'])**2)[ring].sum(0))/1e5,c,label='t = %.2f d'%(float(d['t'])/86400))
ax.set_xlim(0.706,0.80); ax.set_xlabel(r'$r/R_{\rm ph,IC}$'); ax.set_ylabel('ring rms $v_r$ [km/s]'); a2=ax.twinx(); a2.plot(P['x'],np.sqrt(np.maximum(P['N2'],0))*1e4,'k--',lw=1.2); a2.axhline(1.64,color='0.5',ls=':'); a2.set_ylim(0,4); a2.set_ylabel(r'$N$ [$10^{-4}$ rad/s] (dashed)'); ax.legend(loc='upper left')
a2.text(0.99,0.03,'dotted grey: measured 1.64e-4 rad/s',transform=a2.transAxes,ha='right',fontsize=9); fig.savefig('fig_ring_profile.png',dpi=100)
