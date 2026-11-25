import sys, os, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/wt_note_a0dc/docs/handover/plot_style')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, idort_style
import cmasher as cmr, cmocean; from cmcrameri import cm as cr
from matplotlib.ticker import MaxNLocator
idort_style.use(font=14); plt.rcParams.update({'axes.linewidth':1.2,'figure.dpi':100,'savefig.dpi':100,'savefig.bbox':None,'savefig.format':'png','axes.titlesize':15,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12})
RS=6.957e10; z=np.load('bsg114.npz'); rf=z['rf']/RS; th=z['th']; ph=np.degrees(z['ph']); thd=np.degrees(z['th']); rmax=62.0; ni=int(np.searchsorted(rf,rmax)); TT,RR=np.meshgrid(th,rf[:ni+1],indexing='ij'); X=RR*np.sin(TT); Z=RR*np.cos(TT)
OUT=sys.argv[1]; os.makedirs(OUT,exist_ok=True)
def sheet(fn,kind,field,label,maps,vmin,vmax,sym=False,cbfmt=None):
    fig,ax=plt.subplots(2,2,figsize=(12,9.2),constrained_layout=True)
    for a,(nm,cm) in zip(ax.flat,maps):
        if kind=='rt':
            pc=a.pcolormesh(Z,X,field[:,:ni],cmap=cm,vmin=vmin,vmax=vmax,rasterized=True,shading='flat'); a.set_aspect('equal'); a.set_xlabel(r'$\rm z\ [R_\odot]$'); a.set_ylabel(r'$\rm x\ [R_\odot]$'); a.set_xlim(-rmax*0.53,rmax*0.53); a.set_ylim(15,rmax)
        else:
            pc=a.pcolormesh(ph,thd,field.T,cmap=cm,vmin=vmin,vmax=vmax,rasterized=True,shading='flat'); a.set_aspect('equal'); a.set_ylim(120,60); a.set_xlabel(r'$\rm \phi\ [deg]$'); a.set_ylabel(r'$\rm \theta\ [deg]$')
        a.set_title(nm); cb=fig.colorbar(pc,ax=a,shrink=0.92,extend='both'); cb.set_label(label); cb.locator=MaxNLocator(6); cb.update_ticks()
    fig.savefig(os.path.join(OUT,fn)); plt.close(fig)
T=np.log10(z['T']); d=np.log10(z['rho']); v=z['vr']/1e5; q=z['fr']/z['fr'].mean()
G=[('cmasher ember',cmr.ember),('cmasher sunburst',cmr.sunburst),('cmasher cosmic',cmr.cosmic),('cmocean thermal',cmocean.cm.thermal)]
S=[('cmocean balance',cmocean.cm.balance),('cmasher prinsenvlag',cmr.prinsenvlag),('cmasher fusion',cmr.fusion),('crameri vik',cr.vik)]
F=[('crameri lajolla',cr.lajolla),('crameri batlowW',cr.batlowW),('cmocean solar (added)',cmocean.cm.solar),('cmasher amber (added)',cmr.amber)]
pT=np.percentile(T[:,:ni],[0.5,99.5]); pd=np.percentile(d[:,:ni],[0.5,99.5]); vm=float(np.percentile(abs(v),99.5)); qm=np.percentile(q,[1,99])
sheet('palette_temperature.png','rt',T,r'$\rm \log_{10} T\ [K]$',G,*pT)
sheet('palette_density.png','rt',d,r'$\rm \log_{10} \rho\ [g\ cm^{-3}]$',G,*pd)
sheet('palette_vr.png','sh',v,r'$\rm v_r\ [km\ s^{-1}]$',S,-vm,vm)
sheet('palette_flux.png','sh',q,r'$\rm F_r/\langle F_r\rangle$ at $\rm \tau_R=1$',F,*qm)
print('limits: T',pT,'rho',pd,'v +-',vm,'F',qm,'shell r',float(z['rs'])/RS)
