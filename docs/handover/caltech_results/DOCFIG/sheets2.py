import sys, os, json, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/wt_note_a0dc/docs/handover/plot_style')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, idort_style
import cmasher as cmr, cmocean; from cmcrameri import cm as cr
from matplotlib.ticker import MaxNLocator
from copy import copy
idort_style.use(font=13); plt.rcParams.update({'axes.linewidth':1.1,'figure.dpi':100,'savefig.dpi':100,'savefig.bbox':None,'savefig.format':'png','axes.titlesize':13,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'xtick.major.size':4,'ytick.major.size':4,'xtick.minor.size':2,'ytick.minor.size':2})
z=np.load('r2_vol.npz'); r=z['r'][:200]; rf=z['rf'][:201]; rho=z['rho'].astype(float); T=z['T'].astype(float); v=z['v'].astype(float)/1e5; rm=z['rm'][:200]; Tm=z['Tm'][:200]; fl=float(z['fl'])
OUT=sys.argv[1]; os.makedirs(OUT,exist_ok=True); GREY='0.78'; CZ=(38.5,44.0); k0=128
th=np.linspace(np.pi/3,2*np.pi/3,257); ph=np.linspace(0,60,257); thd=np.degrees(th)
i0=int(np.searchsorted(rf,34.0)); i1=int(np.searchsorted(rf,54.0)); R,TH=np.meshgrid(rf[i0:i1+1],th); X=R*np.sin(TH); Z=R*np.cos(TH)
cz=(r>=CZ[0])&(r<=CZ[1]); lim={}
mask=rho<=1.5*fl; dr=np.ma.masked_where(mask,rho/rm-1); dT=np.ma.masked_where(mask,T/Tm-1)
def cm_(c): 
    c=copy(c); c.set_bad(GREY); return c
j0=int(np.searchsorted(rf,37.0)); j1=int(np.searchsorted(rf,46.0)); Rb,THb=np.meshgrid(rf[j0:j1+1],th); Xb=Rb*np.sin(THb); Zb=Rb*np.cos(THb)
def slice_(a,q,cmap,vmin,vmax,title,label,crop=False):
    XX,ZZ,qq=(Xb,Zb,q[k0][:,j0:j1]) if crop else (X,Z,q[k0][:,i0:i1])
    pc=a.pcolormesh(ZZ,XX,qq,cmap=cm_(cmap),vmin=vmin,vmax=vmax,rasterized=True,shading='flat'); a.set_aspect('equal'); a.set_xlim(-27.5,27.5); a.set_ylim(*((31.0,47.0) if crop else (28.5,54.5))); a.set_xlabel(r'$\rm z\ [R_\odot]$'); a.set_ylabel(r'$\rm x\ [R_\odot]$'); a.set_title(title); cb=plt.colorbar(pc,ax=a,shrink=0.6,extend='both'); cb.set_label(label); cb.locator=MaxNLocator(5); cb.update_ticks()
def lateral(a,f,cmap,vmin,vmax,title,label):
    pc=a.pcolormesh(ph,thd,f.T,cmap=cm_(cmap),vmin=vmin,vmax=vmax,rasterized=True,shading='flat'); a.set_aspect('equal'); a.set_ylim(120,60); a.set_xlabel(r'$\rm \phi\ [deg]$'); a.set_ylabel(r'$\rm \theta\ [deg]$'); a.set_title(title); cb=plt.colorbar(pc,ax=a,shrink=0.62,extend='both'); cb.set_label(label); cb.locator=MaxNLocator(5); cb.update_ticks()
DIV=[('cmasher fusion',cmr.fusion),('cmocean balance',cmocean.cm.balance),('crameri vik',cr.vik)]
# slice limits: 99 percentile of |fluctuation| of the convection zone cells of this slice
sel=dr[k0][:,cz]; vd=float(np.percentile(abs(sel.compressed()),99)); selT=dT[k0][:,cz]; vT=float(np.percentile(abs(selT.compressed()),99)); lim['slice drho']=vd; lim['slice dT']=vT
lo=np.log10(rho[k0][:,cz]); lT=np.log10(T[k0][:,cz]); bd=np.percentile(lo,[1,99]); bT=np.percentile(lT,[1,99]); lim['log rho']=list(bd); lim['log T']=list(bT)
for nm,q,raw,vv,b,lab,labr,fn in (('density',dr,np.ma.masked_where(mask,np.log10(rho)),vd,bd,r'$\rm \delta\rho/\langle\rho\rangle$',r'$\rm \log_{10}\rho\ [g\ cm^{-3}]$','r2_slice_density.png'),('temperature',dT,np.ma.masked_where(mask,np.log10(T)),vT,bT,r'$\rm \delta T/\langle T\rangle$',r'$\rm \log_{10}T\ [K]$','r2_slice_temperature.png')):
    fig,ax=plt.subplots(3,2,figsize=(12,9.4),constrained_layout=True)
    for j,(n,c) in enumerate(DIV): slice_(ax[j,0],q,c,-vv,vv,'(a) fluctuation, '+n,lab)
    for j,(n,c) in enumerate((('cmocean thermal',cmocean.cm.thermal),('cmasher ember',cmr.ember))): slice_(ax[j,1],raw,c,b[0],b[1],'(b) zoomed log, '+n,labr,crop=True)
    ax[2,1].axis('off'); ax[2,1].text(0.02,0.75,'grey: floor gas (density <= 1.5 x floor)\n(a): 34-54 Rsun, symmetric limits =\n99th percentile of |fluctuation| in 38.5-44.0 Rsun\n(b): 37-46 Rsun, limits = 1-99 percentile\nof the 38.5-44.0 Rsun cells',fontsize=13,va='top',transform=ax[2,1].transAxes)
    fig.savefig(os.path.join(OUT,fn)); plt.close(fig)
S=[int(np.argmin(abs(r-x))) for x in (40.70,43.52)]; lim['shells']=[float(r[s]) for s in S]
for nm,q,lab,fn,vv in (('density',dr,r'$\rm \delta\rho/\langle\rho\rangle$','r2_lateral_density.png',None),('temperature',dT,r'$\rm \delta T/\langle T\rangle$','r2_lateral_temperature.png',None)):
    fig,ax=plt.subplots(2,3,figsize=(12.6,7.6),constrained_layout=True)
    for a,s in enumerate(S):
        f=q[:,:,s]; V=float(np.percentile(abs(f.compressed()),99)); lim['lateral %s r %.2f'%(nm,r[s])]=V
        for j,(n,c) in enumerate(DIV): lateral(ax[a,j],f,c,-V,V,'%s, %.1f'%(n,r[s])+r'$\rm\ R_\odot$',lab)
    fig.savefig(os.path.join(OUT,fn)); plt.close(fig)
# v_r reference
s=S[0]; vm=float(np.percentile(abs(v[:,:,s]),99.5)); lim['v_r shell %.2f'%r[s]]=vm
fig,a=plt.subplots(figsize=(9,7.4),constrained_layout=True); lateral(a,v[:,:,s],cmr.fusion,-vm,vm,'cmasher fusion (settled)',r'$\rm v_r\ [km\ s^{-1}]$'); fig.savefig(os.path.join(OUT,'r2_vr_fusion.png')); plt.close(fig)
# flux
zq=np.load('bsg114.npz'); q=zq['fr']/zq['fr'].mean()
F=[('Blues (low flux light, high flux dark)',plt.get_cmap('Blues'),np.percentile(q,[2,98])),('Blues_r (low flux dark, high flux light)',plt.get_cmap('Blues_r'),np.percentile(q,[2,98])),('cmocean dense (low light, high dark)',cmocean.cm.dense,np.percentile(q,[2,98])),('cmocean ice (low dark, high light)',cmocean.cm.ice,np.percentile(q,[2,98]))]
fig,ax=plt.subplots(2,2,figsize=(12,10.2),constrained_layout=True)
for a,(n,c,l) in zip(ax.flat,F): lateral(a,q,c,l[0],l[1],n,r'$\rm F_r/\langle F_r\rangle$')
fig.savefig(os.path.join(OUT,'r2_flux.png')); plt.close(fig); lim['flux 2-98 percentile']=[float(x) for x in F[0][2]]; lim['flux 1-99 (round 1)']=[float(x) for x in np.percentile(q,[1,99])]
json.dump(lim,open(os.path.join(OUT,'limits.json'),'w'),indent=1); print(json.dumps(lim))
