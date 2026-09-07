import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'); import dhjcs
S='/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'; B='/viper/u2/jinma/ATHENAK/bench/'
eos=dhjcs.EOS(S+'/eostab/eos_table.txt')
runs={'cs hydro':('cs',B+'cs_prod_hyd_rot/bin/dhj.hydro_w.00010.bin'),
      'cs MHD (fixed BC)':('cs',B+'cs_inflow_bcfix/bin/dhj.mhd_w_bcc.00020.bin'),
      'sp hydro':('sp',B+'sp_dhj_hyd/bin/dhj.hydro_w.00010.bin')}
COL={'cs hydro':'#1f5fbf','cs MHD (fixed BC)':'#c8451e','sp hydro':'#5c5c5c'}
D={k:dhjcs.load(p,eos,kind) for k,(kind,p) in runs.items()}
levels=[1e-3,1e-2,1e-1,1.0]; lev=np.log10(np.array(levels)*1e6)
# ---- Fig 1: T maps
fig,ax=plt.subplots(len(levels),3,figsize=(15,3.2*len(levels)),constrained_layout=True)
for c,(k,d) in enumerate(D.items()):
    Tl=dhjcs.level_interp(d['T'],d['lp'],lev); lats,lons,Tg=dhjcs.latlon_bin(Tl,d['lat'],d['lon'],5,5)
    for l in range(len(levels)):
        vmin,vmax={1e-3:(800,3600),1e-2:(1500,3400),1e-1:(2400,3000),1.0:(3300,3500)}[levels[l]]
        im=ax[l,c].pcolormesh(lons,lats,Tg[l],cmap='magma',vmin=vmin,vmax=vmax,shading='nearest'); plt.colorbar(im,ax=ax[l,c],label='T [K]')
        ax[l,c].set_title(f'{k}, rot {d["rot"]:.0f}: T at {levels[l]:g} bar'); ax[l,c].set_xlabel('lon from substellar [deg]'); ax[l,c].set_ylabel('lat [deg]')
fig.savefig(S+'/fig1_Tmaps.png',dpi=110); plt.close(fig)
# ---- Fig 2: T(p) at four points + zonal-mean wind (lat,p) + equatorial u(p)
plev=np.logspace(-6,2,81); llev=np.log10(plev*1e6)
fig,ax=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
pts={'substellar':(0,0),'antistellar':(0,180),'east terminator':(0,90),'west terminator':(0,-90)}
ls={'substellar':'-','antistellar':'--','east terminator':':','west terminator':'-.'}
for k,d in D.items():
    Tl=dhjcs.level_interp(d['T'],d['lp'],llev); lats,lons,Tg=dhjcs.latlon_bin(Tl,d['lat'],d['lon'],10,10)
    for name,(la,lo) in pts.items():
        i=np.argmin(np.abs(lats-la)); j=np.argmin(np.abs(((lons-lo+180)%360)-180))
        ax[0,0].plot(Tg[:,i,j],plev,ls[name],color=COL[k],label=f'{k} {name}' if name in ('substellar','antistellar') else None)
ax[0,0].set_yscale('log'); ax[0,0].invert_yaxis(); ax[0,0].set_xlabel('T [K]'); ax[0,0].set_ylabel('p [bar]'); ax[0,0].legend(fontsize=7); ax[0,0].set_title('T(p): solid substellar, dashed antistellar,\ndotted east / dash-dot west terminator')
for c,(k,d) in enumerate(D.items()):
    ul=dhjcs.level_interp(d['u'],d['lp'],llev); lats,uz=dhjcs.zonal_mean(ul,d['lat'],5)
    a=ax[0,1] if c==0 else (ax[0,2] if c==1 else ax[1,0])
    im=a.pcolormesh(lats,plev,uz/1e5,cmap='RdBu_r',vmin=-3,vmax=3,shading='nearest'); a.set_yscale('log'); a.invert_yaxis(); plt.colorbar(im,ax=a,label='zonal-mean u [km/s]')
    a.set_title(f'{k}, rot {d["rot"]:.0f}: zonal-mean zonal wind'); a.set_xlabel('lat [deg]'); a.set_ylabel('p [bar]')
    eq=np.abs(lats)<12.5
    ax[1,1].plot(np.nanmean(uz[:,eq],axis=1)/1e5,plev,color=COL[k],label=k)
    ax[1,2].plot(np.nanmax(np.abs(ul.reshape(len(llev),-1)),axis=1)/1e5,plev,color=COL[k],label=k)
ax[1,1].set_yscale('log'); ax[1,1].invert_yaxis(); ax[1,1].set_xlabel('equatorial (|lat|<12.5) zonal-mean u [km/s]'); ax[1,1].set_ylabel('p [bar]'); ax[1,1].legend(); ax[1,1].axvline(0,color='k',lw=0.5)
ax[1,2].set_yscale('log'); ax[1,2].invert_yaxis(); ax[1,2].set_xlabel('max |u| on the isobar [km/s]'); ax[1,2].set_ylabel('p [bar]'); ax[1,2].legend()
fig.savefig(S+'/fig2_Tp_winds.png',dpi=110); plt.close(fig)
# ---- Fig 4: field maps for the MHD run
d=D['cs MHD (fixed BC)']; blev=[1e-2,1.0,10.0,100.0]; bl=np.log10(np.array(blev)*1e6)
fig,ax=plt.subplots(len(blev),3,figsize=(15,3.2*len(blev)),constrained_layout=True)
Bl=dhjcs.level_interp(d['Babs'],d['lp'],bl); Brl=dhjcs.level_interp(d['Br'],d['lp'],bl); beta=dhjcs.level_interp(d['p']/(d['Babs']**2/(8*np.pi)+1e-30),d['lp'],bl)
lats,lons,Bg=dhjcs.latlon_bin(Bl,d['lat'],d['lon'],5,5); _,_,Brg=dhjcs.latlon_bin(Brl,d['lat'],d['lon'],5,5); _,_,bg=dhjcs.latlon_bin(np.log10(beta),d['lat'],d['lon'],5,5)
for l in range(len(blev)):
    im=ax[l,0].pcolormesh(lons,lats,np.log10(Bg[l]),cmap='magma',shading='nearest'); plt.colorbar(im,ax=ax[l,0],label='log10 |B| [G]'); ax[l,0].set_title(f'|B| at {blev[l]:g} bar (rot {d["rot"]:.0f})')
    m=np.nanmax(np.abs(Brg[l])); im=ax[l,1].pcolormesh(lons,lats,Brg[l],cmap='RdBu_r',vmin=-m,vmax=m,shading='nearest'); plt.colorbar(im,ax=ax[l,1],label='B_r [G]'); ax[l,1].set_title(f'B_r at {blev[l]:g} bar')
    im=ax[l,2].pcolormesh(lons,lats,bg[l],cmap='magma',shading='nearest'); plt.colorbar(im,ax=ax[l,2],label='log10 plasma beta'); ax[l,2].set_title(f'beta at {blev[l]:g} bar')
    for a in ax[l]: a.set_xlabel('lon from substellar [deg]'); a.set_ylabel('lat [deg]')
fig.savefig(S+'/fig4_field.png',dpi=110); plt.close(fig)
print('done')
