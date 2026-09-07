import sys, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'); import dhjcs
S='/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'; B='/viper/u2/jinma/ATHENAK/bench/'
eos=dhjcs.EOS(S+'/eostab/eos_table.txt'); lev=np.log10(np.array([1e-3,1e-2,1e-1,1.0])*1e6)
series={'cs hydro':('cs',sorted(glob.glob(B+'cs_prod_hyd_rot/bin/*.bin'))[:31]),
        'cs MHD (drifting BC, rot 0-26)':('cs',sorted(glob.glob(B+'cs_prod_mhd_rot/bin/*.bin'))[:17]),
        'cs MHD (fixed BC, rot 8-20)':('cs',sorted(glob.glob(B+'cs_inflow_bcfix/bin/*.bin'))),
        'sp hydro':('sp',sorted(glob.glob(B+'sp_dhj_hyd/bin/*.bin'))[:31])}
COL={'cs hydro':'#1f5fbf','cs MHD (drifting BC, rot 0-26)':'#e8a33c','cs MHD (fixed BC, rot 8-20)':'#c8451e','sp hydro':'#5c5c5c'}
out={}
for k,(kind,fs) in series.items():
    rows=[]
    for f in fs:
        d=dhjcs.load(f,eos,kind); ul=dhjcs.level_interp(d['u'],d['lp'],lev); lats,uz=dhjcs.zonal_mean(ul,d['lat'],5); eq=np.abs(lats)<12.5
        rows.append([d['rot']]+list(np.nanmean(uz[:,eq],axis=1)/1e5)+list(np.nanmax(np.abs(ul.reshape(len(lev),-1)),axis=1)/1e5))
    out[k]=np.array(rows); print(k, out[k].shape)
np.save(S+'/jet_series.npy',out,allow_pickle=True)
fig,ax=plt.subplots(2,4,figsize=(17,7),constrained_layout=True)
for l in range(4):
    for k,a in out.items():
        ax[0,l].plot(a[:,0],a[:,1+l],'-o',ms=3,color=COL[k],label=k); ax[1,l].plot(a[:,0],a[:,5+l],'-o',ms=3,color=COL[k],label=k)
    ax[0,l].set_title(f'equatorial zonal-mean u at {10**lev[l]/1e6:g} bar'); ax[1,l].set_title(f'max |u| on the {10**lev[l]/1e6:g} bar isobar')
    for r in range(2): ax[r,l].set_xlabel('rotations'); ax[r,l].set_ylabel('km/s'); ax[r,l].axhline(0,color='k',lw=0.5)
ax[0,0].legend(fontsize=8)
fig.savefig(S+'/fig3_jet_time.png',dpi=110); print('done')
