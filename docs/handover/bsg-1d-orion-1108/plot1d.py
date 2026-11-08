import numpy as np, glob, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
R='/orion/ptmp/jinma/athpp_repro3d/run1d'; D=0.0669426
runs=[('lhx','1e-6',0,'--'),('lhx_vf','1e-6',0,'-'),('lhx_e8','1e-8',1,'--'),('lhx_vf_e8','1e-8',1,'-'),('lhx_e10','1e-10',2,'--'),('lhx_vf_e10','1e-10',2,'-')]
col=['#eb6834','#2a78d6','#1baf7a']
def rph(args):
    d,n=args; o=f'{R}/{d}/output'
    a=np.vstack([np.loadtxt(f'{o}/bsg.block{b}.out2.{n:05d}.tab',comments='#',ndmin=2) for b in range(32)])
    t=float(open(f'{o}/bsg.block0.out2.{n:05d}.tab').readline().split('time=')[1].split()[0])*D
    f=a[np.argsort(a[:,1])]; r=f[:,1]; sig=f[:,-3]+f[:,-2]
    rf=np.concatenate([[r[0]-(r[1]-r[0])/2],(r[1:]+r[:-1])/2,[r[-1]+(r[-1]-r[-2])/2]]); dr=np.diff(rf)
    tau=(sig*dr)[::-1].cumsum()[::-1]-.5*sig*dr; i=np.nonzero(tau>1)[0][-1]
    return t, r[i]+(r[i+1]-r[i])*np.log(tau[i])/np.log(tau[i]/tau[i+1])
if not os.path.exists('rph1d.npz'):
    out={}
    with ProcessPoolExecutor(16) as ex:
        for d,*_ in runs:
            out[d]=np.array(list(ex.map(rph,[(d,n) for n in list(range(0,402,2))+[401]])))
    np.savez('rph1d.npz',**out)
P=np.load('rph1d.npz')
fig,ax=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
ins=ax[1].inset_axes([0.08,0.08,0.45,0.38])
for d,e,c,ls in runs:
    lab=f"{'fixed v' if 'vf' in d else 'over-kicked v'}, {e}"
    t,r=P[d].T; ax[0].plot(t,r,ls=ls,c=col[c],lw=2,label=lab)
    h=np.loadtxt(f'{R}/{d}/output/bsg.hst',usecols=(0,9))
    keep=np.r_[True, h[1:,0]>np.maximum.accumulate(h[:,0])[:-1]]; h=h[keep]   # drop rerun duplicates
    ax[1].plot(h[:,0]*D,h[:,1]/h[0,1]-1,ls=ls,c=col[c],lw=2,label=lab)
    if c>0: ins.plot(h[:,0]*D,100*(h[:,1]/h[0,1]-1),ls=ls,c=col[c],lw=1.5)
ax[0].set(xlabel='time [d]',ylabel=r'$r_{\rm photo}$ ($\tau=1$) [$R_\odot$]',title='Photospheric radius')
ax[1].set(xlabel='time [d]',ylabel=r'$E_{\rm tot}(t)/E_{\rm tot}(0)-1$',title='Total energy')
ax[1].axhline(0,c='0.6',lw=.8)
ins.set_xlim(0,200); ins.axhline(0,c='0.6',lw=.6); ins.tick_params(labelsize=7); ins.grid(alpha=.25)
ins.set_title('1e-8 and 1e-10 only, in %',fontsize=8)
for a in ax: a.grid(alpha=.25); a.set_xlim(0,200); [a.spines[s].set_visible(False) for s in ('top','right')]
ax[0].legend(loc='lower left',frameon=False,fontsize=8.5,handlelength=2.5)
fig.suptitle('1-D Athena++ BSG (Linhao setup, cfl_rad 0.025): Doppler velocity × error limit',fontsize=11)
fig.savefig('rph1d_six.png',dpi=150)
for d,*_ in runs: print(d, P[d][0], P[d][-1])
