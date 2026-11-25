import numpy as np, sys
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/wt_note_a0dc/docs/handover/plot_style')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, idort_style; idort_style.use(font=11)
import cmocean
z=np.load('f_00114.npz'); ph=(np.arange(256)+.5)*60/256; lat=30-(np.arange(256)+.5)*60/256
# region phi 22-40, lat 0-14
sp=np.where((ph>20)&(ph<42))[0]; st=np.where((lat>-2)&(lat<16))[0]
F=lambda k: z[k][np.ix_(sp,st)]
from numpy import percentile as P
fig,ax=plt.subplots(3,3,figsize=(13,12),constrained_layout=True)
for a,(k,l) in zip(ax.flat,(('Fc0.67','comoving F_r at tau 2/3 (interpolated)'),('Fc1.00','comoving F_r at tau 1 (interp.)'),('Fown','comoving F_r, own tau=1 cell (page)'),('r0.67','r(tau=2/3) [Rsun]'),('r1.00','r(tau=1) [Rsun]'),('iown','index of own tau=1 cell'),('Fshell','lab F_r at r=60 Rsun'),('Fs55','lab F_r at r=55 Rsun'),('rho0.67','log rho at tau=2/3'))):
    q=F(k); pc=a.pcolormesh(ph[sp],lat[st],q.T,cmap='cmo.ice',vmin=P(q,1),vmax=P(q,99)); a.set_title(l,fontsize=9); fig.colorbar(pc,ax=a,shrink=.8)
fig.savefig('rip.png',dpi=80)
# quantify: high-pass rms (minus 5x5 mean) in the ripple region and spectra across the ripple; correlation between fields
from scipy.ndimage import uniform_filter
def hp(q): return q-uniform_filter(q,9,mode='nearest')
R=[k for k in ('Fc0.67','Fc1.00','Fown','Fl0.67','E0.67','r0.67','r1.00','iown','Fshell','Fs55','Fs48','rho0.67','v0.67')]
B=(slice(None))
reg=np.ix_(np.where((ph>26)&(ph<36))[0],np.where((lat>1)&(lat<11))[0])
H={k:hp(z[k].astype(float))[reg] for k in R}
print('relative high-pass rms (9x9 cells) in phi 26-36, lat 1-11:',{k:float('%.2e'%(H[k].std()/abs(z[k][reg].mean()))) for k in R})
for k in R: print('corr of high-pass with Fc0.67 high-pass: %-8s %.2f'%(k,np.corrcoef(H['Fc0.67'].ravel(),H[k].ravel())[0,1]))
for n in (113,110,104):
    y=np.load('f_%05d.npz'%n); a=hp(y['Fc0.67'].astype(float))[reg]; b=H['Fc0.67']; print('dump',n,'t %.2f'%float(y['t']),'high-pass rel rms same region %.2e'%(a.std()/y['Fc0.67'][reg].mean()),'corr with dump 114 pattern %.2f'%np.corrcoef(a.ravel(),b.ravel())[0,1])
