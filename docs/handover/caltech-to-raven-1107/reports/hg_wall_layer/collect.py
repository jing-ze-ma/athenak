#!/usr/bin/env python3
"""READ ONLY on prodG. hst last-writer-wins dedupe, 0.1-d bins; shell profiles from page red_*.npz; stale checks."""
import numpy as np, re, glob, os, time
from concurrent.futures import ThreadPoolExecutor
RUN='/resnick/groups/carnegie_poc/jingze/hgglob_1021/prodG'; RED='/resnick/groups/carnegie_poc/jingze/hgglob_page/out/data'
OUT='/resnick/groups/carnegie_poc/jingze/hg_relax_1106'
def hst(fn):
    hdr=[l for l in open(fn) if l.startswith('#  [')]
    names=re.findall(r'\[\d+\]=(\S+)',hdr[-1])
    rows=np.loadtxt(fn,comments='#',invalid_raise=False) if False else None
    a=[]
    with open(fn) as f:
        for l in f:
            if l.startswith('#'): continue
            p=l.split()
            if len(p)!=len(names): continue
            try: a.append([float(x) for x in p])
            except ValueError: pass
    a=np.array(a)
    # last writer wins: walk backwards keeping a strictly decreasing time
    t=a[:,0]; cm=np.minimum.accumulate(t[::-1])[::-1]
    keep=np.zeros(len(t),bool); m=np.inf
    for i in range(len(t)-1,-1,-1):
        if t[i]<m: keep[i]=True; m=t[i]
    # rewinds
    rew=np.where(np.diff(t)<0)[0]
    print(os.path.basename(fn),'rows',len(t),'kept',keep.sum(),'rewinds',len(rew),[(round(t[i]/86400,2),round(t[i+1]/86400,2)) for i in rew],flush=True)
    return names,a[keep]
for tag,fn in (('hyd','hegiant.hydro.hst'),('usr','hegiant.user.hst')):
    names,a=hst(os.path.join(RUN,fn))
    tb=np.arange(0,a[-1,0]+8640,8640.)
    ib=np.digitize(a[:,0],tb)-1
    n=np.bincount(ib,minlength=len(tb))
    bm=np.array([np.bincount(ib,weights=a[:,k],minlength=len(tb))/np.maximum(n,1) for k in range(a.shape[1])]).T
    bm=bm[n>0]
    np.savez(f'{OUT}/hst_{tag}.npz',names=np.array(names),a=bm,n=n[n>0],last=a[-1])
    np.save(f'{OUT}/hst_{tag}_full.npy',a)
    print(tag,names,'last t d',a[-1,0]/86400,flush=True)
KEYS='t rho mvr vrrms vtrms E Trad Tgas eint Lr Lconv Lconv_gas Lradadv Lrcom Lbulk Lradbulk Lkin dlrho_rms Tfloor_frac Rph_med Rph_p2 Rph_p98 Teff_med Rph_nan nonfinite vet ver'.split()
fh=sorted(glob.glob(RUN+'/bin/hegiant.hydro_w.*.bin'))
def hdr(fn):
    with open(fn,'rb') as f: m=re.search(rb'time=\s*([0-9.eE+-]+)',f.read(4096))
    return int(fn.split('.')[-2]),float(m.group(1)),os.path.getmtime(fn)
with ThreadPoolExecutor(4) as ex: H=list(ex.map(hdr,fh))
H=np.array(H); print('dumps',len(H),'t range d',H[0,1]/86400,H[-1,1]/86400,'non-monotonic',(np.diff(H[:,1])<=0).sum(),'gaps>1.01 cadence',(np.diff(H[:,1])>29000).sum(),flush=True)
def rd(i):
    idx,t,mt=H[i]; o=f'{RED}/red_{int(idx):05d}.npz'
    if not os.path.exists(o): return None
    z=np.load(o); d={k:(z[k] if k in z else np.nan) for k in KEYS}; d['idx']=idx; d['thdr']=t; d['stale']=abs(float(d['t'])-t)>1.0; d['older']=os.path.getmtime(o)<mt
    if i==0: d['r']=z['r']; d['rf']=z['rf']
    return d
with ThreadPoolExecutor(4) as ex: R=[x for x in ex.map(rd,range(len(H))) if x is not None]
print('reduced read',len(R),'stale(time mismatch)',sum(x['stale'] for x in R),'older than dump',sum(x['older'] for x in R),flush=True)
good=[x for x in R if not x['stale']]
S={k:np.array([x[k] for x in good]) for k in KEYS+['idx']}
S['r']=R[0]['r']; S['rf']=R[0]['rf']
np.savez(f'{OUT}/prof.npz',**S)
print('saved',S['rho'].shape,'t',S['t'][0]/86400,S['t'][-1]/86400,'dump mtimes saved',flush=True)
np.save(f'{OUT}/dump_hdr.npy',H)
print('DONE',flush=True)
