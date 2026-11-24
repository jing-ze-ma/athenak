# thin-layer signed sums for M35ZAMS rc3 (orion NOTE 24o): shells i=320..331 (global x1 index), volume-weighted
import sys, glob, types, numpy as np
sys.modules.setdefault('h5py', types.ModuleType('h5py'))
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/vis/python')
import bin_convert as bc
D='/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/m35zams/bin'
I0,I1=320,331
print('%5s %9s %6s | %13s %13s | %13s %13s | %13s %13s' % ('dump','t','cycle','sumFphiV/|Fr|V','sumFthV/|Fr|V','<v_phi>_m','max|v_phi|','<v_r>_m','<|v_r|>_m'))
for f in sorted(glob.glob(D+'/m35zams.vet.*.bin')):
    k=f.split('.')[-2]; h=f.replace('.vet.','.hydro_w.')
    try: v=bc.read_binary(f); w=bc.read_binary(h)
    except Exception as e: print(k,'read error',e); continue
    x1f=np.array(v['mb_x1f']); x2f=np.array(v['mb_x2f']); x3f=np.array(v['mb_x3f'])
    sF3=sF2=sF1=0.0; sm=smv=smr=smar=0.0; vmax=0.0
    for b in range(len(v['mb_data']['vet_f1'])):
        r3=(x1f[b,1:]**3-x1f[b,:-1]**3)/3.0; dc=np.cos(x2f[b,:-1])-np.cos(x2f[b,1:]); dp=x3f[b,1:]-x3f[b,:-1]
        vol=dp[:,None,None]*dc[None,:,None]*r3[None,None,:]
        sl=slice(I0,I1+1)
        f1=np.array(v['mb_data']['vet_f1'][b],dtype=float)[:,:,sl]; f2=np.array(v['mb_data']['vet_f2'][b],dtype=float)[:,:,sl]; f3=np.array(v['mb_data']['vet_f3'][b],dtype=float)[:,:,sl]
        rho=np.array(w['mb_data']['dens'][b],dtype=float)[:,:,sl]; vr=np.array(w['mb_data']['velx'][b],dtype=float)[:,:,sl]; vp=np.array(w['mb_data']['velz'][b],dtype=float)[:,:,sl]
        V=vol[:,:,sl]; sF1+=(np.abs(f1)*V).sum(); sF2+=(f2*V).sum(); sF3+=(f3*V).sum()
        m=rho*V; sm+=m.sum(); smv+=(m*vp).sum(); smr+=(m*vr).sum(); smar+=(m*np.abs(vr)).sum(); vmax=max(vmax,np.abs(vp).max())
    print('%5s %9.4g %6d | %13.4e %13.4e | %13.4e %13.4e | %13.4e %13.4e' % (k, v['time'], v['cycle'], sF3/sF1, sF2/sF1, smv/sm, vmax, smr/sm, smar/sm))
print('(float32 bins; vet_f* = cell radiation flux moments; shells i=%d..%d = r %.4e..%.4e; v in cm/s)' % (I0, I1, x1f[0,I0], x1f[0,I1+1]))
