import sys, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc
Rs=6.957e10; Tu=8.1063567271e-09  # K per code (p/d) unit
D='/viper/ptmp2/jinma/t19l64_1119/prod_viper/bin/t19l64.hydro_w.%05d.bin'
for n in [int(x) for x in sys.argv[1:]]:
    f=bc.read_binary(D%n); md=f['mb_data']; vn=f['var_names']
    nmb=f['n_mbs']
    if 'mb_x1v' in f: r=f['mb_x1v'][0]
    else:
        g=np.loadtxt('/viper/ptmp2/jinma/t19l64_1119/t19l64.x1grid.txt'); r=0.5*(g[1:]+g[:-1]) if g.ndim==1 else g[:,0]
    rows=[]
    vmaxs=[]; 
    acc={}
    for m in range(nmb):
        d=md['dens'][m]; v1=md['velx'][m]; v2=md['vely'][m]; v3=md['velz'][m]; p=md['eint'][m] if 'eint' in md else md['pres'][m]
        vv=np.sqrt(v1**2+v2**2+v3**2)
        for thr in (1e8,3e8):
            k,j,i=np.nonzero(vv>thr)
            acc.setdefault(thr,[]).append(np.c_[r[i], d[k,j,i], vv[k,j,i], v1[k,j,i], p[k,j,i]/d[k,j,i]])
        vmaxs.append(vv.max())
    print(f"== dump {n} t {f['time']/3600:.1f} h cycle {f['cycle']}  vars {vn}  max|v| {max(vmaxs):.2e}")
    for thr,L in acc.items():
        a=np.concatenate(L)
        if len(a)==0: print(f"  |v|>{thr:.0e}: 0 cells"); continue
        rr=a[:,0]/Rs
        print(f"  |v|>{thr:.0e}: {len(a)} cells; r {np.percentile(rr,1):.1f}/{np.median(rr):.1f}/{np.percentile(rr,99):.1f} Rsun (1/50/99%); rho median {np.median(a[:,1]):.1e} (min {a[:,1].min():.1e} max {a[:,1].max():.1e}); frac rho<1e-17 {np.mean(a[:,1]<1e-17):.2f}; v_r>0 frac {np.mean(a[:,3]>0):.2f}; T(p/rho) median {np.median(a[:,4])*Tu*(2/3):.1e} K")
        h,e=np.histogram(rr,bins=[14,20,25,30,32,35,40,50,75,100,150,200,250,300,340]); print('   r-hist', ' '.join(f"{int(e[q])}-{int(e[q+1])}:{h[q]}" for q in range(len(h)) if h[q]))
    sys.stdout.flush()
