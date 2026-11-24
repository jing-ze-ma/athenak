import sys, types, numpy as np
sys.modules.setdefault('h5py', types.ModuleType('h5py'))
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/vis/python')
import bin_convert as bc
D='/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/m35zams'; k=sys.argv[1]
v=bc.read_binary(D+'/bin/m35zams.vet.%s.bin'%k); h=bc.read_binary(D+'/bin/m35zams.hydro_w.%s.bin'%k)
f1=np.array(v['mb_data']['vet_f1'],dtype=float); f3=np.array(v['mb_data']['vet_f3'],dtype=float); f2=np.array(v['mb_data']['vet_f2'],dtype=float); rho=np.array(h['mb_data']['dens'],dtype=float)
x1=np.loadtxt(D+'/m35zams.x1grid.txt',comments='#',max_rows=390)[3:387,2]
print('time %.6g cycle %d (vet bin %s; float32 output)'%(v['time'],v['cycle'],k))
rm=rho.mean(axis=(0,1,2))
for lo,hi,lab in ((1e-16,4e-13,'thin layer rho 1e-16..4e-13 (cells)'),(1e-14,4e-13,'rho 1e-14..4e-13 (cells)')):
    m=(rho>=lo)&(rho<=hi)&(np.abs(f1)>0)
    r=np.abs(f3[m])/np.abs(f1[m]); r2=np.abs(f2[m])/np.abs(f1[m])
    print('%s: n=%d  max|F_phi|/|F_r| %.3e  p99 %.3e  median %.3e | max|F_theta|/|F_r| %.3e'%(lab,m.sum(),r.max(),np.percentile(r,99),np.median(r),r2.max()))
S3=np.abs(f3).sum(axis=(0,1,2)); S1=np.abs(f1).sum(axis=(0,1,2))
band=np.where((rm>=1e-16)&(rm<=4e-13))[0]
print('shell ratio sum|F_phi|/sum|F_r| over the thin-layer shells i=%d..%d (r %.4e..%.4e): max %.3e'%(band.min(),band.max(),x1[band.min()],x1[band.max()],(S3[band]/S1[band]).max()))
