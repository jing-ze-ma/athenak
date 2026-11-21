import sys, numpy as np
sys.path.insert(0, sys.argv[3] if len(sys.argv)>3 else 'vis/python')
import bin_convert as bc
# fraction of cells per shell with off-diagonals exactly 0 (the J_ray < clean_ray_min E isotropic path)
d=bc.read_binary(sys.argv[1]); mb=d['mb_data']
z=None; n=0
for m in range(d['n_mbs']):
    a=(mb['m1_d12'][m]==0)&(mb['m1_d13'][m]==0)&(mb['m1_d23'][m]==0)
    s=a.sum(axis=(0,1)); z=s if z is None else z+s; n+=a.shape[0]*a.shape[1]
f=z/n; st=int(sys.argv[2]) if len(sys.argv)>2 else 16
print('time',d['time'],'total invalid frac %.4f'%f.mean())
print(' '.join('%d:%.2f'%(i,f[i]) for i in range(0,len(f),st)))
