# usage: python3 stats.py ref arm...   (cpu/<run>/m1slab.{user,hydro}.hst, 200 s slab)
import sys, numpy as np
W='/viper/ptmp2/jinma/krylov_0923/cpu/'
FIN=2.475202e15
def load(r):
    u=np.loadtxt(W+r+'/m1slab.user.hst'); h=np.loadtxt(W+r+'/m1slab.hydro.hst')
    t=u[:,0]; late=t>=100.0
    return dict(F1top_Fin_mean=np.mean(u[late,2])/FIN, F1top_Fin_end=u[-1,2]/FIN,
                KE1_end=h[-1,7], KE2_end=h[-1,8], KE1_mean=np.mean(h[late,7]),
                dt_mean=np.mean(u[:,1]), totE_end=h[-1,6], V1max_end=u[-1,5])
ref=load(sys.argv[1]); keys=list(ref)
print('%-10s'%'run'+''.join('%14s'%k[:13] for k in keys))
print('%-10s'%sys.argv[1]+''.join('%14.6e'%ref[k] for k in keys))
for r in sys.argv[2:]:
    d=load(r); print('%-10s'%r+''.join('%14.2e'%(abs(d[k]-ref[k])/abs(ref[k])) for k in keys))
