import sys, numpy as np
sys.path.insert(0,'/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r14')
from dtdiag import rd
RS=2.3717e11; GM=6.674e-8*6.2651e33
for arm in sys.argv[1:]:
    R=rd(arm+'/rt_profile.bin'); r=R[0][1]
    t=np.array([x[0] for x in R]); sel=t<=30.0
    mom=np.array([x[2][2] for x in R])[sel]   # rho v1
    tt=t[sel]
    A=np.vstack([tt,np.ones_like(tt)]).T
    f=np.linalg.lstsq(A,mom,rcond=None)[0][0]     # d(rho v1)/dt per cell
    g=R[0][2][0]*GM/r**2
    ks=[k for k in range(len(r)) if 0.90<r[k]/RS<0.95]
    print("== %s  n=%d t<=%.0f"%(arm,sel.sum(),tt[-1]))
    print("   f/(rho g):"," ".join("%.3f:%+.4f"%(r[k]/RS,f[k]/g[k]) for k in ks))
    m=max(abs(f[k]/g[k]) for k in ks)
    print("   max|f/rho g| over 0.90-0.95R = %.4f"%m)
