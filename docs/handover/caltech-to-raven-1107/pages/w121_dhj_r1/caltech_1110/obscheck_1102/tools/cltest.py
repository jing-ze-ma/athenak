import sys, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/obscheck_1102/tools')
import em3d
m=em3d.load_model(sys.argv[1], '/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt')
X=m['X']; day = X[:,0] < 0
for ext in (1.0, 2.0, np.inf):
    M=em3d.cloud_mass(m, em3d.SETS['all'], ext)
    c=(M.sum(-1)>0)
    print('ext',ext,'cloudy cells frac %.3f; columns with cloud: day %.2f night %.2f; median cloud log p range night %s'%(c.mean(), c[day].any(1).mean(), c[~day].any(1).mean(),
       np.percentile(np.where(c[~day], m['lp'][~day], np.nan)[np.isfinite(np.where(c[~day], m['lp'][~day], np.nan))],[5,50,95]).round(1)))
lam=np.array([1.4,3.6,4.5])
for sp in em3d.COND:
    for R in (0.1,1,5): print(sp,R,'ext',em3d.cloud_kappa(sp,R*1e-4,'ext',lam).round(0),'abs',em3d.cloud_kappa(sp,R*1e-4,'abs',lam).round(1))
