import re, sys; sys.argv=['x']
src=open('/resnick/groups/carnegie_poc/jingze/bsglz_1124/bsglz.py').read().split("for lab, d, logs in")[0]; exec(src)
for lab,d in (('c1','prod_c1_1124/bsg'),('T2','prod_rc4v_1124/bsg')):
    n,a,_=rd(J+d+'/bsg3d.user.hst'); t=a[:,0]; lz=a[:,n.index('Lz_tot')]; lf=a[:,n.index('Lz_hflx')]; dt=np.diff(t)
    r=np.diff(lz)+0.5*(lf[1:]+lf[:-1])*dt; c=np.arange(1,len(t))                     # rows are one per cycle from c0 in these fresh runs (row index = cycle)
    base=np.abs(r[(c>1000)&(c<5000)]); print('%s: per-cycle residual c1000-5000 median %.2e max %.2e'%(lab,np.median(base),base.max()))
    big=np.where(abs(r)>1e3*np.median(base))[0]; print('   first cycle with |residual| > 1000 x that median: %s (n such cycles %d)'%(c[big[0]] if len(big) else None,len(big)))
    for lo,hi in ((5000,6000),(6000,6500),(6500,7000),(7000,7500),(7500,8000),(8000,8500),(8500,9000)):
        q=(c>=lo)&(c<hi)
        if q.sum(): print('   c%d-%d: sum %.3e, max |step| %.3e at c%d, steps > 1000 x median %d'%(lo,hi,r[q].sum(),abs(r[q]).max(),c[q][abs(r[q]).argmax()],(abs(r[q])>1e3*np.median(base)).sum()))
