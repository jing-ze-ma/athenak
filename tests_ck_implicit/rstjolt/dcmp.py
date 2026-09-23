import numpy as np, sys
dA,dB=(sys.argv[1],sys.argv[2]) if len(sys.argv)>2 else ('dA','dB')
order=['u0','w0','bcc0','b0x1','b0x2','b0x3','wder','wtemp','u1','b1x1','b1x2','b1x3','efx1','efx2','efx3','wsaved','bccsaved','phicc0','wbq0','pwb','u0wb','w0wb','eta_b']
sz={}
for l in open(dA+'/log'):
    if l.startswith('RSTJOLT'): _,n,s=l.split(); sz[n]=int(s)
n1,n2,n3,nmb=132,20,20,4; N=nmb*n1*n2*n3
for r in range(6):
    a=np.fromfile(dA+'/dump_%d.bin'%r); b=np.fromfile(dB+'/dump_%d.bin'%r); p=0
    for nm in order:
        n=sz[nm]; x=a[p:p+n]; y=b[p:p+n]; p+=n
        d=(x!=y)&~(np.isnan(x)&np.isnan(y))
        if d.any():
            rel=np.abs(x-y)/np.maximum(np.abs(x),1e-300)
            info=''
            if n%N==0:
                nv=n//N; idx=np.unravel_index(np.nonzero(d)[0],(nmb,nv,n3,n2,n1)); info='vars %s i %d..%d'%(np.unique(idx[1]),idx[4].min(),idx[4].max())
            print('rank',r,nm,'ndiff',d.sum(),'of',n,'maxrel %.3e'%np.nanmax(rel[d]),info)
