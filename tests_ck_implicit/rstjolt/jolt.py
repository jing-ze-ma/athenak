import sys, numpy as np, os
# T = wtemp slab of a dhj restart (single file, 24 MB records); args: A B nv mhd|hyd
nv=5; mhd=sys.argv[3]=='mhd'
n1,n2,n3=132,20,20; N=n1*n2*n3; nmb=24
rec=nv*N+(n3*n2*(n1+1)+n3*(n2+1)*n1+(n3+1)*n2*n1 if mhd else 0)+3*N
toff=rec-3*N
def T(f):
    sz=os.path.getsize(f); a=np.fromfile(f,dtype='<f8',offset=sz-nmb*rec*8).reshape(nmb,rec)
    return a[:,toff:toff+N].reshape(nmb,n3,n2,n1)[:,2:-2,2:-2,2:-2]
a=T(sys.argv[1]); b=T(sys.argv[2]); r=np.abs(a-b)/a
i=np.unravel_index(np.argmax(r),r.shape)
print('max|dT|/T %.3e at (m,k,j,i_act)=%s T=%.1f dT=%.3e K; ncell>1e-6: %d, >1e-9: %d of %d; mean %.2e'%(r.max(),i,a[i],abs(a[i]-b[i]),(r>1e-6).sum(),(r>1e-9).sum(),r.size,r.mean()))
prof=r.max(axis=(0,1,2)); print('max over columns by i_act (every 8):',' '.join('%d:%.1e'%(k,prof[k]) for k in range(0,128,8)), '127:%.1e'%prof[127])
