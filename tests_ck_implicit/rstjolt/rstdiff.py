import sys, numpy as np, os
# per-MB record layout for the dhj MHD (general EOS) restart: u0(nv), b1,b2,b3, wtemp, wder(2)
nv=int(sys.argv[3]) if len(sys.argv)>3 else 5
mhd = (sys.argv[4]=='mhd') if len(sys.argv)>4 else True
n1,n2,n3=132,20,20; N=n1*n2*n3; nmb=24
parts=[('u%d'%i,(n3,n2,n1)) for i in range(nv)]
if mhd: parts+= [('b1',(n3,n2,n1+1)),('b2',(n3,n2+1,n1)),('b3',(n3+1,n2,n1))]
parts+=[('wtemp',(n3,n2,n1)),('wdPR',(n3,n2,n1)),('wdG1',(n3,n2,n1))]
rec=sum(np.prod(s) for _,s in parts)
def load(f):
    sz=os.path.getsize(f); off=sz-nmb*rec*8
    a=np.fromfile(f,dtype='<f8',offset=off).reshape(nmb,rec); return a
A=load(sys.argv[1]); B=load(sys.argv[2])
p=0
for name,s in parts:
    n=int(np.prod(s)); a=A[:,p:p+n].reshape((nmb,)+s); b=B[:,p:p+n].reshape((nmb,)+s); p+=n
    d=np.abs(a-b); r=d/np.maximum(np.abs(a),1e-300)
    nd=np.count_nonzero(d)
    if nd:
        idx=np.unravel_index(np.argmax(r),r.shape)
        # radial distribution of differing cells (active)
        ii=np.nonzero(d)[-1]
        print('%-6s ndiff %8d maxrel %.3e at (m,k,j,i)=%s  i-range %d..%d'%(name,nd,r.max(),idx,ii.min(),ii.max()))
    else: print('%-6s identical'%name)
