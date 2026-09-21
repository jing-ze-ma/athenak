import numpy as np, struct, sys
def rd(p):
    f=open(p,'rb'); recs=[]
    while True:
        b=f.read(8)
        if len(b)<8: break
        t=struct.unpack('d',b)[0]; n1,nv=struct.unpack('ii',f.read(8))
        x=np.frombuffer(f.read(8*n1),dtype='f8').copy(); d=np.frombuffer(f.read(8*nv*n1),dtype='f8').reshape(nv,n1).copy()
        recs.append((t,x,d))
    return recs
def at(p,tt):
    r=rd(p); return min(r,key=lambda q:abs(q[0]-tt))
tt=float(sys.argv[1]); ref=at(sys.argv[2],tt)
print('reference %s at t=%.1f'%(sys.argv[2],ref[0]))
for p in sys.argv[3:]:
    q=at(p,tt)
    T0=ref[2][5]; T1=q[2][5]; F0=ref[2][9]; F1=q[2][9]; E0=ref[2][8]; E1=q[2][8]
    print('%-28s t=%.1f  max|dT/T|=%.3e  max|dF/F|=%.3e  max|dE/E|=%.3e'%(p,q[0],
        np.abs(T1/T0-1).max(), np.abs(F1/F0-1).max(), np.abs(E1/E0-1).max()))
