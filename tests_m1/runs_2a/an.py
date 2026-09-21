import numpy as np, struct, sys
def rd(p):
    f=open(p,'rb'); recs=[]
    while True:
        b=f.read(8)
        if len(b)<8: break
        t=struct.unpack('d',b)[0]; n1,nv=struct.unpack('ii',f.read(8))
        x=np.frombuffer(f.read(8*n1),dtype='f8').copy()
        d=np.frombuffer(f.read(8*nv*n1),dtype='f8').reshape(nv,n1).copy()
        recs.append((t,x,d))
    return recs
FIN=2.475202e15; TU=1.202724e-8; A=7.565733250033928e-15; G0=3.98107e5
VMLT=1.86e4
def report(p,tsel=None):
    r=rd(p); t0,x,d0=r[0]
    print('# %s : %d records, t=%.1f..%.1f'%(p,len(r),r[0][0],r[-1][0]))
    for t,x,d in ([r[-1]] if tsel is None else [q for q in r if abs(q[0]-tsel)<1]):
        rho=d[0]; T=d[5]*TU; E=d[8]; F=d[9]; v=d[1]
        Trad=(E/A)**0.25
        print(' t=%.1f  max|drho/rho|=%.3e  max|dT/T|=%.3e  max|dE/E|=%.3e'%(
            t,np.abs(rho/d0[0]-1).max(),np.abs(T/(d0[5]*TU)-1).max(),np.abs(E/d0[8]-1).max()))
        print('   F/Fin: min %.4f max %.4f ; below top 5 cells: min %.4f max %.4f ; top cell %.4f'%(
            (F/FIN).min(),(F/FIN).max(),(F[:-5]/FIN).min(),(F[:-5]/FIN).max(),F[-1]/FIN))
        print('   max|Tgas/Trad-1| = %.3e (at i=%d)'%(np.abs(T/Trad-1).max(),np.argmax(np.abs(T/Trad-1))))
        print('   max|v1|=%.3e cm/s = %.2f v_MLT (i=%d)'%(np.abs(v).max(),np.abs(v).max()/VMLT,np.argmax(np.abs(v))))
        print('   F/Fin profile:', ' '.join('%.3f'%(F[i]/FIN) for i in range(0,84,6)), '| top', '%.3f'%(F[-1]/FIN))
report(sys.argv[1], float(sys.argv[2]) if len(sys.argv)>2 else None)
