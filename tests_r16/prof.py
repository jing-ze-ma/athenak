import numpy as np, struct, sys
def rd(fn):
    out=[];f=open(fn,'rb')
    while True:
        h=f.read(16)
        if len(h)<16: break
        t,n1,nv=struct.unpack('<dii',h)
        b=f.read(8*(n1+nv*n1))
        if len(b)<8*(n1+nv*n1): break
        a=np.frombuffer(b,'<f8'); out.append((t,a[:n1].copy(),a[n1:].reshape(nv,n1).copy()))
    f.close(); return out
RS=2.3717e11; GM=6.674e-8*1.989e33*4.0   # placeholder mass; g computed from file? use GM below
arms=sys.argv[1:]
P={a:rd(a+'/rt_profile.bin') for a in arms}
for a in arms: print(a,len(P[a]),'t=%.2f..%.2f'%(P[a][0][0],P[a][-1][0]))
# per-cell net force = d(rho v1)/dt between consecutive records, normalised by rho*g
def force(P,i0,i1):
    t0,r,q0=P[i0]; t1,_,q1=P[i1]
    return r,(q1[2]-q0[2])/(t1-t0), 0.5*(q0[0]+q1[0])
for (t_lo,t_hi) in [(0,20),(80,100),(280,300)]:
    print('\n=== d(rho v1)/dt averaged over t=%g..%g s  [g/cm^2/s^2]'%(t_lo,t_hi))
    res={}
    for a in arms:
        ts=np.array([p[0] for p in P[a]])
        i0=int(np.argmin(abs(ts-t_lo))); i1=int(np.argmin(abs(ts-t_hi)))
        r,f,rho=force(P[a],i0,i1); res[a]=(r,f,rho)
    r=res[arms[0]][0]
    g=GM/r**2
    hdr='   i  r/R   '+ ''.join('%12s'%a for a in arms) + '%13s'%'diff' + '%10s'%'d/(rho g)'
    print(hdr)
    fs=[res[a][1] for a in arms]; rho=res[arms[0]][2]
    d=fs[0]-fs[1] if len(fs)>1 else fs[0]
    key=np.argsort(-np.abs(d))[:14]
    for i in sorted(key):
        print('%4d %6.4f '%(i,r[i]/RS)+''.join('%12.4e'%f[i] for f in fs)
              +'%13.4e'%d[i]+'%10.3f'%(d[i]/(rho[i]*g[i])))
    # integrated
    print(' integrated sum over shells of (f_a - f_b)*... [per-cell, unweighted]: %.4e'%d.sum())
