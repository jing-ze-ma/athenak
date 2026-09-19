import struct, sys, numpy as np
def rd(fn):
    out=[];f=open(fn,'rb')
    while True:
        h=f.read(16)
        if len(h)<16: break
        t,n1,nv=struct.unpack('<dii',h); b=f.read(8*(n1+nv*n1))
        if len(b)<8*(n1+nv*n1): break
        a=np.frombuffer(b,'<f8'); out.append((t,a[:n1].copy(),a[n1:].reshape(nv,n1).copy()))
    return out
RS=2.3717e11; GM=6.674e-8*6.2651e33
for arm in sys.argv[1:]:
    R=rd(arm+'/rt_profile.bin'); r=R[0][1]
    rf=np.concatenate(([1.18585e11],0.5*(r[1:]+r[:-1]),[2.4057e11])); vol=4*np.pi/3*(rf[1:]**3-rf[:-1]**3)
    a=int(np.argmin([abs(x[0]-130) for x in R])); b=len(R)-1
    f=(R[b][2][2]-R[a][2][2])/(R[b][0]-R[a][0]); g=R[a][2][0]*GM/r**2; F=f*vol
    print("== %s  t %.0f-%.0f s  total %.3e dyn | 0.90-0.96R %.3e | rest %.3e"%(arm,R[a][0],R[b][0],F.sum(),F[(r>=0.9*RS)&(r<0.96*RS)].sum(),F[(r<0.9*RS)|(r>=0.96*RS)].sum()))
    print("   f/(rho g) at r/R:", " ".join("%.3f:%+.4f"%(r[k]/RS,f[k]/g[k]) for k in range(len(r)) if 0.89<r[k]/RS<0.96))
