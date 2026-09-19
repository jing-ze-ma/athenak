import numpy as np
RS=2.3717e11; r0=1.18585e11
def u_of(c,xi): return xi+(1-xi)*sum(c[k]*xi**(k+1) for k in range(len(c)))
cold=[0.134559,4.942505,-8.950106,4.403042]; r1o=2.4057e11; n=96
xf=np.linspace(0,1,n+1); ro=r0+(r1o-r0)*u_of(cold,xf); dro=np.diff(ro)
print("OLD grid: dr [1e8 cm] at r/R:"," ".join("%.3f:%.1f"%(0.5*(ro[i]+ro[i+1])/RS,dro[i]/1e8) for i in range(0,n,8)),"| last %.1f"%(dro[-1]/1e8))
# target dr(r): as the old grid up to 0.97R, 6e8 over 0.97-1.10R, then geometric growth to r1
import sys
r1=float(sys.argv[1])*RS; N=int(sys.argv[2])
def dr_target(r):
    x=r/RS
    old=np.interp(r,0.5*(ro[1:]+ro[:-1]),dro)
    if x<0.97: return old
    if x<1.10: return 6.0e8
    return 6.0e8*np.exp((x-1.10)/0.045)
# build faces by marching, then rescale cell count to N by uniform factor on dr
faces=[r0]
while faces[-1]<r1: faces.append(faces[-1]+dr_target(faces[-1]))
faces=np.array(faces); print("target spacing needs %d cells"%(len(faces)-1))
# map: cell index fraction -> radius ; resample to N cells
s=np.linspace(0,1,len(faces)); faces[-1]=r1
xi=np.linspace(0,1,2001); ut=(np.interp(xi,s,faces)-r0)/(r1-r0)
# least squares for c: ut - xi = (1-xi) sum c_k xi^k
A=np.stack([(1-xi)*xi**(k+1) for k in range(4)],1); c,res,_,_=np.linalg.lstsq(A,ut-xi,rcond=None)
xfN=np.linspace(0,1,N+1); rn=r0+(r1-r0)*u_of(c,xfN); drn=np.diff(rn)
print("c =",["%+.6f"%v for v in c],"monotonic:",bool((drn>0).all()),"min dr %.2e max dr %.2e"%(drn.min(),drn.max()))
print("NEW grid N=%d: dr [1e8] at r/R:"%N," ".join("%.3f:%.1f"%(0.5*(rn[i]+rn[i+1])/RS,drn[i]/1e8) for i in range(0,N,8)),"| last %.1f"%(drn[-1]/1e8))
m=(rn[:-1]/RS>0.97)&(rn[:-1]/RS<1.10); print("cells in 0.97-1.10 R: %d, dr there %.1f-%.1f e8"%(m.sum(),drn[m].min()/1e8,drn[m].max()/1e8))
