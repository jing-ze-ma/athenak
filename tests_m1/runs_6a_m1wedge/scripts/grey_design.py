# design of the grey point-mass atmosphere (gate A): integrate from the top inward
import numpy as np, sys
def run(gm, c, a, kt, fin, rin, rout, rho_top, q=0.5, n=20000):
    r = np.linspace(rout, rin, n); dr = r[1]-r[0]
    F = lambda x: fin*rin**2/x**2
    E = F(rout)/(c*q); p = rho_top*(E/a)**0.25
    out=[]
    def rhs(x, y):
        E,p = y; T=(E/a)**0.25; rho=p/T
        return np.array([-3*rho*kt*F(x)/c, -rho*(gm/x**2 - kt*F(x)/c)])
    y=np.array([E,p]); tau=0.0
    for i,x in enumerate(r):
        E,p=y; T=(E/a)**0.25; rho=p/T
        out.append((x,rho,T,E,p,tau))
        if i==n-1: break
        k1=rhs(x,y);k2=rhs(x+dr/2,y+dr/2*k1);k3=rhs(x+dr/2,y+dr/2*k2);k4=rhs(x+dr,y+dr*k3)
        y=y+dr/6*(k1+2*k2+2*k3+k4); tau+= -rho*kt*dr
    return np.array(out)
gm,c,a,kt,fin,rin,rout,rt=[float(v) for v in sys.argv[1:9]]
o=run(gm,c,a,kt,fin,rin,rout,rt)
Gam=kt*fin/(c*gm)
print("Gamma",Gam)
for f in [0,0.02,0.05,0.1,0.2,0.4,0.6,0.8,1.0]:
    i=min(int(f*(len(o)-1)),len(o)-1); x,rho,T,E,p,tau=o[i]
    H=T*x**2/(gm*(1-Gam))
    print("r %.4f rho %.3e T %.3e cs %.3e c/cs %.0f Prad/Pg %.3f tau %.3e H %.4f"%(x,rho,T,np.sqrt(5/3*T),c/np.sqrt(5/3*T),E/3/p,tau,H))
