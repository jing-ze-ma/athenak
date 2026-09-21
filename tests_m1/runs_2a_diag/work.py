import numpy as np, struct, sys
CL=2.99792458e10; TU=1.202724e-8; A=7.565733250033928e-15
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
p=sys.argv[1]; t0=float(sys.argv[2]); t1=float(sys.argv[3]); P=float(sys.argv[4])
r=rd(p)
t=np.array([q[0] for q in r]); x=r[0][1]
D=np.array([q[2] for q in r])        # (nt, nvar, nx1)
m=(t>=t0)&(t<=t1); t=t[m]; D=D[m]
rho=D[:,0,:]; v1=D[:,1,:]; T=D[:,5,:]*TU; E=D[:,8,:]; F1=D[:,9,:]; arad=D[:,10,:]
pg=D[:,11,:]
nt,nx=rho.shape
dz=x[1]-x[0]
# optical depth from the TOP, time-averaged
kap=arad*CL/F1                       # kappa_T (cm^2/g)
rk=(rho*kap).mean(axis=0)
tau=np.cumsum((rk*dz)[::-1])[::-1]   # tau at cell i measured from the top
ptot=pg+E/3.0
# ---- cycle-integrated work per unit mass:  W = oint P_tot d(1/rho)
sv=1.0/rho
W=np.zeros(nx)
for i in range(nx):
    W[i]=np.trapezoid(ptot[:,i]*np.gradient(sv[:,i],t),t) if hasattr(np,'trapezoid') \
         else np.trapz(ptot[:,i]*np.gradient(sv[:,i],t),t)
ncyc=(t[-1]-t[0])/P
W/=ncyc
# ---- complex amplitude at the mode frequency (detrended by the growth: use the
#      raw signal; the phase is unaffected by a real exponential envelope)
w=2*np.pi/P
ph=np.exp(-1j*w*t)
def amp(q):
    qq=q-q.mean(axis=0)
    return (qq*ph[:,None]).sum(axis=0)*2.0/len(t)
Av=amp(v1); Ar=amp(rho); Af=amp(F1); Ae=amp(E); Ap=amp(ptot)
print('# i   z/1e8   tau      |dv1|      dv1/cs   |drho/rho|  ph(F1-rho)  ph(P-rho)  W(cycle)  W*rho')
CS=2.1e6
for i in range(nx):
    phf=np.angle(Af[i]/Ar[i])*180/np.pi
    php=np.angle(Ap[i]/Ar[i])*180/np.pi
    print('%3d %7.3f %9.3f %10.3e %8.4f %10.3e %9.1f %9.1f %11.3e %11.3e'%(
        i,x[i]/1e8,tau[i],abs(Av[i]),abs(Av[i])/CS,abs(Ar[i])/rho[:,i].mean(),
        phf,php,W[i],W[i]*rho[:,i].mean()))
print('# total sum W*rho*dz = %.4e ; positive cells (driving) sum = %.4e'%(
    (W*rho.mean(axis=0)*dz).sum(), (np.clip(W,0,None)*rho.mean(axis=0)*dz).sum()))
