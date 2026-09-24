import numpy as np
def tab(g2,g3):
    return np.array([[0,0,0],[1-g2,g2,0],[0.5,0.5-g3,g3]])
def R(A,z):
    # stage amplifications for y'=l y, z=dt l; explicit first stage
    Y=[1.0+0j]
    for i in range(1,3):
        s=1+z*sum(A[i,j]*Y[j] for j in range(i))
        Y.append(s/(1-z*A[i,i]))
    return Y
def astab(A):
    ys=np.linspace(-200,200,40001)
    m=max(abs(R(A,1j*y)[2]) for y in ys)
    return m
def pr(A,lam,dts,T=1.0):
    # Prothero-Robinson y'=lam(y-phi)+phi', phi=sin t ; errors at T
    errs=[]
    for dt in dts:
        n=int(round(T/dt)); y=0.0; t=0.0
        f=lambda t,y: lam*(y-np.sin(t))+np.cos(t)
        c=[0,1,1]
        for k in range(n):
            K=[f(t,y)]; Ys=[y]
            for i in range(1,3):
                s=y+dt*sum(A[i,j]*K[j] for j in range(i)); ti=t+c[i]*dt
                g=A[i,i]; Yi=(s+g*dt*(-lam*np.sin(ti)+np.cos(ti)))/(1-g*dt*lam)
                K.append(f(ti,Yi)); Ys.append(Yi)
            y=Ys[2]; t+=dt
        errs.append(abs(y-np.sin(T)))
    e=np.array(errs); return e, np.log2(e[:-1]/e[1:])
g=1-1/np.sqrt(2)
for g3 in [g,0.25,0.2,0.15,0.1]:
    g2=(1-2*g3)/(2*(1-g3))
    A=tab(g2,g3)
    r2=-(1-g2)/g2
    z=-1e8; Rinf=R(A,z)[2].real
    # error constant: y''' term: sum b_i c_i^2/2 - 1/6 and b A c -1/6
    b=A[2]; c=A.sum(1)
    ec=(abs(b@c**2/2-1/6)+abs(b@A@c-1/6))
    e,o=pr(A,-1e3,[0.1,0.05,0.025,0.0125,0.00625])
    print(f"g3={g3:.4f} g2={g2:.4f} stageRinf={r2:.3f} Rinf={Rinf:.2e} max|R(iy)|={astab(A):.4f} errc={ec:.4f} PR(-1e3) err {e[0]:.2e} ord {np.round(o,2)}")
# alt: stage BE (g2=1)
for g3 in [0.5,1.0]:
    A=tab(1.0,g3); print('g2=1 g3',g3,'Rinf',R(A,-1e8)[2].real,'maxRiy',astab(A))
print("--- dissipation/dispersion on imaginary axis (radwave-like) and PR at several stiffness")
for g3 in [g,0.2,0.15]:
    g2=(1-2*g3)/(2*(1-g3)); A=tab(g2,g3)
    out=[]
    for y in [0.25,0.5,1.0,2.0]:
        Rv=R(A,1j*y)[2]; out.append(f"y={y}: |R|-1={abs(Rv)-1:+.2e} ph/y-1={np.angle(Rv)/y-1:+.2e}")
    print(f"g3={g3:.4f}", "; ".join(out))
    for lam in [-10,-1e2,-1e4,-1e6]:
        e,o=pr(A,lam,[0.1,0.05,0.025,0.0125])
        print(f"   PR lam={lam:g} err {e[0]:.2e} {e[-1]:.2e} ord {np.round(o,2)}")
    # real negative axis: stage and final amplification
    print("   z: stage/final", [(zz, round(R(A,zz)[1].real,3), round(R(A,zz)[2].real,3)) for zz in [-1,-3,-10,-100]])
print("--- SDIRK2 (g, not Heun-pairable: b^T c_explicit = g != 1/2) PR for reference")
def pr_sdirk(lam,dts,T=1.0):
    errs=[]; gg=g
    for dt in dts:
        n=int(round(T/dt)); y=0.0; t=0.0
        f=lambda t,y: lam*(y-np.sin(t))+np.cos(t)
        for k in range(n):
            t1=t+gg*dt; Y1=(y+gg*dt*(-lam*np.sin(t1)+np.cos(t1)))/(1-gg*dt*lam); K1=f(t1,Y1)
            t2=t+dt; s=y+dt*(1-gg)*K1; Y2=(s+gg*dt*(-lam*np.sin(t2)+np.cos(t2)))/(1-gg*dt*lam)
            y=Y2; t+=dt
        errs.append(abs(y-np.sin(T)))
    e=np.array(errs); return e, np.log2(e[:-1]/e[1:])
for lam in [-10,-1e2,-1e4,-1e6]:
    e,o=pr_sdirk(lam,[0.1,0.05,0.025,0.0125]); print(f"   SDIRK2 PR lam={lam:g} err {e[0]:.2e} {e[-1]:.2e} ord {np.round(o,2)}")
