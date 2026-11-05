import numpy as np, glob, sys, os
RS=6.957e10; ru=1e-9; Tu=1e5; prat=56.41681; GM=2.84560241e+03; c=2.99792458e10; a=7.565733250033928e-15
Fu=c*a*Tu**4; Lin=5.581461e38
for d in sys.argv[1:]:
    fs=sorted(glob.glob(d+'/out/bsg.block0.out2.*.tab'))
    print('==',d, 'files',len(fs))
    print(' t[d]  r_ph(tau=2/3)  r(rho=1e-10)  L_top/L  Etot(code)  Egas  Erad  Egrav  M(r>45)  vmax[km/s]  rho_min(40-50)')
    for f in fs:
        t=float(open(f).readline().split('time=')[1].split()[0])
        x=np.loadtxt(f,comments='#'); r=x[:,1]; rho=x[:,2]; p=x[:,3]; v=x[:,4]; E=x[:,7]; F=x[:,8]
        chi=x[:,-3]+x[:,-2]
        dr=np.gradient(r); dV=4*np.pi*r**2*dr/3.0   # wedge-free (per steradian x 4pi/3 irrelevant)
        tau=np.cumsum((chi*dr)[::-1])[::-1]
        iph=np.where(tau>2/3)[0]; rph=r[iph[-1]] if len(iph) else r[0]
        isu=np.where(rho*ru>1e-10)[0]; rsu=r[isu[-1]] if len(isu) else r[0]
        Eg=np.sum(p/(2/3)*dV); Er=np.sum(prat*E*dV); Ek=np.sum(0.5*rho*v**2*dV); Egr=np.sum(-GM*rho/r*dV)
        L=4*np.pi*(r[-1]*RS)**2*F[-1]*Fu/Lin
        m=np.sum((rho*dV)[r>45])
        sel=(r>40)&(r<50)
        vk=np.max(np.abs(v))*np.sqrt(8.314462618e7/0.62*Tu)/1e5
        print('%5.2f %7.3f %7.3f %7.3f %.6e %.4e %.4e %.4e %.4e %7.2f %.3e'%(t/4.547922,rph,rsu,L,Eg+Er+Ek+Egr,Eg,Er,Egr,m,vk,np.min(rho[sel])*ru))
