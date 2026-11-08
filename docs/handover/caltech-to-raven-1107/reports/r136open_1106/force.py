import sys,os
RD='/resnick/groups/carnegie_poc/jingze/r136open_1106/runs/B'
src=open('/resnick/groups/carnegie_poc/jingze/r136noise_1106/an1.py').read().split('res={}')[0]
n0,n1=int(sys.argv[1]),int(sys.argv[2]); sys.argv=[sys.argv[0]]
exec(src.replace("RUN=R+'/runB/prod'","RUN=%r"%RD))
print('mesh r(Rsun): dv/dt measured [cm/s2] (mean vr), -g, kappa F0/c, -dPg/dr/rho, -dPrad/dr/rho(E/3), net -g+arad-gradPg | mean rho T vr Gamma')
def St(n):
    t,w,m,D=load(n);rho=w['dens'].astype(np.float64);v=w['velx'].astype(np.float64)
    T,Pg=eos.temp_pgas(w['dens'],w['eint']);T=T.astype(np.float64);Pg=Pg.astype(np.float64);E_=m['m1_e'].astype(np.float64);Fl=m['m1_f1'].astype(np.float64)
    kr=kap(w['dens'],T.astype(np.float32)).astype(np.float64);F0=Fl-v*(E_+E_/3)
    return t,rho,v,T,Pg,E_,F0,kr
A=St(n0);B=St(n1);dtt=B[0]-A[0]
t,rho,v,T,Pg,E_,F0,kr=B
dmass=(rho*vol).sum(axis=(0,1));M=GM0/GN+12*(np.cumsum(dmass)-.5*dmass)
g=GN*M/r**2
for rr in (40,50,60,70,80,90,95,98,99.3):
    i=np.argmin(abs(r-rr*RSUN));m=lambda q:q[...,i].mean()
    dPg=np.gradient(Pg.mean(axis=(0,1)),r)[i];dPr=np.gradient((E_/3).mean(axis=(0,1)),r)[i]
    arad=(kr*F0/C)[...,i].mean();ro=m(rho)
    print('%6.1f dvdt %8.1f  g %7.1f arad %7.1f -dPg/rho %8.1f -dPr/rho %8.1f net %8.1f | rho %.2e T %.3g vr %.1f Gam %.3f kap %.3f F0/(L/4pir2) %.3f'%(r[i]/RSUN,(m(B[2])-m(A[2]))/dtt,g[i],arad,-dPg/ro,-dPr/ro,-g[i]+arad-dPg/ro,ro,m(T),m(v)/1e5,arad/g[i],m(kr),m(F0)/(1.92e39*12/(4*np.pi*r[i]**2))))
