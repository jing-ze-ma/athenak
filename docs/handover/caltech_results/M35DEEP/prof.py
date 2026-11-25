import numpy as np
RP=9.554e11; GM=4.6449e27; A=7.565733e-15; kB=1.380649e-16; mu=1.66053906660e-24*0.6215; c=2.99792458e10
z=np.load('/resnick/groups/carnegie_poc/jingze/pages_c1_1124/m35mms/out/data/red_00300.npz'); r=z['r']; rho=z['rho']; T=z['Tgas']; E=z['E']; x=r/RP
P=rho*kB*T/mu; Pr=E/3; Pt=P+Pr; b=P/Pt
G=GM/r**2; 
d=lambda f: np.gradient(np.log(f),r)
Gam1=b+ (4-3*b)**2*(5/3-1)/(b+12*(5/3-1)*(1-b))
N2=-G*(np.gradient(np.log(rho),r)-np.gradient(np.log(Pt),r)/Gam1); Hp=-Pt/np.gradient(Pt,r)
cs=np.sqrt(Gam1*Pt/rho); ca=np.sqrt(5/3*P/rho); fac=cs/(2*Hp)   # acoustic cutoff c/2Hp
# gradients
nab=d(T)/d(Pt)
nabad=np.gradient(np.log(T),np.log(Pt))   # placeholder
nad=(b**2*0+0)
nad=(1-b)*(4-3*b); 
# adiabatic gradient for gas+radiation (gamma 5/3): (Gam2-1)/Gam2
num=(4-3*b)*(5/3-1); den=b+12*(5/3-1)*(1-b)
g3m1=num/den; gam2=1/(1-g3m1/Gam1)
nad=g3m1/Gam1
import re
def kap(rho,T,OP='/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'):
    g=[l for l in open(OP) if re.match(r'#\s*\d',l)][0][1:].split(); nT,nD=int(g[0]),int(g[1]); lT0,dlT,lD0,dlD=map(float,g[2:6]); tab=np.loadtxt(OP,comments='#').reshape(nT,nD)
    xx=np.clip((np.log10(T)-lT0)/dlT,0,nT-1-1e-6); yy=np.clip((np.log10(rho)-lD0)/dlD,0,nD-1-1e-6); i=xx.astype(int); j=yy.astype(int); fx=xx-i; fy=yy-j
    return 10**((1-fx)*(1-fy)*tab[i,j]+fx*(1-fy)*tab[i+1,j]+(1-fx)*fy*tab[i,j+1]+fx*fy*tab[i+1,j+1])
k=kap(rho,T)
L=3.828e33*10**5.444
nrad=3*k*L*Pt/(16*np.pi*A*c*GM*T**4)
# taulocal depth
tau=np.cumsum((k*rho*np.gradient(r))[::-1])[::-1]
eps=1e-3
kT=np.log(kap(rho,T*(1+eps))/kap(rho,T*(1-eps)))/(2*eps); kr=np.log(kap(rho*(1+eps),T)/kap(rho*(1-eps),T))/(2*eps)
omc=lambda l: None
vth=4*A*T**3/(3*0)if False else None
Kth=4*A*c*T**3/(3*k*rho) ; cP=(5/2*kB/mu)*(b*0+1)+ 4*4*A*T**3/(3*rho)
# thermal diffusivity (rad conduction) chi = K/(rho c_p); c_p gas+rad: use kB/mu*(2.5 + (4-3b)^2/b^2...) 
cp=(kB/mu)*( (2.5*b**2*0) + 1/b*((32-24*b-3*b**2)/(2*b)) )
chi=Kth/(rho*cp); Ttherm=Hp**2/chi
print('  x=r/Rp     rho       T       beta   N2[1/s2]   N[1/s]   c_s[km/s]  Hp/r   dz/Hp  w_ac[1/s]  nab   nad  nrad  kappa  dlnk/dlnT dlnk/dlnrho tau    Hp^2/chi[s]  g  Gam1')
dz=np.gradient(r)
for xx in (0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.80,0.85,0.90):
    i=int(np.argmin(abs(x-xx))); print('%6.3f %.2e %.3e %.3f %+.2e %+.2e %6.1f %.4f %.3f %.2e %.3f %.3f %.3f %.3f %+.2f %+.2f %.2e %.2e %.2e %.3f'%(x[i],rho[i],T[i],b[i],N2[i],np.sign(N2[i])*np.sqrt(abs(N2[i])),cs[i]/1e5,Hp[i]/r[i],dz[i]/Hp[i],cs[i]/(2*Hp[i]),nab[i],nad[i],nrad[i],k[i],kT[i],kr[i],tau[i],Ttherm[i],G[i],Gam1[i]))
np.savez('prof.npz',x=x,N2=N2,cs=cs,Hp=Hp,rho=rho,T=T,kap=k,nab=nab,nad=nad,nrad=nrad,tau=tau,kT=kT,kr=kr,r=r,beta=b,chi=chi)
