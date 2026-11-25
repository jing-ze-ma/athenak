"""N2 profile and scales from shell means of rho, e_gas (hydro_w) and E (vet) of one dump; usage prof2.py RUN"""
import sys, re, glob, numpy as np
sys.argv0=sys.argv[:]; run=sys.argv[1]; sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d'); import bin_convert as bc
J='/resnick/groups/carnegie_poc/jingze/'; A=7.565733e-15; kB=1.380649e-16; c_=2.99792458e10
GM=dict(m35zams=4.6449e27,t9l62=7.4033817e27,t9mlt=7.4033817e27,bsg=2.62502e27,m35mms=4.6449e27); MU=dict(m35zams=.6220,t9l62=.674,t9mlt=.674,bsg=.62,m35mms=.6215)
DIR=dict(m35zams=('prod_c1_1124/m35zams','m35zams'),t9l62=('prod_c4_1124/t9l62','t9l62'),t9mlt=('prod_mlt_1125/t9mlt','t9l62'),bsg=('prod_c1_1124/bsg','bsg3d'),m35mms=('prod_c1_1124/m35mms','m35mms'))
OP=dict(bsg=J+'fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt'); OPD=J+'fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
def kap(rho,T,f):
    g=[l for l in open(f) if re.match(r'#\s*\d',l)][0][1:].split(); nT,nD=int(g[0]),int(g[1]); lT0,dlT,lD0,dlD=map(float,g[2:6]); tab=np.loadtxt(f,comments='#').reshape(nT,nD)
    x=np.clip((np.log10(T)-lT0)/dlT,0,nT-1-1e-6); y=np.clip((np.log10(rho)-lD0)/dlD,0,nD-1-1e-6); i=x.astype(int); j=y.astype(int); fx=x-i; fy=y-j
    return 10**((1-fx)*(1-fy)*tab[i,j]+fx*(1-fy)*tab[i+1,j]+(1-fx)*fy*tab[i,j+1]+fx*fy*tab[i+1,j+1])
d,b=DIR[run]; base=b
S=sorted(glob.glob('stars/%s/s_*.npz'%run)); z=np.load(S[len(S)//2]); n=int(z['n']); f=J+d+'/bin/%s.vet.%05d.bin'%(base,n)
m=bc.read_binary(f); lg=np.asarray(m['mb_logical']); ni=len(z['r']); n3,n2,n1=m['mb_data']['vet_e'][0].shape; P=np.empty((n3*(lg[:,2].max()+1),n2*(lg[:,1].max()+1),ni))
for q in range(len(lg)): P[lg[q,2]*n3:(lg[q,2]+1)*n3,lg[q,1]*n2:(lg[q,1]+1)*n2]=np.asarray(m['mb_data']['vet_e'][q])[:,:,:ni]
ex=z['ext']; N3,N2=z['N']; th=np.linspace(ex[0],ex[1],N2+1); w=np.cos(th[:-1])-np.cos(th[1:]); w=w/w.sum()/N3; E=np.einsum('kji,j->i',P,w)
r=z['r']; rho=z['rho']; T=(2/3.)*z['ei']/rho*MU[run]*1.66053906660e-24/kB; p=(2/3.)*z['ei']; pt=p+E/3; be=p/pt; G=GM[run]/r**2; g3m1=((4-3*be)*(2/3.))/(be+12*(2/3.)*(1-be)); gam1=be+(4-3*be)**2*(2/3.)/(be+12*(2/3.)*(1-be))
d1=lambda q: np.gradient(np.log(q),r); N2_=G*(d1(pt)/gam1-d1(rho)); Hp=-pt/np.gradient(pt,r); cs=np.sqrt(gam1*pt/rho); k=kap(rho,T,OP.get(run,OPD)); cp=(kB/(MU[run]*1.66053906660e-24))/be*(32-24*be-3*be**2)/(2*be); chi=4*A*c_*T**3/(3*k*rho)/(rho*cp)
np.savez('stars/prof_%s.npz'%run,r=r,rho=rho,T=T,N2=N2_,Hp=Hp,cs=cs,kap=k,beta=be,chi=chi,tau=np.cumsum((k*rho*np.gradient(r))[::-1])[::-1],n=n,t=float(z['t']))
print(run,'dump',n,'t %.3e'%float(z['t']),'cells',ni,'ext',ex,'N',z['N']); 
for i in range(0,ni,max(1,ni//14)): print('  r %.4e (%.3f of r0) N2 %+.2e rho %.2e T %.2e beta %.2f tau>.. Hp/r %.3f dr/Hp %.3f'%(r[i],r[i]/r[0],N2_[i],rho[i],T[i],be[i],Hp[i]/r[i],np.gradient(r)[i]/Hp[i]))
