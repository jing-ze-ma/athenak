"""wbforce.py n0 n1 [run]: predicted x1 momentum source of the pgen's hs_grav (walked pressures, code's polytropic stencil) above wb_rmax vs measured dv/dt.
Offline, ideal-gas chi's (chi_rho=chi_T=1, p/(rho T) from the table per cell), Phi_eff' = g - a_ref with a_ref taken = shell mean kappa F/c (phi=1)."""
import sys,os
RUNN=sys.argv[3] if len(sys.argv)>3 else 'B'
RD='/resnick/groups/carnegie_poc/jingze/r136open_1106/runs/'+RUNN
src=open('/resnick/groups/carnegie_poc/jingze/r136noise_1106/an1.py').read().split('res={}')[0]
n0,n1=int(sys.argv[1]),int(sys.argv[2]); sys.argv=[sys.argv[0]]
exec(src.replace("RUN=R+'/runB/prod'","RUN=%r"%RD))
def St(n):
    t,w,m,D=load(n);rho=w['dens'].astype(np.float64);v=w['velx'].astype(np.float64)
    T,Pg=eos.temp_pgas(w['dens'],w['eint']);T=T.astype(np.float64);Pg=Pg.astype(np.float64);E_=m['m1_e'].astype(np.float64);Fl=m['m1_f1'].astype(np.float64)
    kr=kap(w['dens'],T.astype(np.float32)).astype(np.float64);F0=Fl-v*(E_+E_/3)
    return dict(t=t,rho=rho,v=v,T=T,P=Pg,F0=F0,kr=kr)
A=St(n0);B=St(n1);dtt=B['t']-A['t']
Am,Ap=x1f[:-1]**2,x1f[1:]**2   # area/(dOmega)
V=dvr
def src_acc(S):
    rho,T,P=S['rho'],S['T'],S['P'];c=P/(rho*T)
    dmass=(rho*vol).sum(axis=(0,1));M=GM0/GN+12*(np.cumsum(dmass)-.5*dmass)
    aref=(S['kr']*S['F0']/C).mean(axis=(0,1))
    # Phi_eff on cells and faces: integrate g - a_ref (piecewise constant between centres)
    g_eff=GN*M/r**2-aref
    PhiC=np.zeros_like(r);PhiF=np.zeros_like(x1f)
    PhiF[0]=0
    for i in range(len(r)):
        PhiC[i]=PhiF[i]+g_eff[i]*(r[i]-x1f[i]);PhiF[i+1]=PhiC[i]+g_eff[i]*(x1f[i+1]-r[i])
    # local polytrope a = (T_ip1-T_im1)/(Phi_ip1-Phi_im1), edge cells one-sided
    Tm=np.concatenate([T[...,:1],T[...,:-1]],axis=2);Tp=np.concatenate([T[...,1:],T[...,-1:]],axis=2)
    Pm=np.concatenate([PhiC[:1],PhiC[:-1]]);Pp=np.concatenate([PhiC[1:],PhiC[-1:]])
    a=(Tp-Tm)/(Pp-Pm)[None,None,:]
    def seg(d0,t0,k0,dph,cc):
        dpred=d0*np.exp(k0*dph);t1=t0+a*dph
        tk=np.where(t1>0,t1,np.nan)
        kend=-(1/(cc*tk)+a/tk)
        d1=d0*np.exp(.5*(k0+kend)*dph);p1=d1*cc*tk
        return d1,tk,p1,-(1/(cc*tk)+a/tk)
    k0=-(1/(c*T)+a/T)
    dmh=(PhiF[:-1]-PhiC)[None,None,:];dph=(PhiF[1:]-PhiC)[None,None,:]
    d1,t1,pl,k1=seg(rho,T,k0,dmh,c);d3,t3,pr,k3=seg(rho,T,k0,dph,c)
    # outer half cells for the guard
    dm1=(Pm-PhiF[:-1])[None,None,:];dp1=(Pp-PhiF[1:])[None,None,:]
    d2,t2,p2,_=seg(d1,t1,k1,dm1,c);d4,t4,p4,_=seg(d3,t3,k3,dp1,c)
    bad=lambda q,q0:~np.isfinite(q)|(q<=0)|(q<1e-6*q0)|(q>1e6*q0)
    gd=bad(pl,P)|bad(pr,P)|bad(p2,P)|bad(p4,P)|bad(d1,rho)|bad(d2,rho)|bad(d3,rho)|bad(d4,rho)
    pl2=np.where(gd,P,pl);pr2=np.where(gd,P,pr)
    acc=(Ap[None,None,:]*(pr2-P)+Am[None,None,:]*(P-pl2))/V[None,None,:]/rho
    plain=-GN*M/r**2+aref
    return acc,plain,gd,aref,g_eff,PhiC,c
acc0,plain,gd0,aref,ge,_,c0=src_acc(A);acc1,_,gd1,_,_,_,_=src_acc(B)
acc=.5*(acc0+acc1)
print('r  dr  dvdt_meas  a_src(mean over cells) a_src(median) -g_eff  ratio(mean a_src/g_eff) guard-frac H/dr(g_eff) arad-aref(shell) dPg/rho')
for rr in (28,30,32,35,38,40,45,50,55,60,65,70,75,80,85,90,95,98,99.3):
    i=np.argmin(abs(r-rr*RSUN));H=(c0*A['T'])[...,i].mean()/ge[i]
    print('%6.1f %6.3f dvdt %9.1f  a_src %9.1f med %9.1f  -geff %8.1f  ratio %7.2f guard %.3f H/dr %.3f sum %.1f'%(r[i]/RSUN,dr[i]/RSUN,(B['v'][...,i].mean()-A['v'][...,i].mean())/dtt,np.nanmean(acc[...,i]),np.nanmedian(acc[...,i]),-ge[i],np.nanmean(acc[...,i])/-ge[i],gd0[...,i].mean(),H/dr[i],np.nanmean(acc[...,i])))
