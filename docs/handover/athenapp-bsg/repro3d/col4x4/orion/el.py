#!/usr/bin/env python3
"""error_limit comparison of the 4x4 column: KE, energies, profiles, Jacobi stats, cost, Er deficit per solve."""
import numpy as np, h5py, glob, re, os
exec(open('units.py').read())
Eu=Pu*lu**3; GM=260.795*lu**3/tu**2; gam=5/3; D=tu/86400
BASE='/orion/ptmp/jinma/athpp_repro3d/col4x4/'
RUNS=[('1e-6','prof3d'),('1e-7','el1e-7'),('1e-8','el1e-8'),('1e-9','el1e-9')]
WANT=[0.05,0.1,0.25,0.5,1,1.5,2,3]
def ix(f):
    n=[x.decode() for x in f.attrs['VariableNames']];return {v:i for i,v in enumerate(n)}
def shell(g):
    with h5py.File(g,'r') as f:
        I=ix(f);h=f['hydro'];x1=f['x1v'][:];x1f=f['x1f'][:];x2f=f['x2f'][:];x3f=f['x3f'][:];x2=f['x2v'][:];T=f.attrs['Time']*tu/86400
        rs=np.unique(x1);n=len(rs);Q={k:np.zeros(n) for k in('rho','P','Er','F1','F01','rv','PR','S','w','vol')};E=dict(dU=0.,Eg=0.,KE=0.,KEr=0.)
        for b in range(x1.shape[0]):
            ir=np.searchsorted(rs,x1[b]);sn=np.sin(x2[b])
            dO=(np.cos(x2f[b,:-1])-np.cos(x2f[b,1:]))[None,:,None]*(np.diff(x3f[b]))[:,None,None]
            dr3=((x1f[b,1:]**3-x1f[b,:-1]**3)/3)[None,None,:]; dV=dO*dr3
            rho=h[I['rho'],b];P=h[I['press'],b];v1=h[I['vel1'],b];v2=h[I['vel2'],b];v3=h[I['vel3'],b];Er=h[I['Er'],b];F1=h[I['Fr1'],b];F0=h[I['Fr01'],b];S=h[I['Sigma_s_0'],b]+h[I['Sigma_a_0'],b]
            rc=x1[b][None,None,:]
            E['dU']+=np.sum((P/(gam-1))*dV)*Eu+np.sum(Er*dV)*Eru*lu**3
            E['Eg']+=np.sum(-rho*dV/rc)*GM*ru*lu**2
            E['KE']+=np.sum(.5*rho*(v1**2+v2**2+v3**2)*dV)*Eu;E['KEr']+=np.sum(.5*rho*v1**2*dV)*Eu
            wt=np.broadcast_to(np.sin(x2[b])[None,:,None],rho.shape).ravel();idx=np.broadcast_to(ir[None,None,:],rho.shape).ravel()
            for k,q in (('rho',rho),('P',P),('Er',Er),('F1',F1),('F01',F0),('rv',rho*v1),('PR',P/rho),('S',S)):
                Q[k]+=np.bincount(idx,wt*q.ravel(),n)
            Q['w']+=np.bincount(idx,wt,n)
    for k in Q:
        if k not in('w','vol'): Q[k]/=Q['w']
    return T,rs,Q,E
res={}
for name,d in RUNS:
    R=BASE+d+'/';fs=sorted(glob.glob(R+'output/bsg.out2.*.athdf'));od='el%s/'%name;os.makedirs(od,exist_ok=True)
    h=np.loadtxt(R+'output/bsg.hst');t=h[:,0]*D;ke=(h[:,6]+h[:,7]+h[:,8])*Eu
    with open(od+'col4x4_ke.txt','w') as o:
        o.write(f'# 4x4 column, error_limit {name}: t[d] KE KE_r KE_lat [erg, pi/3 sr], hst cadence\n')
        for i in range(len(t)): o.write('%.6e %.6e %.6e %.6e\n'%(t[i],ke[i],h[i,6]*Eu,(h[i,7]+h[i,8])*Eu))
    Ts=np.array([h5py.File(g,'r').attrs['Time']*tu/86400 for g in fs]);S=[];E0=None;prof={};v={}
    with open(od+'col4x4_energy.txt','w') as o:
        o.write(f'# 4x4 column, error_limit {name}: t[d] dE dU dEgrav KE KE_r [erg, pi/3 sr] U=sum(Pg/(g-1)+Er_lab)dV, Egrav=-sum GM rho/r dV, E=U+KE+Egrav\n')
        for g,TT in zip(fs,Ts):
            T,rs,Q,E=shell(g)
            if E0 is None: E0=dict(E)
            dU,dEg=E['dU']-E0['dU'],E['Eg']-E0['Eg'];dE=dU+dEg+E['KE']-E0['KE']
            o.write('%.5f %.6e %.6e %.6e %.6e %.6e\n'%(T,dE,dU,dEg,E['KE'],E['KEr']))
            S.append((T,dE,dU,dEg,E['KE']))
            for w in WANT:
                if abs(T-w)<0.0126 and (w not in prof or abs(T-w)<abs(prof[w][0]-w)):
                    prof[w]=(T,rs,{k:Q[k].copy() for k in Q})
    for w,(T,rs,Q) in sorted(prof.items()):
        Sg=Q['S']/lu;rh=Q['rho']*ru
        with open(od+'col4x4_prof_%g.txt'%w,'w') as o:
            o.write(f'# 4x4 column, error_limit {name}, requested {w} d, dump t={T:.5f} d. Shell means of 256 cells, cgs\n# r[Rsun] rho Pgas Er(lab) Fr1(lab) Fr01(comoving) vr(rho-weighted) T kappa_Ross(<Sigma>/<rho>) Sigma_s+Sigma_a[1/cm]\n')
            for i in range(len(rs)): o.write(' '.join('%.17e'%x for x in (rs[i],rh[i],Q['P'][i]*Pu,Q['Er'][i]*Eru,Q['F1'][i]*Fu,Q['F01'][i]*Fu,Q['rv'][i]/Q['rho'][i]*vu,Q['PR'][i]*Tu,Sg[i]/rh[i],Sg[i]))+'\n')
    log=open(R+'run.log',errors='ignore').read();nit=np.array(list(map(int,re.findall(r'niter: (\d+)',log))))
    cyc=int(re.findall(r'cycle=(\d+)',log)[-1]);cpu=float(re.search(r'cpu time used\s*=\s*([\d.e+-]+)',log).group(1))
    res[name]=dict(t=t,ke=ke,S=np.array(S),prof=prof,nit=nit,cyc=cyc,cpu=cpu)
print('error_limit | Jacobi mean max solves | cycles | wall s (s/cycle)')
for n in res:
    r=res[n];print(n,'%.2f %d %d | %d | %.0f (%.4f)'%(r['nit'].mean(),r['nit'].max(),len(r['nit']),r['cyc'],r['cpu'],r['cpu']/r['cyc']))
print('\nKE [erg] at t, dE dU dEgrav [1e44]')
for T in (0.5,1,1.5,2,2.5,3):
    print('t=%g'%T,end=' ')
    for n in res:
        r=res[n];i=np.argmin(abs(r['t']-T));s=r['S'][np.argmin(abs(r['S'][:,0]-T))];print('| %s KE %.2e dE %.3f dU %.3f dEg %.3f'%(n,r['ke'][i],s[1]/1e44,s[2]/1e44,s[3]/1e44),end=' ')
    print()
print('\nv_r [m/s] at 22,25,27 Rsun, 1 d (Lr: -29,-109,-199)')
for n in res:
    T,rs,Q=res[n]['prof'][1];v=Q['rv']/Q['rho']*vu/100;print(n,[round(float(np.interp(x,rs,v)),1) for x in (22,25,27)])
print('\nEr deficit relative to 1e-9: ln(Er_x/Er_1e-9)/N_steps, shell mean, r bands; dt from cycle count')
for w in (0.05,0.1,0.25):
    T0,rs,Q9=res['1e-9']['prof'][w];N=w*86400/10.57
    for n in ('1e-6','1e-7','1e-8'):
        T,rs,Q=res[n]['prof'][w];x=np.log(Q['Er']/Q9['Er'])/N
        print('t=%g %s per step: r24-35 %.2e  r35-45 %.2e  r20-24 %.2e'%(w,n,np.mean(x[(rs>=24)&(rs<35)]),np.mean(x[(rs>=35)&(rs<45)]),np.mean(x[rs<24])))
