#!/usr/bin/env python3
"""Full radial shell means from the dense 0-3 d rerun; check it reproduces the original column."""
import numpy as np, h5py, glob, re, sys
exec(open('units.py').read())
RUN='/orion/ptmp/jinma/athpp_repro3d/col4x4/'
fs=sorted(glob.glob(RUN+'prof3d/output/bsg.out2.*.athdf')); print('dumps',len(fs))
f0=h5py.File(fs[0],'r');print(dict(f0.attrs)['Time'],list(f0.keys())[:12])
names=[n.decode() for n in f0.attrs['VariableNames']];print(names)
ds=list(f0.attrs['DatasetNames']);nv=list(f0.attrs['NumVariables']);print([d.decode() for d in ds],nv)
ix={};k=0
for d,n in zip(ds,nv):
    for i in range(n): ix[names[k]]=(d.decode(),i);k+=1
def load(f,v): d,i=ix[v];return f[d][i]
tt=[];
for g in fs:
    with h5py.File(g,'r') as f: tt.append(f.attrs['Time']*tu/86400)
tt=np.array(tt)
want=[0,.05,.1,.25,.375,.5,.75,1,1.25,1.5,2,2.5,3]
for T in want:
    j=int(np.argmin(abs(tt-T)))
    with h5py.File(fs[j],'r') as f:
        x1=f['x1v'][:];x2=f['x2v'][:];x3=f['x3v'][:];nb=x1.shape[0]
        rs=np.unique(x1);ii=[];w=[]
        rho=load(f,'rho');P=load(f,'press');v=load(f,'vel1');Er=load(f,'Er');F1=load(f,'Fr1');F01=load(f,'Fr01')
        Ss=load(f,'Sigma_s_0') if 'Sigma_s_0' in ix else load(f,'Sigma_s');Sa=load(f,'Sigma_a_0') if 'Sigma_a_0' in ix else load(f,'Sigma_a')
        out=np.zeros((len(rs),9));mw=np.zeros(len(rs))
        for b in range(nb):
            th=np.sin(x2[b])[None,:,None];wt=np.broadcast_to(th,rho[b].shape)
            ir=np.searchsorted(rs,x1[b])
            for kk,q in enumerate((rho[b],P[b],Er[b],F1[b],F01[b],rho[b]*v[b],P[b]/rho[b],Ss[b],Sa[b])):
                out[:,kk]+=np.bincount(np.broadcast_to(ir[None,None,:],rho[b].shape).ravel(),(wt*q).ravel(),len(rs))
            mw+=np.bincount(np.broadcast_to(ir[None,None,:],rho[b].shape).ravel(),wt.ravel(),len(rs))
    m=out/mw[:,None]; r0,P0,E0,F,F0,rv,PR,sS,sA=m.T
    Sig=(sS+sA)/lu  # code Sigma per code length -> cm^-1
    rhoc=r0*ru
    name='col4x4_prof_%g.txt'%T
    with open(name,'w') as o:
        o.write(f'# Orion 4x4 column rerun, requested {T} d, dump at t={tt[j]:.5f} d. Shell means (solid-angle weighted) of all 256 cells, cgs\n')
        o.write('# r[Rsun] rho Pgas Er(lab) Fr1(lab) Fr01(comoving) vr(rho-weighted) T kappa_Ross(shell=<Sigma_s+Sigma_a>/<rho>) Sigma_s+Sigma_a[1/cm]\n')
        for i in range(len(rs)):
            o.write(' '.join('%.17e'%x for x in (rs[i],rhoc[i],P0[i]*Pu,E0[i]*Eru,F[i]*Fu,F0[i]*Fu,rv[i]/r0[i]*vu,PR[i]*Tu,Sig[i]/rhoc[i],Sig[i]))+'\n')
    print(name,'t=%.4f'%tt[j],'kappa cell0 %.3f, cell127 %.3f, v(26.41) %.1f m/s'%(Sig[0]/rhoc[0],Sig[127]/rhoc[127],np.interp(26.41,rs,rv/r0*vu)/100))
# reproduce check: hst vs original, 0.5 d dump vs original tab
a=np.loadtxt(RUN+'run/output/bsg.hst');b=np.loadtxt(RUN+'prof3d/output/bsg.hst');n=min(len(a),len(b));print('hst rows',len(b),'max rel diff cols 0-11:',np.nanmax(abs(a[:n,:12]-b[:n,:12])/(abs(a[:n,:12])+1e-30),axis=0).round(9))
