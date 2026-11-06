import numpy as np,glob,os
Rs=6.957e10;vu=np.sqrt(8.314462618e7*1e5/.62);tu=Rs/vu;Eu=1e-9*vu**2*Rs**3;LS=5.581e38;F=12*Eu
Pu=1e-9*vu**2;Eru=7.5657e-15*1e20;W='/viper/ptmp2/jinma/athpp_ke_1106/wall_1106/'
idx=dict(E=2,hyd=3,grav=4,pm=5,bill=6,cons=7,clip=8,reset=9,fE=10,fM=11,radD=12,gas=13)
days=[0,.5,1,2,3,5,7,10]
for a in ['A_af1','A_tc4']:
    print('=====',a)
    fs=sorted(glob.glob(f'{a}/out/ledger_prof_*.txt')); P=[]
    for f in fs:
        d=np.loadtxt(f)
        if not np.isnan(d[:40]).any(): P.append((float(os.path.basename(f)[12:-4]),d))
    h=np.loadtxt(f'{a}/out/bsg.hst')
    D=sorted(glob.glob(f'{W}{a}/out/bsg.block0.out2.*.tab'))
    print('cell budget per interval in L* (full sphere). B=billing, C=consistent rad transport, H=hydro flux+pointmass, G=grav mass x Phi, R=reset+floors, dU=dE-G')
    for (t0,d0),(t1,d1) in zip(P[:-1],P[1:]):
        dt=(t1-t0)*tu; r=(d1-d0)*F/dt/LS; s=np.median(h[(h[:,0]>=t0)&(h[:,0]<t1),1])*tu
        line=f'{t0*tu/86400:5.2f}-{t1*tu/86400:5.2f} step {s:4.0f}s'
        for i in [0,1,2]:
            B=r[i,6]-r[i,7];C=r[i,7];H=r[i,3]+r[i,5];G=r[i,4]+r[i,11];R=r[i,9]+r[i,10]
            line+=f' | c{i} B{B:+7.1f} C{C:+7.1f} H{H:+7.1f} R{R:+5.1f} dU{r[i,2]-G:+7.1f}'
        print(line)
    print('cell0 per step, B erg and fraction of cell U (gas+rad, from dumps):')
    d0=P[0][1]
    for k,(t1,d1) in enumerate(P[1:],1):
        pass
    for (t0,d0_),(t1,d1) in zip(P[:-1],P[1:]):
        dt=(t1-t0)*tu; dd=round(t1*tu/86400*2); q=np.loadtxt(D[dd]); r=q[0,1]*Rs; dr=(q[1,1]-q[0,1])*Rs
        V=4*np.pi/3*((r+dr/2)**3-(r-dr/2)**3); U=(1.5*q[0,3]*Pu+q[0,7]*Eru)*V
        s=np.median(h[(h[:,0]>=t0)&(h[:,0]<t1),1])*tu; bs=((d1-d0_)[0,6]-(d1-d0_)[0,7])*F/dt*s
        print(f'  {t0*tu/86400:5.2f}-{t1*tu/86400:5.2f} B per step {bs:+.2e} erg  U_cell0(end)={U:.3e} erg  frac/step {bs/U:+.1e}  rate {bs/s/LS:+.1f} L*')
    d0=P[0][1]; t,d=P[-1]
    print('cum to',round(t*tu/86400,2),'d [1e45]: all B',round(((d[:,6]-d0[:,6]-d[:,7]+d0[:,7]).sum())*F/1e45,2),'cons',round((d[:,7]-d0[:,7]).sum()*F/1e45,2),'dE',round((d[:,2]-d0[:,2]).sum()*F/1e45,2))
    for dd in (5,):
        i=[k for k,p in enumerate(P) if abs(p[0]*tu/86400-dd)<.1][0]; t,d=P[i]; b=(d[:,6]-d0[:,6]-d[:,7]+d0[:,7])*F/1e45; c=(d[:,7]-d0[:,7])*F/1e45
        print(f'5 d: B all {b.sum():.2f} cells0-2 {b[:3].sum():.2f} cells3-10 {b[3:11].sum():.2f} cells>10 {b[11:].sum():.2f}; C all {c.sum():.2f}; dE {(d[:,2]-d0[:,2]).sum()*F/1e45:.2f}')
    # ledger.txt to the end (regions r<30, cells0-10, all)
    l=np.loadtxt(f'{a}/out/ledger.txt'); l=l[~np.isnan(l).any(axis=1)]
    n=1+11*0; 
    def col(reg,k): return 1+reg*12+1+k
    print('ledger.txt: t[d] all B, cells0-10 B, all C [1e45]:')
    for dd in [5,7,9,9.9]:
        j=np.argmin(abs(l[:,0]*tu/86400-dd))
        print(f'  {l[j,0]*tu/86400:.2f}: {(l[j,col(2,3)]-l[j,col(2,4)]-l[0,col(2,3)]+l[0,col(2,4)])*F/1e45:8.2f} {(l[j,col(1,3)]-l[j,col(1,4)]-l[0,col(1,3)]+l[0,col(1,4)])*F/1e45:8.2f} {(l[j,col(2,4)]-l[0,col(2,4)])*F/1e45:8.2f}')
