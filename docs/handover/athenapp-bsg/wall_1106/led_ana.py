import numpy as np,glob,os
Rs=6.957e10;vu=np.sqrt(8.314462618e7*1e5/.62);tu=Rs/vu;Eu=1e-9*vu**2*Rs**3;LS=5.581e38;F=12*Eu
names='hydflux gravmass pointmass radT_bill radT_cons clip reset floorE floorM radD gasbill'.split()
W='/viper/ptmp2/jinma/athpp_ke_1106/wall_1106/'
for a in ['A_af1','A_tc4']:
    print('=====',a)
    # bitwise check vs the runs already sent
    nb=nd=0
    for f in sorted(glob.glob(f'{W}{a}/out/bsg.block0.out2.*.tab')):
        g=f.replace(W,W+'led/')
        if os.path.exists(g):
            nd+=1; nb+= open(f,'rb').read()==open(g,'rb').read()
    print('dumps byte-identical to wall_1106 run:',nb,'of',nd)
    P=[]
    for f in sorted(glob.glob(f'{a}/ledger_prof_*.txt'))+sorted(glob.glob(f'{a}/out/ledger_prof_*.txt')):
        t=float(os.path.basename(f)[12:-4]); d=np.loadtxt(f)
        if d.size and not np.isnan(d[:40]).any(): P.append((t,d))
    print('profiles t[d]:',[round(p[0]*tu/86400,2) for p in P])
    d0=P[0][1]
    # closure candidates at last valid profile
    t,d=P[-1]; dE=(d[:,2]-d0[:,2]); T=lambda k:d[:,3+k]-d0[:,3+k]
    c=T(0)+T(1)+T(2)+T(3)+T(6)+T(7)+T(8)
    print('closure dE vs hyd+grav+pm+radT_bill+reset+floors, cells 0,1,2,all [1e45]:',
      [f'{dE[i]*F/1e45:.3f}/{c[i]*F/1e45:.3f}' for i in range(3)],f'{dE.sum()*F/1e45:.3f}/{c.sum()*F/1e45:.3f}')
    print('t[d]  | cumulative billing loss B=radT_bill-radT_cons [1e45 erg full sphere] cells 0,1,2,3-10, all | dE cells 0,1,2,all')
    for t,d in P[1:]:
        B=(d[:,6]-d0[:,6]-(d[:,7]-d0[:,7]))*F/1e45; dE=(d[:,2]-d0[:,2])*F/1e45
        print(f'{t*tu/86400:5.2f} | {B[0]:8.3f} {B[1]:8.3f} {B[2]:8.3f} {B[3:11].sum():8.3f} {B.sum():8.3f} | {dE[0]:8.3f} {dE[1]:8.3f} {dE[2]:8.3f} {dE.sum():8.3f}')
    print('interval rates in L*=5.58e38: billing cell0,1,2,3-10,all | cons transport(boundary+net) all | hydro cell0 | radT_cons cell0 ')
    for (t0,d0_),(t1,d1) in zip(P[:-1],P[1:]):
        dt=(t1-t0)*tu; r=(d1-d0_)*F/dt/LS
        B=r[:,6]-r[:,7]
        print(f'{t0*tu/86400:5.2f}-{t1*tu/86400:5.2f} | {B[0]:8.1f} {B[1]:8.1f} {B[2]:8.1f} {B[3:11].sum():8.1f} {B.sum():8.1f} | {r[:,7].sum():8.1f} | {r[0,3]:8.1f} | {r[0,7]:8.1f}   per-step cell0 [erg] {B[0]*LS*dt/((t1-t0)/0.0103):.2e}')
    # Rough cell energy context
    t,d=P[3]; print('cell0 |E_meas| [erg]:',f'{abs(d[0,2])*F:.3e}')
