import glob, numpy as np
BB=3; RB=dict(m35mms=(0.924*9.554e11,0.74*9.554e11),t9l62=(3.05e12,1.966e12),t9mlt=(2.95e12,1.966e12))
for run,(rb,rr) in RB.items():
    S=[dict(np.load(f)) for f in sorted(glob.glob('stars3/%s/b_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); r=S[0]['r']; t=np.array([float(z['t']) for z in S]); ib=int(np.argmin(abs(r-rb))); ii=int(np.argmin(abs(r-rr)))
    B=np.array([np.sqrt(z['pv'][ib,BB])/1e5 for z in S]); R=np.array([np.sqrt(z['pv'][ii,BB])/1e5 for z in S]); U=[]
    print('=====',run,'dumps',len(S),'base at r %.3e (r/r_in %.3f), ring at r %.3e'%(r[ib],r[ib]/r[0],r[ii]))
    for nm,x in (('base of convection zone',B),('ring',R)):
        o=[]
        for lv in (1e-5,1e-4,1e-3,1e-2,1e-1,1.0):
            k=np.where(x>lv)[0]; o.append('%g: %s'%(lv,'%.2e'%t[k[0]] if len(k) else '-'))
        print('  first time band 16-40 rms > level [km/s],',nm,'|',' ; '.join(o),'| last',x[-1])
    # intermediate radii: first time > 1e-3 km/s on the way between base and ring
    prog=[]
    for x in np.linspace(rr,rb,9):
        i=int(np.argmin(abs(r-x))); y=np.array([np.sqrt(z['pv'][i,BB])/1e5 for z in S]); k=np.where(y>1e-3)[0]; prog.append('r/r_in %.3f: %s'%(r[i]/r[0],'%.2e'%t[k[0]] if len(k) else 'never'))
    print('  first time >1e-3 km/s from ring to base:',' | '.join(prog))
    m=(t>0); lB=np.log(B+1e-12); lR=np.log(R+1e-12)
    g=lambda y: np.gradient(y,t)
    for nm,x in (('d ln base/dt',g(lB)),('d ln ring/dt',g(lR))): print('  ',nm,'range %.2e..%.2e 1/s'%(np.percentile(x[2:],5),np.percentile(x[2:],95)))
    # lead/lag of growth: correlation of detrended ln amplitudes after both exceed 1e-4
    ok=np.where((B>1e-4)&(R>1e-4))[0]; s=ok[0] if len(ok) else 0
    a=lB[s:]-np.polyval(np.polyfit(t[s:],lB[s:],3),t[s:]); b=lR[s:]-np.polyval(np.polyfit(t[s:],lR[s:],3),t[s:]); L=len(a); c=[]
    for lag in range(-8,9):
        x=a[max(0,-lag):L-max(0,lag)]; y=b[max(0,lag):L-max(0,-lag)]; c.append(np.corrcoef(x,y)[0,1])
    print('   corr of detrended ln base vs ln ring (dumps, dt %.2e s, n %d): lag0 %.2f; max |r| %.2f at lag %d dumps (positive = base later); white-noise 95%% ~%.2f'%(np.median(np.diff(t)),L,c[8],max(c,key=abs),int(np.argmax(np.abs(c)))-8,2/np.sqrt(L)))
