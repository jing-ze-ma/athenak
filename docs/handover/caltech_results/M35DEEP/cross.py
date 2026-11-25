import glob, numpy as np
BB=3
for run in ('m35mms','t9l62','t9mlt','m35zams','bsg'):
    S=[dict(np.load(f)) for f in sorted(glob.glob('stars3/%s/b_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); L=S[-6:]; r=S[0]['r']; pv=np.mean([z['pv'][:,BB] for z in L],0); av=np.mean([z['av'][:,BB] for z in L],0); frac=S[0]['cax'][BB]/S[0]['count'][BB]; x=(av/np.maximum(pv,1e-99))/frac; pT=np.mean([z['pT'][:,BB] for z in L],0); aT=np.mean([z['aT'][:,BB] for z in L],0); xT=(aT/np.maximum(pT,1e-99))/frac; rr=r/r[0]
    # coming from above: last radius (downward) where the ratio exceeds 2; mean in ranges
    o=[]; 
    for a,b in ((1.0,1.05),(1.05,1.10),(1.10,1.15),(1.15,1.20),(1.20,1.25),(1.25,1.30),(1.30,1.40),(1.40,1.60),(1.6,2.0),(2.0,3.0)):
        m=(rr>=a)&(rr<b)
        if m.any(): o.append('%.2f-%.2f: v %.1f T %.1f'%(a,b,np.median(x[m]),np.median(xT[m])))
    i=np.where(x>2)[0]; i=i[rr[i]<3]; print(run,'x = axis share / counting share, medians:',' | '.join(o)); 
    j=[k for k in range(len(rr)-1) if x[k]>2 and x[k+1]<=2 and rr[k]<2.5]; print('   crossings of x=2 (outward): r/r_in',[round(float(rr[k]),3) for k in j][:6])
