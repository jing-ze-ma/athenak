import numpy as np, sys, os
VMLT=1.86e4; CS=2.1e6
def load(d):
    f=os.path.join(d,'m1col.user.hst')
    a=np.loadtxt(f)
    return a
def analyse(d,t0=None,t1=None,col=8):
    a=load(d); t=a[:,0]
    if len(t)<50: return None
    v=a[:,col]          # 8 = V1mid (signed), 9 = KEcol
    ke=a[:,9]
    # window
    lo = t0 if t0 is not None else 100.0
    hi = t1 if t1 is not None else t[-1]
    m=(t>=lo)&(t<=hi)
    tw=t[m]; vw=v[m]; kw=ke[m]
    # ---- period from zero crossings of the detrended signal
    # detrend with a running mean over ~1 window
    vv=vw-np.mean(vw)
    sg=np.sign(vv); idx=np.where(np.diff(sg)!=0)[0]
    P=np.nan; Perr=np.nan
    if len(idx)>4:
        # linear interp of crossing times
        tc=tw[idx]+(tw[idx+1]-tw[idx])*(-vv[idx])/(vv[idx+1]-vv[idx])
        halfP=np.diff(tc)
        P=2*np.mean(halfP); Perr=2*np.std(halfP)/np.sqrt(len(halfP))
    # ---- period from FFT of the detrended signal
    dt=np.median(np.diff(tw))
    n=len(vv); w=np.hanning(n)
    ft=np.abs(np.fft.rfft(vv*w)); fr=np.fft.rfftfreq(n,dt)
    k=np.argmax(ft[1:])+1
    Pf=1.0/fr[k] if fr[k]>0 else np.nan
    # ---- growth rate from ln KE / 2  (KE ~ exp(2 gamma t))
    g=np.nan; gerr=np.nan
    good=kw>0
    if good.sum()>20:
        y=0.5*np.log(kw[good]); x=tw[good]
        A=np.vstack([x,np.ones_like(x)]).T
        c,res,_,_=np.linalg.lstsq(A,y,rcond=None)
        g=c[0]
        yfit=A@c; s=np.std(y-yfit)
        gerr=s/np.sqrt(np.sum((x-x.mean())**2))
    return dict(P=P,Perr=Perr,Pf=Pf,g=g,gerr=gerr,gP=g*P,t0=lo,t1=hi,
                vend=a[-1,5],vend_cs=a[-1,5]/CS,tend=t[-1])
if __name__=='__main__':
    for d in sys.argv[1:]:
        # allow dir:t0:t1
        p=d.split(':'); dd=p[0]
        t0=float(p[1]) if len(p)>1 else None
        t1=float(p[2]) if len(p)>2 else None
        r=analyse(dd,t0,t1)
        if r is None: print('%-16s too short'%dd); continue
        print('%-16s t=%.0f..%.0f  P_zc=%6.2f+-%.2f  P_fft=%6.2f  gamma=%9.3e+-%.1e  1/gamma=%7.1f  gamma*P=%7.4f  |v1|max_end=%.2e (%.3f cs)'%(
            dd,r['t0'],r['t1'],r['P'],r['Perr'],r['Pf'],r['g'],r['gerr'],1/r['g'] if r['g'] else 0,r['gP'],r['vend'],r['vend_cs']))
