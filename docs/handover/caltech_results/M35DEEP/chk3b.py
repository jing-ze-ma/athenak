import re, numpy as np
z=np.load('ring.npz'); tr=z['t']; r=z['r']; R=z['R']; mean=z['mean']; i=int(np.argmin(abs(r-0.74)))
nm=None; H=[]
for l in open('/resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35mms/m35mms.user.hst'):
    if l.startswith('#'):
        if '[1]=' in l: nm=re.findall(r'\]=(\S+)',l)
        continue
    p=l.split()
    if len(p)==len(nm): H.append([float(x) for x in p])
H=np.array(H); g={n:H[:,k] for k,n in enumerate(nm)}; t=g['time']
m=(tr>3.6e5)&(tr<=7.78e5); tt=tr[m]; ln=np.log(R[m,i]); q=ln-np.polyval(np.polyfit(tt,ln,5),tt)
S={'ring ln rms':q,'L_midC':None,'Mr_int':None,'L_top':None,'mean v_r(0.74)':None,'PV_int':None}
for n in ('L_midC','Mr_int','L_top','PV_int'):
    y=np.array([g[n][(t>a-900)&(t<=a+900)].mean() for a in tt]); S[n]=y-np.polyval(np.polyfit(tt,y,5),tt); S[n]=S[n]/S[n].std()
mv=mean[m,i]; S['mean v_r(0.74)']=(mv-np.polyval(np.polyfit(tt,mv,5),tt))/mv.std(); S['ring ln rms']=q/q.std()
def demod(x,w0,sig):
    out=[]
    for a in tt: 
        k=np.exp(-0.5*((tt-a)/sig)**2); out.append((x*np.exp(-1j*w0*tt)*k).sum()/k.sum())
    return np.array(out)
def peak(x):
    f=np.linspace(2.2e-4,5e-4,561); return f,[abs((x*np.hanning(len(x))*np.exp(-1j*w*tt)).sum()) for w in f]
for n,x in S.items():
    f,p=peak(x); p=np.array(p); o=np.argsort(p)[::-1][:1]; print('%-15s peak in 2.2e-4-5e-4: %.3e rad/s (resolution ~%.1e); power near 3.2e-4: %.2f, near 3.36e-4: %.2f, near 4.0e-4: %.2f (relative to peak)'%(n,f[o][0],2*np.pi/(tt[-1]-tt[0]),p[np.argmin(abs(f-3.2e-4))]/p.max(),p[np.argmin(abs(f-3.36e-4))]/p.max(),p[np.argmin(abs(f-4.0e-4))]/p.max()))
w0=3.36e-4; sig=4e4
A={n:demod(x,w0,sig) for n,x in S.items()}
sel=np.arange(0,len(tt),44)  # about one independent sample per 8e4 s
for n in ('L_midC','L_top','Mr_int','PV_int','mean v_r(0.74)'):
    d=np.angle(A['ring ln rms']*np.conj(A[n])); R_=abs(np.exp(1j*d[sel]).mean()); print('phase difference ring rms - %s at 3.36e-4 (window 4e4 s): samples %s rad; resultant length %.2f (%d samples, random ~%.2f)'%(n,' '.join('%.1f'%x for x in d[sel]),R_,len(sel),np.sqrt(-np.log(0.05)/len(sel))))
