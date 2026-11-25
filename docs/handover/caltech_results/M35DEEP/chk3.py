import re, glob, numpy as np
RP=9.554e11; z=np.load('ring.npz'); tr=z['t']; r=z['r']; R=z['R']; mean=z['mean']; ir=lambda x:int(np.argmin(abs(r-x)))
nm=None; H=[]
for l in open('/resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35mms/m35mms.user.hst'):
    if l.startswith('#'):
        if '[1]=' in l: nm=re.findall(r'\]=(\S+)',l)
        continue
    p=l.split()
    if len(p)==len(nm): H.append([float(x) for x in p])
H=np.array(H); g={n:H[:,k] for k,n in enumerate(nm)}; t=g['time']; print('hst columns',nm[:20])
# PSD of hst series in 4.0e5-7.8e5 after removing a cubic trend (log for KE)
a,b=4.0e5,t[-1]; u=np.arange(a,b,12.0)
def line(y,lo,hi,label,logy=False):
    yi=np.interp(u,t,np.log(y) if logy else y); p=np.polyfit(u-u[0],yi,5); q=yi-np.polyval(p,u-u[0]); F=np.fft.rfft(q*np.hanning(len(q))); w=2*np.pi*np.fft.rfftfreq(len(q),12.0); P=abs(F)**2; m=(w>=5e-5)&(w<1e-3)
    k=np.argsort(P[m])[::-1][:6]; ks=np.sort(k); print('%-8s total rms %.3e; strongest omegas below 1e-3 rad/s:'%(label,q.std()),' '.join('%.2e(%.0f%%)'%(w[m][i],100*P[m][i]/P[m].sum()) for i in k)); return w,P
for n in ('L_top','L_midC','Mr_int','KEr_int','KE_int','PV_int','Etot','dt'):
    if n in g: line(g[n],5e-5,1e-3,n)
print('resolution of this window: %.2e rad/s'%(2*np.pi/(b-a)))
# ring rms modulation frequency
m=(tr>3.5e5)&(tr<=7.8e5); ln=np.log(R[:,ir(0.74)])[m]; tt=tr[m]; p=np.polyfit(tt,ln,4); q=ln-np.polyval(p,tt); dt=np.diff(tt).mean(); F=np.fft.rfft(q*np.hanning(len(q)),4096); w=2*np.pi*np.fft.rfftfreq(4096,dt); P=abs(F)**2; mm=(w>5e-5); k=np.argmax(P[mm]); print('ln ring rms(0.74) modulation: rms %.3f, peak omega %.3e rad/s (resolution %.1e, Nyquist %.1e)'%(q.std(),w[mm][k],2*np.pi/(tt[-1]-tt[0]),np.pi/dt))
print('   next:',' '.join('%.2e'%x for x in w[mm][np.argsort(P[mm])[::-1][:5]]))
# shell-mean v_r at 0.74 (radial pulsation) and its frequency
mv=mean[:,ir(0.74)][m]; print('shell mean v_r(0.74): rms %.3e cm/s'%mv.std()); F=np.fft.rfft((mv-np.polyval(np.polyfit(tt,mv,4),tt))*np.hanning(len(mv)),4096); P=abs(F)**2; print('   strongest omegas:',' '.join('%.2e(%.0f%%)'%(w[mm][x],100*P[mm][x]/P[mm].sum()) for x in np.argsort(P[mm])[::-1][:6]))
# growth rate series vs L and KE
g1=np.gradient(np.log(R[:,ir(0.74)]),tr); m2=(tr>4.0e5)&(tr<7.7e5); ts=tr[m2]; gs=np.convolve(g1,np.ones(7)/7,'same')[m2]
for n in ('L_top','L_midC','KE_int','KEr_int','L_intC','dt'):
    if n not in g: continue
    y=g[n]/(g['L_in'] if n.startswith('L') else 1.0); yi=np.array([y[(t>a0-900)&(t<=a0+900)].mean() for a0 in ts]); yd=yi-np.polyval(np.polyfit(ts,yi,4),ts); gd=gs-gs.mean(); r0=np.corrcoef(gd,yd)[0,1]
    lag=[np.corrcoef(gd[max(0,-s):len(gd)-max(0,s)],yd[max(0,s):len(yd)-max(0,-s)])[0,1] for s in range(-10,11)]; print('corr of d ln R/dt (7-dump mean) with %s (detrended): r=%.2f at lag 0; largest |r| over +-10 dumps %.2f at lag %d dumps; n %d, 95%% |r| for white noise ~%.2f'%(n,r0,max(lag,key=abs),int(np.argmax(np.abs(lag)))-10,len(gd),2/np.sqrt(len(gd)/7)))
# coherence of the ring-rms modulation with the hst KE_r line
v=np.interp(tt,t,g['KEr_int']); v=np.log(v)-np.polyval(np.polyfit(tt,np.log(v),4),tt); q2=q
F1=np.fft.rfft(q2*np.hanning(len(q2))); F2=np.fft.rfft(v*np.hanning(len(v))); ww=2*np.pi*np.fft.rfftfreq(len(q2),dt)
# window-averaged cross spectrum (4 overlapping segments)
def seg(x,n=60): 
    s=[]; 
    for k in range(0,len(x)-n+1,n//2): s.append(np.fft.rfft(x[k:k+n]*np.hanning(n)))
    return np.array(s)
A=seg(q2); B=seg(v); Cxy=(A*np.conj(B)).mean(0); coh=abs(Cxy)/np.sqrt((abs(A)**2).mean(0)*(abs(B)**2).mean(0)); wA=2*np.pi*np.fft.rfftfreq(60,dt)
for k in range(1,len(wA)):
    if 2e-4<wA[k]<5e-4: print('coherence ring rms vs ln KEr_int at omega %.2e: %.2f (segments %d; 95%% level for independent ~%.2f) phase %.2f rad'%(wA[k],coh[k],len(A),np.sqrt(1-0.05**(1/(len(A)-1))),np.angle(Cxy[k])))
np.savez('chk3.npz')
