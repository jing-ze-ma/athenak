import glob, numpy as np
RP=9.554e11; fs=sorted(glob.glob('out/q_*.npz')); D=[]
for f in fs:
    d=np.load(f); D.append(dict(t=float(d['t']),n=int(d['n']),r=d['r'][int(d['i0']):]/RP,C=d['C'],ky=d['ky'],kx=d['kx'],rms=d['rms'][int(d['i0']):],mean=d['mean'][int(d['i0']):]))
D.sort(key=lambda d:d['t']); t=np.array([d['t'] for d in D]); r=D[0]['r']; ky=D[0]['ky']; kx=D[0]['kx']; ir=lambda x:int(np.argmin(abs(r-x)))
KX,KY=np.meshgrid(kx,ky); kk=np.hypot(KX,KY); C=np.array([d['C'] for d in D]); P=abs(C)**2; i=ir(0.74)
bands={'1-2':(0.5,2.5),'3-8':(2.5,8.5),'9-18':(8.5,18.5),'19-27':(18.5,27.5),'28-30':(27.5,31)}
print('rms [km/s] at r=0.74 per |k| band')
for j in list(range(0,len(t),10)):
    print('n %d t %.3e'%(D[j]['n'],t[j]),' '.join('%s %.2e'%(b,np.sqrt(P[j][(kk>=lo)&(kk<hi)][:,i].sum())/1e5) for b,(lo,hi) in bands.items()))
ring=(kk>=18.5)&(kk<27.5); R=np.sqrt(P[:,ring,:].sum(1)); th=np.degrees(np.arctan2(abs(KY),abs(KX)))
print('switch at 1.8967e5 s; ring rms at 0.74 every 4 dumps from 1.2e5 to 3.0e5:'); m=(t>1.1e5)&(t<3.1e5); print(' '.join('%.3e:%.2e'%(a,b) for a,b in zip(t[m][::5],R[m,i][::5]/1e5)))
for nn in ('0:all',):
  for j in [int(np.argmin(abs(t-tt))) for tt in (3e5,4e5,5e5,6e5,7e5,7.7e5)]:
    Pj=P[j][:,:,i]*ring; tot=Pj.sum(); ax=Pj[(abs(KX)<=1)|(abs(KY)<=1)].sum(); iso=((abs(KX)<=1)|(abs(KY)<=1))[ring].sum()/ring.sum(); h=[Pj[(th>=a)&(th<a+15)].sum()/tot for a in range(0,90,15)]; hi=[(ring&(th>=a)&(th<a+15)).sum()/ring.sum() for a in range(0,90,15)]
    print('t %.2e ring power share within |kx|<=1 or |ky|<=1: %.2f (isotropic counting %.2f); by angle 0..90 in 15 deg bins (power / mode count):'%(t[j],ax/tot,iso),' '.join('%.2f/%.2f'%(a,b) for a,b in zip(h,hi)),'|k| median %.1f'%np.sort(kk[ring])[0])
# temporal spectra of strongest ring modes at 0.74 (t >= 5.8e5), detrended by the exponential growth
m=t>=5.8e5; tt=t[m]; dt=np.diff(tt).mean(); Pm=P[m][:,:,:,i].mean(0)*ring; q=np.argsort(Pm.ravel())[::-1][:8]
for qi in q:
    a,b=np.unravel_index(qi,Pm.shape); c=C[m][:,a,b,i]; g=np.polyfit(tt,np.log(abs(c)),1); z=c*np.exp(-np.polyval(g,tt)); F=np.fft.fft(z*np.hanning(len(z)),4096); f=np.fft.fftfreq(4096,dt); k=np.argmax(abs(F)**2); pk=np.sort(abs(F)**2)[::-1]
    f2=2*np.pi*f; print('mode (ky %d,kx %d): |c| %.4f km/s at end; growth e-fold %.2e s; peak omega %.2e rad/s (bin res %.1e, Nyquist %.2e, alias shift 2pi/dt=%.2e); share of top 3 resolution bins %.2f'%(ky[a],kx[b],abs(c[-1])/1e5,1/g[0],f2[k],2*np.pi/(len(z)*dt),np.pi/dt,2*np.pi/dt,(abs(F)**2)[max(0,k-2):k+3].sum()/ (abs(F)**2).sum()))
print('cadence',dt)
np.savez('ring.npz',t=t,r=r,R=R,mean=np.array([d['mean'] for d in D]))
