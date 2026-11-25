import glob, numpy as np
RP=9.554e11; fs=sorted(glob.glob('out/q_*.npz')); D=[]
for f in fs:
    d=np.load(f); D.append(dict(t=float(d['t']),n=int(d['n']),r=d['r'][int(d['i0']):]/RP,C=d['C'],ky=d['ky'],kx=d['kx'],rms=d['rms'][int(d['i0']):],mean=d['mean'][int(d['i0']):],vt=d['vtrms'][int(d['i0']):],drr=d['drr'][int(d['i0']):]))
D.sort(key=lambda d:d['t']); t=np.array([d['t'] for d in D]); r=D[0]['r']; ky=D[0]['ky']; kx=D[0]['kx']; ir=lambda x:int(np.argmin(abs(r-x)))
KX,KY=np.meshgrid(kx,ky); kk=np.hypot(KX,KY); ring=(kk>=19)&(kk<=27); low=(kk>=0.5)&(kk<=6)
C=np.array([d['C'] for d in D])    # [t,ky,kx,r]
P=abs(C)**2; R=np.sqrt(P[:,ring,:].sum(1)); Lw=np.sqrt(P[:,low,:].sum(1)); tot=np.sqrt(P[:, :, :,:].reshape(len(t),-1,len(r)).sum(1))
print('dumps',len(t),t[0],t[-1],'spacing',np.diff(t)[-1])
for x in (0.72,0.73,0.74,0.75,0.76,0.77,0.80):
    print('r',x,'ring (19<=|k|<=27) rms [km/s] at t=',' '.join('%.1e:%.4f'%(t[j],R[j,ir(x)]/1e5) for j in [int(np.argmin(abs(t-tt))) for tt in (1.8e5,2.4e5,3.0e5,3.6e5,4.2e5,4.8e5,5.4e5,6e5,6.6e5,7.2e5,7.7e5)]),'| low (1..6):',' '.join('%.4f'%(Lw[j,ir(x)]/1e5) for j in [int(np.argmin(abs(t-tt))) for tt in (1.8e5,3.6e5,5.4e5,7.2e5)]))
i=ir(0.74); lnR=np.log(R[:,i])
for a,b in ((1.6e5,2.6e5),(2.4e5,3.2e5),(3.0e5,4.0e5),(4e5,5e5),(5e5,6e5),(6e5,7e5),(6.8e5,7.8e5),(3.0e5,7.8e5),(3.6e5,7.8e5)):
    m=(t>=a)&(t<=b); p=np.polyfit(t[m],lnR[m],1); print('ring rms at 0.74, t %.1e-%.1e n %d: e-fold %.3e s, ln residual rms %.3f'%(a,b,m.sum(),1/p[0],(lnR[m]-np.polyval(p,t[m])).std()))
print('ring growth per radius, 3.6e5-7.0e5 s and 7.0e5-7.8e5:')
for x in (0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.80):
    j=ir(x); s='%.2f '%x
    for a,b in ((3.6e5,7.0e5),(7.0e5,7.8e5)):
        m=(t>=a)&(t<=b); p=np.polyfit(t[m],np.log(R[m,j]),1); s+=' e-fold %.3e s'%(1/p[0])
    print(s)
# dominant axis modes
sel=[(22,0),(0,21),(22,1),(21,0),(23,0)]
for (a,b) in sel:
    ia=np.where(ky==a)[0][0]; ib=np.where(kx==b)[0][0]; c=C[:,ia,ib,i]; 
    m=(t>=6.0e5); ph=np.unwrap(np.angle(c[m])); am=np.log(abs(c[m])); p=np.polyfit(t[m],ph,1); dph=np.diff(ph); print('mode (ky %d,kx %d) at 0.74, t>=6.0e5: amplitude %.4f->%.4f km/s; phase step per 1800 s: mean %+.3f rad std %.3f; linear slope %.3e rad/s (mod 2pi/1800 = 3.49e-3), growth e-fold %.3e s'%(a,b,abs(c[m])[0]/1e5,abs(c[m])[-1]/1e5,dph.mean(),dph.std(),p[0],1/np.polyfit(t[m],am,1)[0]))
# phase vs r for the mode (22,0) at last dumps, and its variation over the layer
for (a,b) in ((22,0),(0,21)):
    ia=np.where(ky==a)[0][0]; ib=np.where(kx==b)[0][0]
    for j in (-1,-5,-20):
        c=C[j,ia,ib,:]; ph=np.unwrap(np.angle(c)); print('phase vs r mode (%d,%d) t %.3e:'%(a,b,t[j]),' '.join('%.3f:%.4f/%+.2f'%(r[q],abs(c[q])/1e5,ph[q]) for q in range(ir(0.708),ir(0.80),4)))
np.savez('ring.npz',t=t,r=r,R=R,low=Lw)
