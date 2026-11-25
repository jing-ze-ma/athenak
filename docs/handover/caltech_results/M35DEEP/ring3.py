import glob, numpy as np, sys
REF=dict(m35mms=[(0.72,9.554e11),(0.74,9.554e11),(0.76,9.554e11)],t9l62=[(1.75e12,1),(1.85e12,1),(1.97e12,1),(2.08e12,1)],t9mlt=[(1.75e12,1),(1.85e12,1),(1.97e12,1),(2.08e12,1)],m35zams=[(4.85e11,1),(4.95e11,1),(5.05e11,1),(5.15e11,1)],bsg=[(1.45e12,1),(1.8e12,1),(2.1e12,1),(2.4e12,1)])
for run in REF:
    S=[dict(np.load(f)) for f in sorted(glob.glob('stars/%s/s_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); r=S[0]['r']; t=np.array([float(z['t']) for z in S]); N3,N2=S[0]['N']; ex=S[0]['ext']; thc=(ex[1]-ex[0])/N2; P=np.load('stars/prof_%s.npz'%run)
    print('=====',run)
    for x,u in REF[run]:
        i=int(np.argmin(abs(r-x*u))); z=S[-1]; sp=z['sp'][:,i]; k=np.arange(182); N=np.sqrt(max(np.interp(r[i],P['r'],P['N2']),0)); kb=np.argmax(sp[12:60])+12
        tot=np.sqrt(sp[1:].sum()); rg=np.sqrt(sp[16:41].sum()); ax=z['sa'][16:41,i].sum()/sp[16:41].sum()
        ms=[np.sqrt(s['sp'][16:41,i].sum())/1e5 for s in S]; chi=np.interp(r[i],P['r'],P['chi']); Hp=np.interp(r[i],P['r'],P['Hp']); kh=2*np.pi*kb/(256*thc*r[i])
        print(' r %.3e: N %.2e, total rms %.2e km/s, rms in k256 16-40 %.2e km/s (%.1f%% of the total rms^2), peak bin in 12-60: %d (%.1f cells), axis share %.2f (count %.2f) | rms(16-40) vs t:'%(r[i],N,tot/1e5,rg/1e5,100*rg**2/tot**2,kb,256/kb,ax,z['cax'][16:41].sum()/z['cnt'][16:41].sum()),' '.join('%.1e:%.1e'%(t[j],ms[j]) for j in range(0,len(t),max(1,len(t)//9))),'| l %.0f kh*Hp %.1f 1/(chi kh2) %.1e s N*dt*? cells/dr %.2f'%(2*np.pi*kb/(256*thc),kh*Hp,1/(chi*kh**2),thc*r[i]/np.gradient(r)[i]))
