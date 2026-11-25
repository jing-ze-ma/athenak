import glob, numpy as np, sys
RP=dict(m35mms=9.554e11,t9l62=6.5326e12,t9mlt=3.9282e12,m35zams=5.4972e11,bsg=3.4e12)
RING=dict(m35mms=0.74*9.554e11,t9l62=1.966e12,t9mlt=1.966e12)
BB=3  # band index (16,40)
for run in ('m35mms','t9l62','t9mlt','m35zams','bsg'):
    S=[dict(np.load(f)) for f in sorted(glob.glob('stars3/%s/b_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); r=S[0]['r']; t=np.array([float(z['t']) for z in S]); ni=len(r); cnt=S[0]['count'][BB]; cax=S[0]['cax'][BB]; N3,N2=S[0]['N']; frac=cax/cnt
    L=S[-6:]; pv=np.mean([z['pv'][:,BB] for z in L],0); av=np.mean([z['av'][:,BB] for z in L],0); pT=np.mean([z['pT'][:,BB] for z in L],0); aT=np.mean([z['aT'][:,BB] for z in L],0); tv=np.mean([z['rmsv'] for z in L],0); tT=np.mean([z['rmsT'] for z in L],0)
    P=np.load('prof.npz') if run=='m35mms' else np.load('stars/prof_%s.npz'%run); 
    print('=====',run,'cells',ni,'dumps',len(S),'t last %.3e'%t[-1],'count share of band 16-40 (k256; 0.0625-0.156 cycles/cell = 6.4-16 cells): %.3f'%frac,'N3,N2',N3,N2)
    n2=lambda rr: np.interp(rr,P['r'] if 'r' in P else P['x']*RP[run],P['N2'],left=np.nan,right=np.nan)
    for i in [int(x) for x in np.linspace(1,ni-1,22)]:
        ax=av[i]/max(pv[i],1e-99); aTt=aT[i]/max(pT[i],1e-99); print('  r %.3e (r/r_in %.3f, r/RPH %.3f) N2 %s | band rms v %.2e km/s (%.0f%% of rms^2), rms T %.2e (%.0f%%); axis share v %.2f (x%.1f), T %.2f (x%.1f)'%(r[i],r[i]/r[0],r[i]/RP[run],'%+.1e'%n2(r[i]) if np.isfinite(n2(r[i])) else 'n/a',np.sqrt(pv[i])/1e5,100*pv[i]/tv[i]**2,np.sqrt(pT[i]),100*pT[i]/max(tT[i]**2,1e-99),ax,ax/frac,aTt,aTt/frac))
    share=av/np.maximum(pv,1e-99)/frac; stt=np.where(share>2)[0]; print('   radii with axis share > 2x counting:', 'none' if len(stt)==0 else 'r/r_in %.3f-%.3f (n %d of %d); contiguous groups:'%(r[stt[0]]/r[0],r[stt[-1]]/r[0],len(stt),ni), [ (round(r[a]/r[0],3),round(r[b]/r[0],3)) for a,b in zip([stt[0]]+[stt[k+1] for k in range(len(stt)-1) if stt[k+1]!=stt[k]+1],[stt[k] for k in range(len(stt)-1) if stt[k+1]!=stt[k]+1]+[stt[-1]])] if len(stt) else '')
