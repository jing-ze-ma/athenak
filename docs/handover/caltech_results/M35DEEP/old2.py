import glob,numpy as np
RP=9.554e11; ref=[dict(np.load(f)) for f in sorted(glob.glob('stars/m35mms/s_*.npz'))]; ref.sort(key=lambda d: float(d['t']))
q=np.load('ring.npz'); ir=int(np.argmin(abs(q['r']-0.74))); rt=q['t']; rR=q['R'][:,ir]/1e5
for run in ('i1','i1f1k','fD0','fD1','fD2','fE0','fE4','fD1h'):
    fs=sorted(glob.glob('stars/%s/s_*.npz'%run))
    if not fs: print(run,'no output'); continue
    S=[dict(np.load(f)) for f in fs]; S.sort(key=lambda d: float(d['t'])); r=S[0]['r']/RP; i=int(np.argmin(abs(r-0.74))); out=[]
    for z in S[::max(1,len(S)//7)]+[S[-1]]:
        t=float(z['t']); R=np.sqrt(z['sp'][19:29,i].sum())/1e5; c=np.exp(np.interp(t,rt[1:],np.log(rR[1:]+1e-30))) if t>2e4 else np.nan; sp=z['sp'][:,i]; out.append('%.3e: %.2e (c1 %.2e) ax %.2f pk %d'%(t,R,c,z['sa'][16:41,i].sum()/max(sp[16:41].sum(),1e-99),12+np.argmax(sp[12:60])))
    print(run,len(S),'|',' ; '.join(out))
