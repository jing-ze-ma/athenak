import numpy as np
P=np.load('prof.npz')
import glob; q=np.load(sorted(glob.glob('stars/m35mms/s_*.npz'))[-1]); rq=q['r']/9.554e11; R=np.sqrt(q['sp'][19:29].sum(0))/1e5; tot=q['rms']/1e5; 
for x in (0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.91,0.92,0.925,0.93):
    i=int(np.argmin(abs(rq-x))); N2=np.interp(rq[i],P['x'],P['N2']); print('r %.3f ring-band rms %.3e km/s, total rms %.3e, N2 %+.2e, band share %.2f, axis share 19-28 %.2f, K.E.'%(rq[i],R[i],tot[i],N2,(R[i]/tot[i])**2,q['sa'][19:29,i].sum()/q['sp'][19:29,i].sum()))
N2=P['N2']; x=P['x']; s=np.where((N2[1:]<0)&(N2[:-1]>0)&(x[1:]>0.8))[0]; print('first N2<0 above 0.80 at r',x[s[0]+1] if len(s) else None)
