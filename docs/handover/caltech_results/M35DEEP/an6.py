import re, numpy as np
nm=None; H=[]
for l in open('/resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35mms/m35mms.user.hst'):
    if l.startswith('#'):
        if '[1]=' in l: nm=re.findall(r'\]=(\S+)',l)
        continue
    p=l.split()
    if len(p)==len(nm): H.append([float(x) for x in p])
H=np.array(H); g={n:H[:,k] for k,n in enumerate(nm)}; t=g['time']; print('hst to',t[-1])
for a,b in ((6.5e5,7.0e5),(7.2e5,7.6e5),(7.5e5,t[-1])):
    m=(t>a)&(t<b); y=g['KEr_int'][m]; tt=t[m]; u=np.arange(tt[0],tt[-1],12.0); yi=np.interp(u,tt,y); p=np.polyfit(u-u[0],np.log(yi),3); z=np.exp(np.log(yi)-np.polyval(p,u-u[0]))-1
    F=np.fft.rfft(z*np.hanning(len(z))); w=2*np.pi*np.fft.rfftfreq(len(z),12.0); P=abs(F)**2
    def band(lo,hi): return np.sqrt(P[(w>=lo)&(w<hi)].sum()/P[w>1e-5].sum())
    print('KEr_int, t %.2e-%.2e: relative rms after cubic-log detrend %.2e; share of rms in omega bands: 1e-4-2.5e-4 %.2f, 2.5e-4-4e-4 (2N of 0.74) %.2f, 4e-4-1.5e-3 %.2f, 1.5e-3-9e-3 %.2f, 9e-3-1.2e-2 (cs kh) %.3f, 1.9e-2-2.3e-2 (2 cs kh) %.4f, >1.2e-2 %.3f; resolution %.1e'%(a,b,z.std(),*[band(*q) for q in ((1e-4,2.5e-4),(2.5e-4,4e-4),(4e-4,1.5e-3),(1.5e-3,9e-3),(9e-3,1.2e-2),(1.9e-2,2.3e-2),(1.2e-2,1))],w[1]))
    k=np.argsort(P[(w>1e-5)])[::-1][:6]; print('   strongest omegas',w[w>1e-5][k])
