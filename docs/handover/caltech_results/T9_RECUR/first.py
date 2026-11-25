import sys; sys.argv=['x','/tmp/unused_t9']; import os, numpy as np
src=open('/resnick/groups/carnegie_poc/jingze/t9ceil_1124/t9hist.py').read().split('# ------------------------------------------------ per-run series')[0].replace("os.makedirs(OUT, exist_ok=True)","")
src=src.replace("open(OUT + '/runs.txt', 'w').write","(lambda *a: None)")
exec(src)
o=['run | first t with L_top/L_in > 1.1 / > 2.0 (t > 1e5 s) | first t dt < 48.4 (t > 1e5) | max L_top/L_in (t > 1e5) at t | max KE_int at t | last t  [runs starting at t < 1e5 and reaching 4e5]']
for r,(n,a) in INFO.items():
    if len(a)<20 or a[0,0]>1e5 or a[-1,0]<4e5: continue
    a=a[a[:,0]>1e5]; t=a[:,0]; q=a[:,n.index('L_top')]/a[:,n.index('L_in')]; K=a[:,n.index('KE_int')]; d=a[:,n.index('dt')]
    f=lambda m:'%.4e'%t[np.argmax(m)] if m.any() else 'none'
    o.append('%s | %s / %s | %s | %.4f at %.4e | %.4e at %.4e | %.4e'%(r,f(q>1.1),f(q>2),f(d<0.5*QUIET_DT),q.max(),t[q.argmax()],K.max(),t[K.argmax()],t[-1]))
L=['\nold-scheme / clean long run prod_fresh_1121b: which chunk wrote the surviving rows is not marked in the hst; chunks.txt is copied']
open('/resnick/groups/carnegie_poc/jingze/t9ceil_1124/hist/first.txt','w').write('\n'.join(o+L)+'\n'); print('\n'.join(o))
