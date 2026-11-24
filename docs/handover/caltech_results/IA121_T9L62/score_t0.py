import numpy as np, sys, re
Q='/resnick/groups/carnegie_poc/jingze/t9q_1123'
def load(a):
    d=np.loadtxt(Q+'/%s/t9l62.user.hst'%a, comments='#')
    t=d[:,0]; _,i=np.unique(t,return_index=True)  # drop repeated (dt_retry) rows, keep first
    d=d[np.sort(i)]; d=d[np.argsort(d[:,0])]
    return d[:,0], d[:,5]/d[:,6], d[:,2]/d[:,6]
A=load('T0a'); D=load('T0da'); H=load('T0h')
t0,t1=1.5776e6,1.607e6
tg=np.linspace(max(t0,A[0][0],D[0][0],H[0][0]), min(t1,A[0][-1],D[0][-1],H[0][-1]), 2001)
print('grid %.6e - %.6e, %d pts' % (tg[0], tg[-1], tg.size))
for k,name in ((1,'L_top/L_in'),(2,'L_bot/L_in')):
    a=np.interp(tg,A[0],A[k]); d=np.interp(tg,D[0],D[k]); h=np.interp(tg,H[0],H[k])
    for lab,x in (('T0da-T0a',d-a),('T0a-T0h',a-h)):
        print('%-11s %-9s rms %.3e max %.3e  (rel to mean %.4f)' % (name, lab, np.sqrt(np.mean(x**2)), np.max(np.abs(x)), np.mean(a)))
for a in ('T0a','T0da','T0h'):
    L=open(Q+'/%s/log.0'%a).read().splitlines(); el=[l for l in L if l.startswith('elapsed=') and 'cycle=' in l]
    e=[float(re.search(r'elapsed=(\S+)',l).group(1)) for l in el]; c=[int(re.search(r'cycle=(\d+)',l).group(1)) for l in el]; t=[float(re.search(r'time=(\S+)',l).group(1)) for l in el]
    print('%-5s s/cyc %.3f (c%d-%d)  wall per 1e4 s model %.0f s  (t %.4e - %.4e)' % (a,(e[-1]-e[0])/(c[-1]-c[0]),c[0],c[-1],(e[-1]-e[0])/(t[-1]-t[0])*1e4,t[0],t[-1]))
