import numpy as np
X='/resnick/groups/carnegie_poc/jingze/tests_1124'
def ld(a):
    d=np.loadtxt(X+'/%s/m35mms.user.hst'%a,comments='#'); _,i=np.unique(d[:,0],return_index=True); d=d[np.sort(i)]; return d[:,0], d[:,5]/d[:,6], d[:,2]/d[:,6]
P=ld('P0')
for a in ('PR','PRL4','PL4'):
    A=ld(a); t=P[0][(P[0]>=A[0][0])&(P[0]<=A[0][-1])]
    for k,n in ((1,'L_top/L_in'),(2,'L_bot/L_in')):
        d=np.interp(t,A[0],A[k])-np.interp(t,P[0],P[k]); print('%-5s %-10s vs P0 on %d P0 times: rms %.3e max %.3e' % (a,n,t.size,np.sqrt(np.mean(d**2)),np.abs(d).max()))
