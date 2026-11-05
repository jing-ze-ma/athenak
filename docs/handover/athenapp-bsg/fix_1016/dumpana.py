import numpy as np, sys
a=np.loadtxt(sys.argv[1]); nang=int(sys.argv[2]) if len(sys.argv)>2 else 8
s=a[:,0].astype(int)
I=a[:,2:2+2*nang:2]; d=a[:,3:3+2*nang:2]; T=a[:,-1]
for ss in sorted(set(s)):
    if ss not in [1,2,5,10,50,100,200,500,1000,1500,2000,2500,3000]: continue
    m=s==ss
    rel=np.abs(d[m]).sum(1)/I[m].sum(1)
    g=np.abs(d[m]).sum()/I[m].sum()
    top=np.argsort(rel)[::-1][:6]
    print(ss,'global %.2e'%g,'top cells',top,'log rel',np.round(np.log10(rel[top]+1e-300),1), 'frac numer cells0-10 %.2f'%(np.abs(d[m][:11]).sum()/np.abs(d[m]).sum()))
for ss in [100,1000]:
  m=s==ss
  if not m.any(): continue
  rel=np.abs(d[m]).sum(1)/I[m].sum(1)
  print(ss,'log rel every 8 cells:',np.round(np.log10(rel+1e-300),1)[::8])
print('log I sum per cell every 8:',np.round(np.log10(I[m].sum(1)),1)[::8])
print('T every 8:',np.round(T[m],3)[::8])
