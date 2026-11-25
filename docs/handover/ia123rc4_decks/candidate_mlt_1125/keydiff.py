import re,sys
def par(f):
    b=None;o={}
    for l in open(f):
        l=l.split('#')[0].strip()
        if l.startswith('<') : b=l; continue
        if '=' in l and b: k,v=l.split('=',1); o[(b,k.strip())]=v.strip()
    return o
A,B=par(sys.argv[1]),par(sys.argv[2])
for k in sorted(set(A)|set(B)):
    if A.get(k)!=B.get(k): print('%s/%s : %s -> %s'%(k[0],k[1],A.get(k),B.get(k)))
