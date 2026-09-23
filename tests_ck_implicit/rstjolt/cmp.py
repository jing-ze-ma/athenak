import sys,numpy as np
def rdhst(f):
    L=[l for l in open(f) if not l.startswith('#')]
    return np.array([[float(x) for x in l.split()] for l in L])
a=rdhst(sys.argv[1]); b=rdhst(sys.argv[2])
# align by time (col 0)
tb={r[0]:r for r in b}
first=None
for r in a:
    if r[0] in tb:
        d=np.abs(r-tb[r[0]]); rel=d/np.maximum(np.abs(r),1e-300)
        if first is None and d.max()>0: first=r[0]
        print('%.10e maxrel %.3e col %d'%(r[0],rel.max(),rel.argmax()))
print('first differing time',first)
