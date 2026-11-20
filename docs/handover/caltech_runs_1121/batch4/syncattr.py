#!/usr/bin/env python3
# syncattr.py SQLITE [C0 NC PERCYC]: rank-0 host-sync calls and kernel launches per cycle, attributed to the
# innermost enclosing NVTX (Kokkos label) range on the calling thread (window as r0lab.py: cycles C0+1..C0+NC)
import sys, sqlite3, bisect, collections
db=sqlite3.connect(sys.argv[1]); c=db.cursor()
c0=int(sys.argv[2]) if len(sys.argv)>2 else 1; nc=int(sys.argv[3]) if len(sys.argv)>3 else 3; per=float(sys.argv[4]) if len(sys.argv)>4 else 2
S=dict(c.execute("select id,value from StringIds"))
kn=c.execute("select k.start,k.demangledName from CUPTI_ACTIVITY_KIND_KERNEL k").fetchall()
mk=sorted(st for st,dn in kn if 'NewTimeStep' in S.get(dn,''))
cl=[]
for t in mk:
    if not cl or t-cl[-1][-1]>2e8: cl.append([t])
    else: cl[-1].append(t)
ends=[x[0] for x in cl]; w0,w1=ends[c0],ends[min(c0+nc,len(ends)-1)]; ncy=(min(c0+nc,len(ends)-1)-c0)/per
nv=c.execute("select start,end,coalesce(text,''),textId,globalTid from NVTX_EVENTS where end is not null and eventType=59 and end>=? and start<=?",(w0,w1)).fetchall()
rng=collections.defaultdict(list)
for st,en,tx,tid,g in nv: rng[g].append((st,en,(tx or S.get(tid,'')).lstrip(':')))
for g in rng: rng[g].sort()
starts={g:[r[0] for r in rng[g]] for g in rng}
def lab(st,g):
    L=rng.get(g,[]); a=starts.get(g,[]); i=bisect.bisect_right(a,st)-1; best=None
    while i>=0 and i>bisect.bisect_right(a,st)-600:
        s0,e0,n=L[i]
        if e0>=st and (best is None or s0>=best[0]): best=(s0,n); break
        i-=1
    return best[1] if best else '?'
rt=c.execute("select r.start,r.end,r.nameId,r.globalTid from CUPTI_ACTIVITY_KIND_RUNTIME r where r.start>=? and r.start<=?",(w0,w1)).fetchall()
sy=collections.Counter(); syt=collections.Counter(); la=collections.Counter(); kinds=collections.Counter()
for st,en,nid,g in rt:
    n=S.get(nid,'')
    if 'Synchronize' in n or n.startswith('cudaMemcpy') and 'Async' not in n:
        l=lab(st,g); sy[l]+=1; syt[l]+=(en-st); kinds[n.split('_')[0]]+=1
    elif 'LaunchKernel' in n: la[lab(st,g)]+=1
print('window cycles %.1f  syncs/cyc %.0f (%.0f ms)  launches/cyc %.0f'%(ncy,sum(sy.values())/ncy,sum(syt.values())/ncy/1e6,sum(la.values())/ncy))
print('kinds', {k:round(v/ncy) for k,v in kinds.most_common(8)})
print('%-60s %8s %8s %8s'%('label','sync/cy','ms/cy','launch/cy'))
for l,v in sorted(sy.items(),key=lambda x:-x[1])[:int(sys.argv[5]) if len(sys.argv)>5 else 60]:
    print('%-60s %8.0f %8.1f %8.0f'%(l[:60],v/ncy,syt[l]/ncy/1e6,la[l]/ncy))
print('--- top launches'); 
for l,v in la.most_common(40): print('%-60s %8.0f %8.0f'%(l[:60],v/ncy,sy[l]/ncy))
