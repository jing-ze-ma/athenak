#!/usr/bin/env python3
# r0lab.py SQLITE [C0] [NC] [NTOP]: rank-0 NVTX (Kokkos label) attribution of GPU kernel+memcpy time over cycles C0+1..C0+NC
# (cycle ends = NVTX ranges whose label contains 'NewTimeStep' or kernels named Hydro::NewTimeStep); prints JSON-ish table
import sys, sqlite3, bisect, collections, json, re
db=sqlite3.connect(sys.argv[1]); c=db.cursor()
c0=int(sys.argv[2]) if len(sys.argv)>2 else 1; nc=int(sys.argv[3]) if len(sys.argv)>3 else 3; NT=int(sys.argv[4]) if len(sys.argv)>4 else 80
tabs={r[0] for r in c.execute("select name from sqlite_master where type='table'")}
S=dict(c.execute("select id,value from StringIds"))
nv=c.execute("select start,end,coalesce(text,''),textId,globalTid from NVTX_EVENTS where end is not null and eventType=59").fetchall()
rng=collections.defaultdict(list)
for st,en,tx,tid,g in nv: rng[g].append((st,en,(tx or S.get(tid,'')).lstrip(':')))
for g in rng: rng[g].sort()
starts={g:[r[0] for r in rng[g]] for g in rng}
# cycle markers: kernel names containing NewTimeStep
kn=c.execute("select k.start,k.end,k.correlationId,k.shortName,k.demangledName from CUPTI_ACTIVITY_KIND_KERNEL k").fetchall()
mk=sorted(st for st,en,cid,sn,dn in kn if 'NewTimeStep' in S.get(dn,'') )
# cluster markers (several per cycle) -> first of clusters separated by >0.2 s
cl=[]
for t in mk:
    if not cl or t-cl[-1][-1]>2e8: cl.append([t])
    else: cl[-1].append(t)
ends=[x[0] for x in cl]
if len(ends)<c0+nc+1: print('too few cycle markers',len(ends)); 
w0,w1=ends[c0],ends[min(c0+nc,len(ends)-1)]; ncy=(min(c0+nc,len(ends)-1)-c0)/float(sys.argv[5] if len(sys.argv)>5 else 2)
rt={cid:(st,g) for cid,st,g in c.execute("select correlationId,start,globalTid from CUPTI_ACTIVITY_KIND_RUNTIME")}
def label(cid):
    if cid not in rt: return '?'
    st,g=rt[cid]; L=rng.get(g,[]); a=starts.get(g,[]); i=bisect.bisect_right(a,st)-1; j=i
    while i>=0 and i>j-400:
        s0,e0,n=L[i]
        if s0<=st<=e0: return n
        i-=1
    return '?'
tot=collections.Counter(); cnt=collections.Counter(); busy=[]
for st,en,cid,sn,dn in kn:
    if en<w0 or st>w1: continue
    n=label(cid); tot[n]+=en-st; cnt[n]+=1; busy.append((st,en))
for t in ('CUPTI_ACTIVITY_KIND_MEMCPY','CUPTI_ACTIVITY_KIND_MEMSET'):
    if t in tabs:
        for st,en,cid in c.execute(f"select start,end,correlationId from {t}"):
            if en<w0 or st>w1: continue
            n=label(cid)+'['+t[20:23]+']'; tot[n]+=en-st; cnt[n]+=1; busy.append((st,en))
busy.sort(); u=0; cs=ce=None
for s,e in busy:
    if cs is None or s>ce:
        if cs is not None: u+=ce-cs
        cs,ce=s,e
    else: ce=max(ce,e)
if cs is not None: u+=ce-cs
# host NVTX range wall (non-kernel ranges i.e. labels that are region markers), instance counts per cycle of all labels
icnt=collections.Counter(); iwall=collections.Counter()
for g in rng:
    for s,e,n in rng[g]:
        if w0<=s<=w1: icnt[n]+=1; iwall[n]+=e-s
mpi=0; mpic=collections.Counter()
for t in tabs:
    if t.startswith('MPI_'):
        try:
            for row in c.execute(f"select start,end from {t}"):
                s,e=row
                if s>=w0 and e<=w1: mpi+=e-s
        except Exception: pass
# runtime sync time
sy=0; nsy=0
for st,en,nid in c.execute("select r.start,r.end,r.nameId from CUPTI_ACTIVITY_KIND_RUNTIME r"):
    if st>=w0 and en<=w1 and 'ynchronize' in S.get(nid,''): sy+=en-st; nsy+=1
wall=(w1-w0)/ncy/1e6
print(f"cycles {ncy} wall {wall:.1f} ms/cyc busy {u/ncy/1e6:.1f} idle {wall-u/ncy/1e6:.1f} MPI(host) {mpi/ncy/1e6:.1f} sync {sy/ncy/1e6:.1f} ms in {nsy/ncy:.0f} calls; launches {sum(cnt.values())/ncy:.0f}; cycle-ends(ms) {[round((ends[i+1]-ends[i])/1e6) for i in range(len(ends)-1)]}")
print("label  gpu_ms/cyc  n/cyc  nvtx_wall_ms/cyc")
for n,v in tot.most_common(NT): print(f"  {n[:60]:60s} {v/ncy/1e6:8.2f} {cnt[n]/ncy:7.1f} {iwall.get(n,0)/ncy/1e6:8.2f}")
json.dump({'wall':wall,'busy':u/ncy/1e6,'gpu':{k:v/ncy/1e6 for k,v in tot.items()},'n':{k:v/ncy for k,v in cnt.items()},'inst':{k:v/ncy for k,v in icnt.items()},'iwall':{k:v/ncy/1e6 for k,v in iwall.items()}},open(sys.argv[1]+'.lab.json','w'))
