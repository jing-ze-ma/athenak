#!/usr/bin/env python3
# kname.py SQLITE [C0] [NC] [NTOP] [SUBCYC]: rank-0 GPU kernel table keyed by demangled enclosing function + lambda instance
# (works without kp_nvtx_connector). Cycle ends as r0lab.py (Hydro::NewTimeStep kernel clusters); SUBCYC = markers per cycle (default 2).
# Columns: ms/cycle, launches/cycle, us/launch, regs/thread, block, grid(max), local mem B/thread, NVTX label (if present), key
import sys, sqlite3, re, bisect, collections
db = sqlite3.connect(sys.argv[1]); c = db.cursor()
c0 = int(sys.argv[2]) if len(sys.argv) > 2 else 1; nc = int(sys.argv[3]) if len(sys.argv) > 3 else 3
NT = int(sys.argv[4]) if len(sys.argv) > 4 else 80; sub = float(sys.argv[5]) if len(sys.argv) > 5 else 2.0
S = dict(c.execute("select id,value from StringIds"))
def key(dn):
    s = S.get(dn, '?')
    m = re.findall(r'([A-Za-z_][\w:]*)\([^()]*\)::\[lambda\([^\]]*\) \(instance (\d+)\)\]', s)
    m = [x for x in m if not x[0].endswith(('par_for', 'CleanFor', 'par_reduce', 'CleanReduce')) and 'Kokkos' not in x[0]]
    if m: f, i = m[0]; return f.replace('radm1::RadiationM1::', 'M1::') + '#' + i
    m = re.search(r'ParallelFor<([\w:]+)', s) or re.search(r'Parallel\w+<([\w:]+)', s)
    return (m.group(1) if m else s[:60])
kn = c.execute("select start,end,correlationId,demangledName,registersPerThread,blockX*blockY*blockZ,gridX*gridY*gridZ,localMemoryPerThread from CUPTI_ACTIVITY_KIND_KERNEL").fetchall()
mk = sorted(st for st, en, cid, dn, *_ in kn if 'NewTimeStep' in S.get(dn, ''))
cl = []
for t in mk:
    if not cl or t - cl[-1][-1] > 2e8: cl.append([t])
    else: cl[-1].append(t)
ends = [x[0] for x in cl]; j1 = min(c0 + nc, len(ends) - 1); w0, w1 = ends[c0], ends[j1]; ncy = (j1 - c0)/sub
nv = c.execute("select start,end,coalesce(text,''),textId,globalTid from NVTX_EVENTS where end is not null and eventType=59").fetchall()
rng = collections.defaultdict(list)
for st, en, tx, tid, g in nv: rng[g].append((st, en, (tx or S.get(tid, '')).lstrip(':')))
for g in rng: rng[g].sort()
sts = {g: [r[0] for r in rng[g]] for g in rng}
rt = {cid: (st, g) for cid, st, g in c.execute("select correlationId,start,globalTid from CUPTI_ACTIVITY_KIND_RUNTIME where start between ? and ?", (w0 - 10**9, w1))} if nv else {}
def label(cid):
    if cid not in rt: return ''
    st, g = rt[cid]; L = rng.get(g, []); i = bisect.bisect_right(sts.get(g, []), st) - 1; lo = i - 300
    while i >= 0 and i > lo:
        if L[i][0] <= st <= L[i][1]: return L[i][2]
        i -= 1
    return ''
tot = collections.Counter(); cnt = collections.Counter(); info = {}; lab = collections.defaultdict(collections.Counter); busy = []
for st, en, cid, dn, rg, bl, gr, lm in kn:
    if st < w0 or st >= w1: continue
    k = key(dn); tot[k] += en - st; cnt[k] += 1; busy.append((st, en))
    o = info.get(k); info[k] = (rg, bl, max(gr, o[2] if o else 0), lm)
    if nv: lab[k][label(cid)] += en - st
busy.sort(); u = 0; e0 = w0
for s, e in busy:
    if e > e0: u += e - max(s, e0); e0 = e
T = sum(tot.values())
print(f'{sys.argv[1]}: cycles {ncy:g} wall {(w1-w0)/ncy/1e6:.1f} ms/cyc  busy {u/ncy/1e6:.1f}  kernel sum {T/ncy/1e6:.1f}  launches {sum(cnt.values())/ncy:.0f}/cyc  keys {len(tot)}')
print(f'{"ms/cyc":>8} {"%":>5} {"n/cyc":>7} {"us/l":>7} {"reg":>4} {"blk":>4} {"grid":>8} {"lmem":>4}  label | key')
for k, v in tot.most_common(NT):
    rg, bl, gr, lm = info[k]; L = ','.join(x for x, _ in lab[k].most_common(3) if x) if nv else ''
    print(f'{v/ncy/1e6:8.2f} {100*v/T:5.1f} {cnt[k]/ncy:7.0f} {v/cnt[k]/1e3:7.1f} {rg:4d} {bl:4d} {gr:8d} {lm:4d}  {L[:40]} | {k}')
