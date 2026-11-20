#!/usr/bin/env python3
# gaps.py SQLITE C0 NC SUBCYC: GPU idle time in the cycle window and what the host thread was doing meanwhile
# (MPI call / cuda sync / memcpy API / kernel launch API / no traced call = host compute); gap-size histogram; memcpy totals.
import sys, os, collections, bisect; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import open_db, key, window
c, S = open_db(sys.argv[1]); c0, nc, sub = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
w0, w1, ncy = window(c, S, c0, nc, sub)
tabs = {r[0] for r in c.execute("select name from sqlite_master where type='table'")}
act = [(st, en, 'K', S.get(dn, '?')) for st, en, dn in c.execute("select start,end,demangledName from CUPTI_ACTIVITY_KIND_KERNEL where end>? and start<?", (w0, w1))]
mc = []
if 'CUPTI_ACTIVITY_KIND_MEMCPY' in tabs:
    mc = c.execute("select start,end,bytes,copyKind from CUPTI_ACTIVITY_KIND_MEMCPY where end>? and start<?", (w0, w1)).fetchall()
    act += [(st, en, 'M', str(k)) for st, en, b, k in mc]
if 'CUPTI_ACTIVITY_KIND_MEMSET' in tabs:
    act += [(st, en, 'S', '') for st, en in c.execute("select start,end from CUPTI_ACTIVITY_KIND_MEMSET where end>? and start<?", (w0, w1))]
act.sort()
# host intervals: CUDA runtime API (by name) and MPI (tables or NVTX ranges named MPI_*)
host = []
for st, en, nid in c.execute("select start,end,nameId from CUPTI_ACTIVITY_KIND_RUNTIME where end>? and start<?", (w0, w1)):
    n = S.get(nid, '?')
    cat = 'sync' if 'Synchronize' in n or 'EventQuery' in n else ('memcpy-api' if 'emcpy' in n or 'emset' in n else ('launch' if 'Launch' in n else 'cuda-other'))
    host.append((st, en, cat, n))
for t in tabs:
    if t.startswith('MPI_') and t.endswith('EVENTS'):
        try: host += [(st, en, 'MPI', t) for st, en in c.execute(f"select start,end from {t} where end>? and start<?", (w0, w1))]
        except Exception: pass
if 'NVTX_EVENTS' in tabs:
    for st, en, tx, tid in c.execute("select start,end,text,textId from NVTX_EVENTS where end>? and start<? and end is not null", (w0, w1)):
        n = tx or S.get(tid, '') or ''
        if n.startswith('MPI_'): host.append((st, en, 'MPI', n))
host.sort(); hs = [h[0] for h in host]
gaps = []; e0 = w0; prev = None
for st, en, ty, nm in act:
    if st > e0: gaps.append((e0, st, prev, nm if ty == 'K' else ty))
    if en > e0: e0 = en; prev = nm if ty == 'K' else ty
if w1 > e0: gaps.append((e0, w1, prev, 'end'))
prio = ['MPI', 'memcpy-api', 'sync', 'launch', 'cuda-other']
bycat = collections.Counter(); bysz = collections.Counter(); nsz = collections.Counter(); big = []
for g0, g1, a, b in gaps:
    L = g1 - g0; ov = collections.Counter(); i = max(bisect.bisect_left(hs, g0) - 200, 0)
    while i < len(host) and host[i][0] < g1:
        h0, h1, cat, n = host[i]
        o = min(h1, g1) - max(h0, g0)
        if o > 0: ov[cat] += o
        i += 1
    cat = next((p for p in prio if ov[p] > 0.5*L), 'host-compute' if not ov else max(ov, key=ov.get) + '(part)')
    bycat[cat] += L
    sz = '<10us' if L < 1e4 else ('10-100us' if L < 1e5 else ('0.1-1ms' if L < 1e6 else '>1ms'))
    bysz[sz] += L; nsz[sz] += 1; big.append((L, cat, key(a)[0] if a else '-', key(b)[0] if b else '-'))
wall = (w1 - w0)/ncy/1e6; idle = sum(bysz.values())/ncy/1e6
print(f'{sys.argv[1]}: cycles {ncy:g} wall {wall:.1f} ms/cyc  GPU idle {idle:.1f} ms/cyc ({100*idle/wall:.1f} %)  gaps {len(gaps)/ncy:.0f}/cyc')
print('idle by host activity (ms/cyc): ' + '  '.join(f'{k} {v/ncy/1e6:.1f}' for k, v in bycat.most_common()))
print('idle by gap size   (ms/cyc, n/cyc): ' + '  '.join(f'{k} {bysz[k]/ncy/1e6:.1f} ({nsz[k]/ncy:.0f})' for k in ['<10us', '10-100us', '0.1-1ms', '>1ms']))
mk = collections.defaultdict(lambda: [0, 0, 0])
for st, en, b, k in mc: x = mk[k]; x[0] += 1; x[1] += b; x[2] += en - st
print('memcpy (kind: n/cyc, MB/cyc, ms/cyc; 1=HtoD 2=DtoH 8=DtoD 10=PtoP): ' + '  '.join(f'{k}: {v[0]/ncy:.0f}, {v[1]/ncy/1e6:.1f}, {v[2]/ncy/1e6:.2f}' for k, v in sorted(mk.items())))
hc = collections.Counter()
for h0, h1, cat, n in host:
    if cat == 'MPI' or cat == 'sync': hc[cat] += min(h1, w1) - max(h0, w0)
print('host time in (ms/cyc, overlapping GPU work or not): ' + '  '.join(f'{k} {v/ncy/1e6:.1f}' for k, v in hc.items()))
print('largest gaps (us, host activity, previous kernel -> next):')
for L, cat, a, b in sorted(big, reverse=True)[:15]: print(f'  {L/1e3:9.1f} {cat:14s} {a} -> {b}')
