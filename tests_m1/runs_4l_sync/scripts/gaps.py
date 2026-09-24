# usage: python3 gaps.py <prof dir> <rank> <ncyc_total> [skip=10] [top=25]
# GPU-idle gaps per cycle, attributed to (kernel before -> kernel after) and to the host
# sync/memcpy API call that covers most of the gap.  Also: host sync time that overlaps
# GPU work (free) vs GPU-idle (real cost), per API function.
import csv, re, sys, collections, glob, bisect
sys.path.insert(0, __file__.rsplit('/', 1)[0])
d, rk, ntot = sys.argv[1], sys.argv[2], int(sys.argv[3])
skip = int(sys.argv[4]) if len(sys.argv) > 4 else 10
top = int(sys.argv[5]) if len(sys.argv) > 5 else 25
def short(n):
    if 'M1PCRX' in n: return 'M1PCRX'
    m = re.search(r'([A-Za-z_][\w:]*)\(\)[^:]*::\{lambda', n) or \
        re.search(r'(\w[\w:]*)\([^()]*\)::\{lambda', n)
    s = m.group(1) if m else n[:50]
    lam = re.search(r'::\{lambda\([^)]*\)#(\d+)\}', n)
    s = s.replace('radm1::RadiationM1::', 'M1::')
    return s + ('#' + lam.group(1) if lam else '')
kt = glob.glob(f'{d}/**/rank_{rk}_kernel_trace.csv', recursive=True)[0]
K = sorted((int(r['Start_Timestamp']), int(r['End_Timestamp']), short(r['Kernel_Name']))
           for r in csv.DictReader(open(kt)))
T = [a for a, b, s in K if 'NewTimeStep' in s]
per = max(1, round(len(T)/ntot))
t0, t1 = T[per*skip], T[-1]
ncyc = len([t for t in T if t0 <= t < t1])/per
f = 1e6*ncyc
W = [k for k in K if t0 <= k[0] < t1]
# busy intervals (merged)
bus = []
for a, b, s in W:
    if bus and a <= bus[-1][1]: bus[-1][1] = max(bus[-1][1], b)
    else: bus.append([a, b])
gaps = []   # (start, end, prev kernel, next kernel)
last = W[0][1]; prev = W[0][2]; lastend = W[0][1]
for a, b, s in W[1:]:
    if a > lastend: gaps.append((lastend, a, prev, s))
    if b >= lastend: lastend = b; prev = s
SY = {'hipStreamSynchronize', 'hipDeviceSynchronize', 'hipEventSynchronize',
      'hipMemcpyAsync', 'hipMemcpy', 'hipMemcpyDtoH', 'hipMemcpyHtoD'}
ht = glob.glob(f'{d}/**/rank_{rk}_hip_api_trace.csv', recursive=True)[0]
A = sorted((int(r['Start_Timestamp']), int(r['End_Timestamp']), r['Function'])
           for r in csv.DictReader(open(ht)) if r['Function'] in SY)
A = [x for x in A if t0 <= x[0] < t1]
bs = [x[0] for x in bus]
def busy_in(a, b):
    i = max(0, bisect.bisect_right(bs, a) - 1); tot = 0
    while i < len(bus) and bus[i][0] < b:
        tot += max(0, min(b, bus[i][1]) - max(a, bus[i][0])); i += 1
    return tot
tot = collections.Counter(); idl = collections.Counter(); n = collections.Counter()
for a, b, fn in A:
    tot[fn] += b - a; idl[fn] += (b - a) - busy_in(a, b); n[fn] += 1
print(f'{d} rank {rk}: {ncyc:.1f} cyc, span {(t1-t0)/f:.1f} ms/cyc, '
      f'idle {sum(g[1]-g[0] for g in gaps)/f:.2f} ms/cyc in {len(gaps)/ncyc:.0f} gaps/cyc')
print('  host API          calls/cyc  host-wait ms  of which GPU idle ms')
for fn in sorted(tot, key=lambda k: -idl[k]):
    print(f'  {fn:24s} {n[fn]/ncyc:7.1f} {tot[fn]/f:10.2f} {idl[fn]/f:10.2f}')
As = [x[0] for x in A]
g2 = collections.Counter(); gn = collections.Counter()
for a, b, p, q in gaps:
    # the API sync covering most of the gap
    i = max(0, bisect.bisect_left(As, a) - 50); best = ('launch/host', 0)
    while i < len(A) and A[i][0] < b:
        ov = min(b, A[i][1]) - max(a, A[i][0])
        if ov > best[1]: best = (A[i][2].replace('hip', '').replace('Synchronize', 'Sync'), ov)
        i += 1
    key = (p, q, best[0])
    g2[key] += b - a; gn[key] += 1
print(f'  top GPU-idle gaps: ms/cyc  n/cyc  kernel before -> kernel after [host call in gap]')
for k, v in g2.most_common(top):
    print(f'  {v/f:6.2f} {gn[k]/ncyc:6.1f}  {k[0][:38]} -> {k[1][:38]} [{k[2]}]')
