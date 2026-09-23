# usage: python3 vetprof.py <kernel_trace.csv> : per-call split of the vet_sc formal solution
import csv, sys, collections, re
R = list(csv.DictReader(open(sys.argv[1])))
R.sort(key=lambda r: int(r['Start_Timestamp']))
def tag(n):
    if 'Vet' not in n: return None
    m = re.search(r'radm1::(?:RadiationM1::|\(anonymous namespace\)::)(\w+)', n)
    base = m.group(1) if m else 'vet?'
    l = re.search(r'\{lambda\(([^)]*)\)#(\d+)\}', n)
    return base + (('#' + str(l.group(1).count(',')+1) + 'a' + l.group(2)) if l else '')
calls = []; cur = None
for r in R:
    t = tag(r['Kernel_Name'])
    s, e = int(r['Start_Timestamp']), int(r['End_Timestamp'])
    if t is None:
        if cur is not None and cur['open']: cur['other'] += e - s
        continue
    if t in ('VetShortChar#3a1', 'VetShortChar#1a1') or (cur is None):
        cur = dict(t0=s, t1=e, k=collections.Counter(), n=collections.Counter(), other=0,
                   open=True); calls.append(cur)
    cur['t1'] = e; cur['k'][t] += e - s; cur['n'][t] += 1
    if t.startswith('VetFullTensor#3') or t.startswith('VetShortChar#2'):
        pass
calls = calls[2:]  # skip warm-up
nc = len(calls)
span = sum(c['t1'] - c['t0'] for c in calls)/nc/1e6
tot = collections.Counter(); cnt = collections.Counter()
for c in calls:
    tot.update(c['k']); cnt.update(c['n'])
busy = sum(tot.values())/nc/1e6
print(f'calls {nc}: span first..last vet kernel {span:.2f} ms, vet kernel busy {busy:.2f} ms, '
      f'non-vet kernels inside {sum(c["other"] for c in calls)/nc/1e6:.2f} ms')
for k, v in sorted(tot.items(), key=lambda x: -x[1]):
    print(f'  {k:28s} {v/nc/1e6:8.3f} ms/call  {cnt[k]/nc:7.1f} launches/call  '
          f'{v/max(cnt[k],1)/1e3:8.2f} us/launch')
