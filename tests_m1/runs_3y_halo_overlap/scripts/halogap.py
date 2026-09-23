# GPU timeline around implicit_halo_mpi exchanges (rocprofv3 kernel trace, one rank):
# for each m1_impl_hm_pack ... m1_impl_hm_unpack, the GPU idle time between the last
# kernel before the unpack and the unpack, and the duration of the next kernels (the
# operator) after the unpack.  usage: halogap.py <kernel_trace.csv> [t0frac]
import sys, csv, re, collections
def short(n):
    m = re.findall(r'RadiationM1::(\w+)', n)
    return m[-1] if m else n[:40]
rows = []
for r in csv.DictReader(open(sys.argv[1])):
    rows.append((int(r['Start_Timestamp']), int(r['End_Timestamp']), short(r['Kernel_Name'])))
rows.sort()
n = len(rows); lo = n//2   # second half: steady state
gaps = []; ops = collections.Counter(); opt = collections.Counter(); busy = 0
tot_idle = 0
for a in range(lo+1, n):
    idle = rows[a][0] - max(r[1] for r in rows[max(lo, a-3):a])
    if idle > 0: tot_idle += idle
for a in range(lo, n):
    if rows[a][2] == 'ImplicitHaloMPI' and rows[a-1][2] == 'ImplicitHaloDirect':
        prev = rows[a-1]
        gaps.append((rows[a][0] - prev[1])/1e3)
        nx = rows[a+1]
        ops[nx[2]] += 1; opt[nx[2]] += (nx[1]-nx[0])/1e3
span = (rows[-1][1] - rows[lo][0])/1e6
print('exchanges %d, span %.1f ms, total idle %.1f ms' % (len(gaps), span, tot_idle/1e6))
gaps.sort()
if gaps:
    print('idle before unpack us: median %.1f mean %.1f p90 %.1f sum %.1f ms' % (
        gaps[len(gaps)//2], sum(gaps)/len(gaps), gaps[int(.9*len(gaps))], sum(gaps)/1e3))
for k, c in ops.most_common(6):
    print('  after unpack: %-40s n=%d mean %.1f us' % (k, c, opt[k]/c))
# idle classified by (previous kernel -> next kernel)
cls = collections.Counter(); cnt = collections.Counter()
for a in range(lo+1, n):
    idle = rows[a][0] - max(r[1] for r in rows[max(lo, a-3):a])
    if idle > 0:
        key = rows[a-1][2][:28] + ' -> ' + rows[a][2][:28]
        cls[key] += idle; cnt[key] += 1
print('largest idle classes (ms over the span, count):')
for k, v in cls.most_common(int(sys.argv[2]) if len(sys.argv) > 2 else 12):
    print('  %6.1f %6d  %s' % (v/1e6, cnt[k], k))
