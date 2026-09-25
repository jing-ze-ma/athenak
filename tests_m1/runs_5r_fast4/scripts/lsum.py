"""lsum.py <prof dir> [ncyc=40] [skip=10] [top=40]: kernel time per Kokkos label (kokkos-trace
renamed kernels), ms/cycle and us/call over the last ncyc-skip cycles (NewTimeStep marks)."""
import csv, glob, sys, collections
d = sys.argv[1]
nt = int(sys.argv[2]) if len(sys.argv) > 2 else 40
sk = int(sys.argv[3]) if len(sys.argv) > 3 else 10
top = int(sys.argv[4]) if len(sys.argv) > 4 else 40
kt = glob.glob(d + '/**/rank_0_kernel_trace.csv', recursive=True)[0]
rows = sorted((int(r['Start_Timestamp']), int(r['End_Timestamp']), r['Kernel_Name'])
              for r in csv.DictReader(open(kt)))
MK = sys.argv[5] if len(sys.argv) > 5 else 'newdt'
ts = [r[0] for r in rows if MK in r[2]]
per = max(1, round(len(ts)/nt))
t0 = ts[len(ts) - 1 - per*(nt - sk)]
w = [r for r in rows if r[0] > t0]
nc = nt - sk
T = collections.Counter(); N = collections.Counter(); I = collections.Counter()
last = w[0][0]
for a, b, n in w:
    T[n] += b - a; N[n] += 1; I[n] += max(0, a - last); last = max(last, b)
tot = (w[-1][1] - w[0][0])/nc/1e6
print(f'{d}: {tot:.2f} ms/cycle, busy {sum(T.values())/nc/1e6:.2f}')
for n, t in T.most_common(top):
    print(f'{t/nc/1e6:7.3f} {I[n]/nc/1e6:6.3f} {N[n]/nc:6.1f} {t/N[n]/1e3:8.1f}  {n[:60]}')
