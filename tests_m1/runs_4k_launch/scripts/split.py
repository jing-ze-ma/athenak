# usage: python3 split.py <prof dir> <rank> <ncyc_total> [skipcyc=2]
# ms/cycle per category: kernel time + the idle gap BEFORE each kernel (attributed to the
# kernel that follows the gap).  Window = NewTimeStep of cycle `skipcyc` .. last NewTimeStep.
# bvals Pack/Unpack kernels count as M1 halo when the last non-bvals kernel was an
# ImplicitHaloCopy (the Krylov/transverse exchange), else as hydro.
import csv, re, sys, collections, glob
d, rk, ntot = sys.argv[1], sys.argv[2], int(sys.argv[3])
skip = int(sys.argv[4]) if len(sys.argv) > 4 else 2
def short(n):
    if 'M1PCRX<float>' in n: return 'M1PCRX<float>'
    if 'M1PCRX<double>' in n: return 'M1PCRX<double>'
    m = re.search(r'([A-Za-z_][\w:]*)\(\)[^:]*::\{lambda', n) or \
        re.search(r'(\w[\w:]*)\([^()]*\)::\{lambda', n)
    s = m.group(1) if m else n[:60]
    lam = re.search(r'::\{lambda\([^)]*\)#(\d+)\}', n)
    return s + ('#' + lam.group(1) if lam else '')
CAT = [('K:precond', r'ImplicitPCRSolve|ImplicitThomasSolve|ImplicitPrecond|M1PCRX'),
       ('K:operator', r'ImplicitApplyOp|ImplicitOffDiagOp|ImplicitODCache|ImplicitStencil|ImplicitOpX'),
       ('K:vec+red', r'ImplicitBiCGStab|ImplicitBcg|M1Bcg|ImplicitFloat'),
       ('K:halo', r'ImplicitHaloCopy|ImplicitKrylovHalo|ImplicitHaloDirect'),
       ('SC/tensor', r'Vet|ShortChar|vet'),
       ('Picard', r'RadiationM1::'),
       ('hydro+rest', r'.')]
kt = glob.glob(f'{d}/**/rank_{rk}_kernel_trace.csv', recursive=True)[0]
rows = []
for r in csv.DictReader(open(kt)):
    rows.append((int(r['Start_Timestamp']), int(r['End_Timestamp']), short(r['Kernel_Name'])))
rows.sort()
T = [r[0] for r in rows if 'NewTimeStep' in r[2]]
per = max(1, round(len(T)/ntot))
t0, t1 = T[per*skip], T[-1]
w = [r for r in rows if t0 <= r[0] < t1]
ncyc = len([t for t in T if t0 <= t < t1])/per
span = t1 - t0
kern = collections.Counter(); idle = collections.Counter(); cnt = collections.Counter()
byk = collections.Counter(); byn = collections.Counter()
last = w[0][0]; lastnb = ''
for a, b, s in w:
    if 'PackAndSend' in s or 'RecvAndUnpack' in s:
        c = 'K:halo' if 'ImplicitHaloCopy' in lastnb else 'hydro+rest'
    else:
        c = next(k for k, p in CAT if re.search(p, s))
        lastnb = s
    kern[c] += b - a; cnt[c] += 1
    gap = a - last
    if gap > 0: idle[c] += gap
    last = max(last, b)
    byk[s] += b - a; byn[s] += 1
f = 1e6*ncyc
print(f'{d}: rank {rk}, {ncyc:.1f} cycles, {span/f:.1f} ms/cycle, launches/cycle {len(w)/ncyc:.0f}, '
      f'busy {sum(kern.values())/f:.1f}, idle {sum(idle.values())/f:.1f}')
print(f'  {"category":12s} {"kernel ms":>9s} {"idle-before":>11s} {"total":>7s} {"launch/cyc":>10s}')
for k, _ in CAT:
    print(f'  {k:12s} {kern[k]/f:9.2f} {idle[k]/f:11.2f} {(kern[k]+idle[k])/f:7.2f} {cnt[k]/ncyc:10.1f}')
print('  top kernels (ms/cycle, calls/cycle, us/call):')
for s, v in byk.most_common(25):
    print(f'    {v/f:7.2f} {byn[s]/ncyc:7.1f} {v/byn[s]/1e3:8.1f}  {s[:90]}')
print('  top kernels by launches (calls/cycle, ms/cycle):')
for s, v in byn.most_common(30):
    print(f'    {v/ncyc:7.1f} {byk[s]/f:7.2f}  {s[:90]}')
ht = glob.glob(f'{d}/**/rank_{rk}_hip_api_trace.csv', recursive=True)
if ht:
    SY = ('hipStreamSynchronize', 'hipDeviceSynchronize', 'hipEventSynchronize',
          'hipMemcpy', 'hipMemcpyAsync', 'hipMemcpyDtoH', 'hipMemcpyHtoD', 'hipStreamWaitEvent',
          'hipEventQuery', 'hipStreamQuery', 'hipMemsetAsync', 'hipMemset')
    c = collections.Counter(); tt = collections.Counter()
    for r in csv.DictReader(open(ht[0])):
        a = int(r['Start_Timestamp'])
        if t0 <= a < t1:
            fn = r['Function']
            c[fn] += 1; tt[fn] += int(r['End_Timestamp']) - a
    sy = sum(c[k] for k in ('hipStreamSynchronize', 'hipDeviceSynchronize', 'hipEventSynchronize'))
    print(f'  host syncs/cycle (Stream+Device+Event) {sy/ncyc:.1f}')
    for k, v in c.most_common(20):
        print(f'    {k:32s} {v/ncyc:8.1f}/cyc {tt[k]/f:7.2f} ms/cyc')
