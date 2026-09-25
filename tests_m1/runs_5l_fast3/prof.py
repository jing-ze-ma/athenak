# usage: python3 prof.py <prof dir> <rank> <ncyc_total> <skipcyc>
# ms/cycle per category of the implicit M1 cycle (kernel time + GPU idle before each
# kernel, attributed to the kernel that follows the gap), top kernels, HIP sync calls.
# Window = NewTimeStep of cycle `skipcyc` .. last NewTimeStep (as split.py).
import csv, re, sys, collections, glob
d, rk, ntot, skip = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])


def encl(n):
    # the function enclosing the first lambda: walk back over its balanced (...)
    p = n.find('::{lambda')
    if p < 0:
        return n[:50]
    q = p - 1
    while q >= 0 and n[q] != ')':
        q -= 1
    dep = 0
    while q >= 0:
        if n[q] == ')':
            dep += 1
        elif n[q] == '(':
            dep -= 1
            if dep == 0:
                break
        q -= 1
    m = re.search(r'([A-Za-z_][\w:]*)(?:<[^<>]*(?:<[^<>]*>[^<>]*)*>)?$', n[:q])
    return m.group(1).split('::')[-1] if m else n[:50]


def short(n):
    if 'M1PCRX<float>' in n:
        return 'M1PCRX<float>'
    if 'M1PCRX<double>' in n:
        return 'M1PCRX<double>'
    m = re.search(r'RadiationM1::(\w+)\(', n)
    s = m.group(1) if m else encl(n)
    lam = re.search(r'::\{lambda\([^)]*\)#(\d+)\}', n)
    red = 'R' if 'ParallelReduce' in n else ''
    return s + ('#' + lam.group(1) if lam else '') + red


CAT = [('KRY:precond', r'PCRSolve|ThomasSolve|ImplicitPrecond|M1PCRX'),
       ('KRY:matvec+dot', r'ImplicitApplyOp|ImplicitOffDiagOp|ImplicitStencilOp|'
                          r'ImplicitOpX|ImplicitHaloOp'),
       ('KRY:vec+red', r'ImplicitBiCGStab|ImplicitBcg|M1Bcg|ImplicitFloat'),
       ('KRY:halo', r'ImplicitHaloCopy|ImplicitKrylovHalo|ImplicitHaloDirect|'
                    r'ImplicitHaloMPI|ImplicitHaloExchange|ImplicitTransverseHalo'),
       ('VET build', r'Vet|vet|ShortChar'),
       ('Picard setup+rest M1', r'RadiationM1::|radm1'),
       ('hydro+rest', r'.')]
kt = glob.glob(f'{d}/**/rank_{rk}_kernel_trace.csv', recursive=True)[0]
rows = []
for r in csv.DictReader(open(kt)):
    n = r['Kernel_Name']
    rows.append((int(r['Start_Timestamp']), int(r['End_Timestamp']), short(n), n))
rows.sort()
T = [r[0] for r in rows if 'NewTimeStep' in r[2]]
per = max(1, round(len(T)/ntot))
# the NewTimeStep kernels end a cycle (per of them per cycle): the window is the last
# ntot-skip cycles, from the last NewTimeStep of cycle `skip` to the last one of the run
t0, t1 = T[len(T) - 1 - per*(ntot - skip)], T[-1] + 1
w = [r for r in rows if t0 < r[0] < t1]
ncyc = ntot - skip
span = t1 - t0
kern = collections.Counter(); idle = collections.Counter(); cnt = collections.Counter()
byk = collections.Counter(); byn = collections.Counter(); cat_of = {}
redwait = 0.0; nredwait = 0
last = w[0][0]; lastnb = ''; lastname = ''
for a, b, s, n in w:
    if 'PackAndSend' in s or 'RecvAndUnpack' in s:
        c = 'KRY:halo' if 'Halo' in lastnb else 'hydro+rest'
    else:
        c = next(k for k, p in CAT if re.search(p, n))
        lastnb = s
    cat_of[s] = c
    g = max(0, a - last)
    if lastname.startswith('ImplicitStencilOp') and lastname.endswith('R') \
       and c.startswith('KRY'):
        redwait += g; nredwait += 1
    kern[c] += b - a; idle[c] += g; cnt[c] += 1
    byk[s] += b - a; byn[s] += 1
    last = max(last, b); lastname = s
busy = sum(kern.values())
print(f'{d}: rank {rk}, {ncyc:.1f} cycles, {span/ncyc/1e6:.2f} ms/cycle, launches/cycle '
      f'{len(w)/ncyc:.0f}, busy {busy/ncyc/1e6:.2f}, idle {(span-busy)/ncyc/1e6:.2f}')
print(f'  {"category":22s} kernel_ms idle_before  total launch/cyc')
for c, _ in CAT:
    print(f'  {c:22s} {kern[c]/ncyc/1e6:8.2f} {idle[c]/ncyc/1e6:10.2f} '
          f'{(kern[c]+idle[c])/ncyc/1e6:7.2f} {cnt[c]/ncyc:8.1f}')
print(f'  idle after a reducing stencil op before the next Krylov kernel: '
      f'{redwait/ncyc/1e6:.2f} ms/cycle over {nredwait/ncyc:.1f} gaps')
print('  top kernels (ms/cycle, calls/cycle, us/call, category):')
for s, t in byk.most_common(30):
    print(f'   {t/ncyc/1e6:7.3f} {byn[s]/ncyc:7.1f} {t/byn[s]/1e3:8.1f}  {s}  [{cat_of[s]}]')
ht = glob.glob(f'{d}/**/rank_{rk}_hip_api_trace.csv', recursive=True)
if ht:
    hs = collections.Counter(); hn = collections.Counter()
    for r in csv.DictReader(open(ht[0])):
        a = int(r['Start_Timestamp'])
        if t0 <= a < t1:
            f = r['Function']
            hs[f] += int(r['End_Timestamp']) - a; hn[f] += 1
    print('  HIP API (host ms/cycle, calls/cycle):')
    for f, t in hs.most_common(8):
        print(f'   {t/ncyc/1e6:7.2f} {hn[f]/ncyc:7.1f}  {f}')
