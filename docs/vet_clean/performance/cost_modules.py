#!/usr/bin/env python3
"""Module cost per cycle from the profiles of the four star runs -> cost_modules.csv.  cost_modules.py DATADIR [OUT.csv]
nsys stars (bsg, m35mms, m35zams): kname_r*.txt, gaps_r*.txt, mt.*.txt (weights), stats_r0.txt (split of the key "void").
t9l62: kt.*.txt (Kokkos kernel timer, a fence after every kernel) and mt.*.txt; nsys did not finish."""
import csv, glob, re, sys, collections
D = sys.argv[1]; OUT = sys.argv[2] if len(sys.argv) > 2 else None
C = collections.Counter
MODS = ['hydro', 'rays', 'rad-vet', 'others']
# cycles in the profile window: the header of kname / gaps counts two NewTimeStep kernel clusters per cycle and splits
# clusters at 200 ms, which fails for short cycles.  The true number = launches of CleanGdFaceSums#6 (one per cycle).
CELLS = {'bsg': 16777216, 'm35mms': 22020096, 'm35zams': 6291456, 't9l62': 67108864}

def kmod(k):
    """nsys kernel key -> (module, class)"""
    if k.startswith('MeshBoundaryValuesCC'): return 'others', 'comm'
    if k.startswith('fofc_flags::Exchange'): return 'hydro', 'comm'
    if k == 'bool' or re.search(r'VetGdHl(Begin|End)|VetGdHx|VetGdGbHalo|VetGdHalo', k): return 'rays', 'comm'
    if re.search(r'VetGd|VgdRing|VgdScan|CleanGdFaceSums|CleanClosure|CleanRayMoments|CleanRayCut|CleanRays|VetLatSweep', k): return 'rays', 'calc'
    if 'RadVetBCs' in k or k.startswith('radvet::'): return 'rad-vet', 'calc'
    if re.match(r'hydro::|fluid::|Coordinates|Star|HistoryOutput|GeneralHydro|FofcTrial|GravEtot|MeshBoundaryValues::HydroBCs', k): return 'hydro', 'calc'
    if k.startswith('Kokkos::Impl::View'): return 'others', 'unattr'
    return None

def lmod(s):
    """Kokkos label (kernel timer name, or the last kernel before an MPI call) -> (module, class of a kernel)"""
    if re.match(r'SendBuff|RecvBuff|Seam', s): return 'others', 'comm'
    if re.match(r'FOFC-owner|fofc_x', s): return 'hydro', 'comm'
    if re.match(r'cln_hd_|cln_ghost', s): return 'rad-vet', 'comm'
    if re.match(r'vet_vgd_(h[lx]|hc|halo|unpack|ipc|xff)|vet_vgd(gb|cs)_(?!shell|ring|estamp|eshell|list|rmask)', s): return 'rays', 'comm'
    if re.match(r'vet_vgd|cln_gd|cln_clositer|cln_raymom|cln_extrap|cln_raycut|cln_rays|cln_raynan|cln_fill|vet_vlat', s): return 'rays', 'calc'
    if re.match(r'cln_|cls_|vet_|radvetbc', s): return 'rad-vet', 'calc'
    if re.match(r'hyd|Hydro|hflux|h_update|dfl_|duale|FOFC|etotgrav|star_|Hist|mlat|wbcache|spsrc|hydrobc', s): return 'hydro', 'calc'
    return None

def key(s):     # the key of kname.py
    m = re.findall(r'([A-Za-z_][\w:]*)\([^()]*\)::\[lambda\([^\]]*\) \(instance (\d+)\)\]', s)
    m = [x for x in m if not x[0].endswith(('par_for', 'CleanFor', 'par_reduce', 'CleanReduce')) and 'Kokkos' not in x[0]]
    if m: return m[0][0] + '#' + m[0][1]
    m = re.search(r'ParallelFor<([\w:]+)', s) or re.search(r'Parallel\w+<([\w:]+)', s)
    return m.group(1) if m else s[:60]

def void_split(st):
    """parts of the key "void" (templated kernels), rank 0, ms per cycle over the 25 cycles; VetGdGbHaloC flags: start-up only"""
    o = C()
    L = open(f'{D}/{st}/stats_r0.txt').read().split('=== ')[1].split('\n')[2:]
    for r in csv.reader(L):
        if len(r) < 9 or key(r[8]) != 'void' or 'VetGdGbHaloC' in r[8]: continue
        v = int(r[1])/25e6; s = r[8]
        o['halo' if 'CleanHalo' in s else 'ghost' if 'CleanGhostPhys' in s else 'hydro' if 'hydro::' in s
          else 'dots' if 'CleanDots' in s else 'x:' + s[:200]] += v
    assert not [k for k in o if k.startswith('x:')], o
    return o

def startup(site, n):   # MPI site that is start-up or output: not a cycle cost
    if re.search(r'd_out_var|rst-', site): return True
    return site.startswith('Kokkos::View::initialization') and n < 26 or site == '(none)'

def mpi(st):
    """per rank: host MPI ms per cycle by module (26 cycles) and Isend bytes by module; unassigned sites listed"""
    T, B, una, drop = [], [], C(), C()
    for f in sorted(glob.glob(f'{D}/{st}/mt.*.txt')):
        if 'mt_rc' in f: continue
        t, b = C(), C()
        for l in open(f):
            p = l.rstrip('\n').split('\t')
            if len(p) < 3: continue
            s, n, by = p[0].split(); s = float(s); n = int(n); by = float(by)
            if p[1].startswith(('fence', 'deep_copy')): continue
            if startup(p[2], n): drop[p[1] + ' ' + p[2][:60]] += s/4; continue
            q = p[2]; m = re.search(r'\[(\w+)\]', q)
            if q.startswith('Kokkos::View::initialization') and m: q = m.group(1)
            a = lmod(q)
            if a is None: una[p[1] + ' ' + q[:60]] += s/26e-3/4
            a = a[0] if a else 'others'
            t[a] += s/26e-3
            if p[1] == 'Isend': b[a] += by/26
        T.append(t); B.append(b)
    return T, B, una, drop

def norm(c): s = sum(c.values()); return {k: c[k]/s for k in MODS}

def nsys(st):
    R = []; chk = []
    T, B, una, drop = mpi(st); vs = void_split(st); v0 = None; det = C()
    for r in range(4):
        L = open(f'{D}/{st}/kname_r{r}.txt').read().split('\n')
        h = re.search(r'cycles ([\d.]+) wall ([\d.]+) ms/cyc  busy ([\d.]+)  kernel sum ([\d.]+)', L[0])
        c, wall, busy, ks = map(float, h.groups())
        rows = []
        for l in L:
            p = l.split()
            if len(p) > 8 and re.match(r'^[\d.]+$', p[0]): rows.append((l.split('|')[-1].strip(), float(p[0]), float(p[2]), float(p[3])))
        g6 = next(x for x in rows if x[0].endswith('CleanGdFaceSums#6'))
        N = round(g6[1]*c*1e3/g6[3]); f = c/N                          # true cycles; factor to ms per true cycle
        wall, busy, ks = wall*f, busy*f, ks*f
        G = open(f'{D}/{st}/gaps_r{r}.txt').read().split('\n')
        idle = float(re.search(r'GPU idle ([\d.]+)', G[0]).group(1))*f
        act = {k: float(v)*f for k, v in re.findall(r'([a-z][\w-]*(?:\(part\))?|MPI(?:\(part\))?) ([\d.]+)', G[1].split(':', 1)[1])}
        mc = {int(k): float(v)*f for k, v in re.findall(r'(\d+): \d+, [\d.]+, ([\d.]+)', G[3].split(':', 1)[1].split('PtoP)')[1])}
        o = C(); det['ncyc'] = N
        for k, ms, n, us in rows:
            ms *= f
            if k == 'void':
                if v0 is None: v0 = ms
                w = ms/v0
                o['rad-vet', 'comm'] += (vs['halo'] + vs['ghost'])*w; o['hydro', 'calc'] += vs['hydro']*w
                o['rad-vet', 'calc'] += vs['dots']*w; rem = ms - sum(vs.values())*w; o['others', 'unattr'] += rem
                det['void: halo'] += vs['halo']*w/4; det['void: ghost'] += vs['ghost']*w/4; det['void: hydro'] += vs['hydro']*w/4
                det['void: dots'] += vs['dots']*w/4; det['void: not resolved'] += rem/4
                continue
            a = kmod(k)
            if a is None: print('UNMAPPED', st, k, ms); a = ('others', 'unattr')
            o[a] += ms
            if a[1] != 'calc': det[f'{a}: {k}'] += ms/4
        lost = ks - sum(o.values()); o['others', 'unattr'] += lost; det['kernel sum minus listed'] += lost/4
        # GPU copies outside kernels (wall = kernels + copies not under a kernel + idle)
        x = wall - busy - idle; pp = mc.get(10, 0)/sum(mc.values()); wb = norm(B[r]); wt = norm(T[r])
        im = act.get('MPI', 0) + act.get('MPI(part)', 0)
        for m in MODS: o[m, 'comm'] += x*pp*wb[m] + im*wt[m]
        o['others', 'unattr'] += x*(1 - pp); o['others', 'idle'] += idle - im
        det['copies outside kernels'] += x/4; det['of which peer (comm)'] += x*pp/4; det['idle total'] += idle/4; det['idle under MPI (comm)'] += im/4
        for k, v in act.items(): det['idle ' + k] += v/4
        det['busy'] += busy/4
        chk.append((wall, sum(o.values()))); R.append(o)
    return R, chk, det, T, B, una, drop

def ktimer(st, W):
    T, B, una, drop = mpi(st); R = []; det = C()
    det['ncyc'] = 25
    for r, f in enumerate(sorted(glob.glob(f'{D}/{st}/kt.*.txt'))):
        o = C()
        for l in open(f):
            m = re.match(r'\s*([\d.]+)\s+(\d+) ([FRSCX]) (.*)', l)
            if not m: continue
            s, n, kd, nm = float(m.group(1)), int(m.group(2)), m.group(3), m.group(4)
            ms = s*1e3/25
            if kd in 'CX': continue                      # copies contain their kernels, fences are the timer's own
            if nm.startswith('Kokkos::') or n < 25 and nm not in ('vet_t2_ffct_e', 'vet_vgd_shl'):
                det['dropped (start-up, output, generic copies): ' + (nm if ms > 2 else 'rest')] += ms/4; continue
            a = lmod(nm)
            if a is None: print('UNMAPPED', st, nm, ms); a = ('others', 'unattr')
            o[a] += ms; det['launches'] += n/25/4
            if a[1] != 'calc': det[f'{a}: {nm}'] += ms/4
        K = sum(o.values()); M = sum(T[r].values()); det['fenced kernels'] += K/4; det['host MPI'] += M/4
        for m in MODS: o[m, 'comm'] += T[r][m]
        det['plain wall minus (fenced kernels + host MPI)'] += (W - K - M)/4
        R.append(o)
    return R, [(sum(o.values()),)*2 for o in R], det, T, B, una, drop

def plain(st):
    t = open(f'{D}/{st}/plain.txt').read(); c = [float(x) for x in re.search(r'per cycle wall \[s\]: (.*)', t).group(1).split()][:-1]
    return 1e3*sum(c)/len(c), 1e3*sorted(c)[len(c)//2]        # last cycle (restart dump) left out

rowsout = []
for st in CELLS:
    pm, pmed = plain(st)
    R, chk, det, T, B, una, drop = nsys(st) if st != 't9l62' else ktimer(st, pm)
    W = sum(c[0] for c in chk)/4; S = sum(c[1] for c in chk)/4
    av = C()
    for o in R:
        for k, v in o.items(): av[k] += v/4
    print(f'== {st}: cycles {det["ncyc"]}  wall {W:.1f}  sum {S:.1f}  residual {W-S:+.2f} ms | plain mean {pm:.1f} median {pmed:.1f}')
    for m in MODS:
        tot = sum(v for k, v in av.items() if k[0] == m)
        print(f'  {m:8s} {tot:8.1f} {100*tot/W:5.1f}% ', '  '.join(f'{c} {av[m,c]:.1f} ({100*av[m,c]/W:.1f}%)' for c in ('calc', 'comm', 'unattr', 'idle') if abs(av[m, c]) > 0.005),
              ' ranks', ' '.join('%.1f' % sum(v for k, v in o.items() if k[0] == m) for o in R))
        for c in ('calc', 'comm', 'unattr', 'idle'):
            if c in ('calc', 'comm') or m == 'others':
                rowsout.append([st, CELLS[st], 'nsys' if st != 't9l62' else 'kernel timer (fenced)', m, {'unattr': 'unattributed', 'idle': 'GPU idle'}.get(c, {'calc': 'calculation', 'comm': 'communication'}.get(c)), f'{av[m,c]:.2f}', f'{100*av[m,c]/W:.2f}', f'{W:.1f}', f'{pm:.1f}'])
    for k, v in det.items(): print(f'     {v:9.2f}  {k}')
    print('   MPI host ms/cycle (mean of ranks):', {m: round(sum(t[m] for t in T)/4, 2) for m in MODS}, ' Isend MB/cycle:', {m: round(sum(b[m] for b in B)/4e6, 1) for m in MODS})
    print('   MPI sites not assigned (others), ms/cycle:', {k: round(v, 2) for k, v in una.most_common(6)})
    print('   MPI sites dropped (start-up / output), s:', {k: round(v, 3) for k, v in drop.most_common(4)})
if OUT:
    with open(OUT, 'w', newline='') as fo:
        w = csv.writer(fo); w.writerow(['star', 'cells', 'method', 'module', 'part', 'ms_per_cycle', 'percent', 'wall_ms_per_cycle', 'plain_run_ms_per_cycle']); w.writerows(rowsout)
