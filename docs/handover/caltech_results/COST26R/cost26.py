# cost26.py: 26r/26s tables per star from kt (named kernels, fenced, 25 cycles, mean of 4 ranks), mt (MPI per site, 26 cycles), nsys kname/gaps (3 stars), plain timing
import glob, re, os, sys, json, numpy as np, collections
S = '/resnick/groups/carnegie_poc/jingze/b1_1124/prof26r/'
RAY = re.compile(r'^(vet_vgd|vet_vlat|vet_lay|cln_gd|cln_raymom|cln_rayj|cln_raycut|vet_t2_avg|cln_clositer)')
HYD = re.compile(r'^(hflux|hyd|Hydro|h_update|hydro_|FOFC|dfl|duale|wbcache|spsrc|star_|etotgrav|HistSums|mlat|fofc|hst_)')
VET = re.compile(r'^(cln_|cls_|vet_)')
COMM = re.compile(r'^(cln_hd_|cln_ghost|vet_vgd_hl_|vet_vgd_hxb|vet_vgd_halo|vet_vgd_unpack|vet_vgd_shl|SendBuff|RecvBuff|Pack|Unpack|fofc_xc)')
AMB = re.compile(r'^(cln_cslat|cln_cstie|cln_tiet|cln_fwb|vet_vgd_ring|cln_slot|cln_fs3_sh|cln_apply_st_sh|cln_cslat_sh)')
def mod(n):
    if n.startswith('Kokkos::ViewCopy') or n.startswith('Kokkos::View') or n.startswith('Kokkos::Impl') or n.startswith('rst-') or '->' in n: return 'others'
    if RAY.match(n): return 'rays'
    if HYD.match(n): return 'hydro'
    if VET.match(n): return 'rad-vet'
    return 'others'
def rdkt(st, d='kt'):
    R = []
    for f in sorted(glob.glob(S + st + '/' + d + '/kt.*.txt')):
        o = []
        for l in open(f):
            p = l.split(None, 3)
            if l[0] == '#' or len(p) < 4: continue
            o.append((float(p[0]), int(p[1]), p[2], p[3].strip()))
        R.append(o)
    return R
for st in sys.argv[1:]:
    R = rdkt(st); R6 = rdkt(st, 'kt6'); print('\n### %s' % st); sub = len(R6) == len(R) and len(R) > 0
    if sub: R = [[(t - d6.get((k, nm), (0.0, 0))[0], n, k, nm) for t, n, k, nm in o] for o, d6 in zip(R, [{(k, nm): (t, n) for t, n, k, nm in o} for o in R6])]
    # kernels (F and R rows), C = copies, X = fences
    tot = collections.defaultdict(list); amb = collections.Counter(); names = collections.Counter()
    for o in R:
        acc = collections.Counter()
        for t, n, k, nm in o:
            if k in 'FRS':
                m = mod(nm); c = 'comm' if COMM.match(nm) else ('copy' if nm.startswith('Kokkos::View') else 'calc'); acc[(m, c)] += t
                names[(m, c, nm)] += t
                if AMB.match(nm): amb[nm] += t
            elif k == 'C': acc[('others', 'copy(C rows)')] += t
        for x in set(acc): tot[x].append(acc[x])
    nc = 19.0 if sub else 25.0; K = sum(v for x, v in [(x, np.mean(tot[x])) for x in tot] if x[1] != 'copy(C rows)')
    tab = {x: 1e3*np.mean(v)/nc for x, v in tot.items()}
    print('kernel timer (fenced, %d ranks, %s): ms per cycle, mean of the ranks' % (len(R), '25-cycle minus 6-cycle run, per 19 cycles (start-up and end removed)' if sub else 'all 25 cycles incl. the restart set-up kernels'))
    mods = ('hydro', 'rays', 'rad-vet', 'others')
    print('| module | calculation | pack / copy / unpack / ghost / request | Kokkos view copies, init | sum | share of kernel time |\n|---|---|---|---|---|---|')
    s = {m: sum(v for (a, b), v in tab.items() if a == m and b != 'copy(C rows)') for m in mods}; SS = sum(s.values())
    for m in mods: print('| %s | %.2f | %.2f | %.2f | %.2f | %.1f %% |' % (m, tab.get((m, 'calc'), 0), tab.get((m, 'comm'), 0), tab.get((m, 'copy'), 0), s[m], 100*s[m]/SS))
    print('- kernel sum %.1f ms per cycle (fenced run wall per cycle includes the fence: see below); ambiguous names (cannot be split with certainty; counted as calculation above): %s' % (SS, '; '.join('%s %.2f' % (k, 1e3*v/nc/len(R)) for k, v in amb.most_common(12))))
    print('- largest "others" names: ' + '; '.join('%s %.2f' % (k[2], 1e3*v/nc/len(R)) for k, v in sorted(((k, v) for k, v in names.items() if k[0] == 'others'), key=lambda x: -x[1])[:8]))
    # MPI per site
    M = collections.Counter(); kinds = collections.Counter(); W = []; msub = len(glob.glob(S + st + '/mt6/mt.*.txt')) == 4
    files = sorted(glob.glob(S + st + '/mt/mt.*.txt'))
    for f in files:
        for l in open(f):
            p = l.rstrip('\n').split('\t')
            if l[0] == '#' or len(p) < 3: continue
            q = p[0].split(); kind = p[1]; site = p[2]
            if kind.startswith('fence') or kind == 'deep_copy': kinds['non-MPI ' + ('fence' if kind.startswith('fence') else 'deep_copy')] += float(q[0]); continue
            mm = mod(site); M[(mm, kind, site)] += float(q[0])
    if msub:
        for f in sorted(glob.glob(S + st + '/mt6/mt.*.txt')):
            for l in open(f):
                p = l.rstrip('\n').split('\t')
                if l[0] == '#' or len(p) < 3: continue
                q = p[0].split(); kind = p[1]; site = p[2]
                if kind.startswith('fence') or kind == 'deep_copy': kinds['non-MPI ' + ('fence' if kind.startswith('fence') else 'deep_copy')] -= float(q[0]); continue
                M[(mod(site), kind, site)] -= float(q[0])
    nr = max(len(files), 1)
    pm = collections.Counter(); sites = collections.Counter()
    for (m, k, s2), v in M.items(): pm[m] += v; sites[(m, k, s2)] += v
    ncm = 20.0 if msub else 26.0
    print('MPI timer (' + ('26 minus 6 cycles, per 20 cycles' if msub else '26 cycles') + ', %d ranks, no fences): host time in MPI, ms per cycle (mean of the ranks), by module of the last kernel before the call: ' % nr + ', '.join('%s %.2f' % (m, 1e3*pm[m]/ncm/nr) for m in ('hydro', 'rays', 'rad-vet', 'others')))
    print('  - sites (module, call): ' + '; '.join('%s %s %s %.2f' % (k[0], k[1], k[2], 1e3*v/ncm/nr) for k, v in sites.most_common(14)))
    print('  - non-MPI host waits in the same timer (fences %.2f, deep copies %.2f ms per cycle, these overlap GPU work)' % (1e3*kinds.get('non-MPI fence', 0)/ncm/nr, 1e3*kinds.get('non-MPI deep_copy', 0)/ncm/nr))
    # nsys
    g = S + st + '/n/'
    if os.path.exists(g + 'kname_r0.txt'):
        f = None; rows = []
        for r in range(4):
            h = open(g + 'kname_r%d.txt' % r).readline(); L = open(g + 'gaps_r%d.txt' % r).read().split('\n')
            m = re.search(r'cycles (\S+) wall (\S+) ms/cyc\s+busy (\S+)', h); i = re.search(r'GPU idle (\S+) ms/cyc', L[0]); idl = dict(re.findall(r'(\S+?) +([0-9.]+)(?=  |$)', L[1].split(':', 1)[1])) if len(L) > 1 else {}
            t = open(g + 'kname_r%d.txt' % r).read(); z = re.search(r'^ *([0-9.]+) +\S+ +(\d+) +([0-9.]+).*CleanGdFaceSums#6', t, re.M); f = float(z[2]) if z else 1.0
            fx = float(z[1])/ (float(z[3])/1e3) if z else 1.0
            rows.append((float(m[2])/fx, float(m[3])/fx, float(i[1])/fx, {k: float(v)/fx for k, v in idl.items()}, fx, m[1]))
        print('nsys (kernels unfenced; per real cycle: values divided by %s = launches of CleanGdFaceSums#6 per printed cycle; printed cycles %s): mean of 4 ranks: wall %.1f, GPU busy %.1f, GPU idle %.1f ms per cycle' % ('%.2f' % np.mean([x[4] for x in rows]), rows[0][5], np.mean([x[0] for x in rows]), np.mean([x[1] for x in rows]), np.mean([x[2] for x in rows])))
        ks = sorted(set(k for x in rows for k in x[3])); print('  - idle by host activity (ms per cycle, mean of 4 ranks): ' + ', '.join('%s %.2f' % (k, np.mean([x[3].get(k, 0) for x in rows])) for k in ks))
    else: print('nsys: the run hung, no kname / gaps (kernel timer and MPI timer only)')
    print(open(S + st + '/plain.txt').read().split('\n')[0][:430])
