# analyse.py <runs/JOBID>: TASK-2026-10-11-deltaai-prodcand-agcarA report numbers
import sys, os, re, glob, hashlib, warnings
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, '/u/jma20/ATHENAK/builds/wt/hes_gpu/vis/python')
import bin_convert as bc
S = sys.argv[1]; RPH = 388.3 * 6.957e10; ROT = 10; DT = 1239.996
ARMS = [a for a in ('cand1', 'cand2', 'ref') if os.path.exists(f'{S}/{a}/out.txt')]
x1v = np.loadtxt('/work/nvme/bivj/jma20/prodcand_1011/agcar3d.x1grid.txt', max_rows=486)[3:483, 2]; r = x1v / RPH

def hst(a):
    f = f'{S}/{a}/agcar3d.user.hst'; hdr = [l for l in open(f) if l.startswith('#')][-1]
    return np.loadtxt(f), {m.group(2): int(m.group(1)) - 1 for m in re.finditer(r'\[(\d+)\]=(\S+)', hdr)}
print('== per run (s/cycle: steady = median excl. rotation cycles and the last; rotation = mean of rotation cycles; average = mean of all)')
for a in ARMS:
    o = open(f'{S}/{a}/out.txt', errors='replace').read()
    el = [(int(m.group(2)), float(m.group(1))) for m in re.finditer(r'elapsed=(\S+) cycle=(\d+)', o)]
    last = el[-1][0]; d = [(c1, t1 - t0) for (c0, t0), (c1, t1) in zip(el[:-1], el[1:]) if c1 == c0 + 1]
    st = [x for c, x in d if c % ROT and c != last]; ro = [x for c, x in d if c % ROT == 0 and c != last]; al = [x for c, x in d]
    p = re.findall(r'solves=(\S+) Picard iterations mean=(\S+) max=(\S+) NON-CONVERGED=(\S+)', o)[-1]
    print(f'{a:6s} cycles {last} FATAL {o.count("FATAL")} NONCONV {float(p[3]):g} Picard {float(p[1]):.2f}/{float(p[2]):.0f} (solves {float(p[0]):.0f})'
          f'  steady {np.median(st):.3f} rotation {np.mean(ro):.3f} average {np.mean(al):.3f} s/cycle'
          f'  cost s per 1e4 s sim: {1e4*np.median(st)/DT:.3f} / {1e4*np.mean(al)/DT:.3f} (steady / average)  setup elapsed@c0 {el[0][1]:.1f} s')
    for l in o.splitlines():
        if re.search(r'ic_rad: done|ic_rad = ', l): print('    ' + l.strip()[:200])
H = {a: hst(a) for a in ARMS}
if 'cand1' in H and 'ref' in H:
    print('\n== L_top (erg/s) per cycle, cand1 vs ref')
    hc, cc = H['cand1']; hr, cr = H['ref']
    for c in range(0, min(len(hc), len(hr), 61)):
        print(f'{c:3d} {hc[c, cc["L_top"]]:.6e} {hr[c, cr["L_top"]]:.6e}')
    l = hc[:61, cc['L_top']]; print(f'cand1 L_top over cycles 0-60: min/L0 {l.min()/l[0]:.6f} max/L0 {l.max()/l[0]:.6f}')

def last(a, kind):
    bs = [bc.read_binary(f) for f in glob.glob(f'{S}/{a}/bin/agcar3d.{kind}.*.bin')]
    return max(bs, key=lambda b: b['cycle'])
def cells(b, v): return np.stack([np.array(b['mb_data'][v][m]) for m in range(b['n_mbs'])])   # (mb, k, j, i)
def smean(b, v):
    num = 0.; den = 0.
    for m in range(b['n_mbs']):
        g = b['mb_geometry'][m]; x = np.array(b['mb_data'][v][m]); n2 = x.shape[1]
        w = np.sin(g[2] + (np.arange(n2) + .5) * (g[3] - g[2]) / n2)[None, :, None]
        num = num + (x * w).sum(axis=(0, 1)); den = den + w.sum() * x.shape[0]
    return num / den

if 'cand1' in ARMS and 'cand2' in ARMS:
    fs = sorted(os.path.basename(x) for x in glob.glob(f'{S}/cand1/bin/*.bin')); eq = 0
    for f in fs:
        A, B = bc.read_binary(f'{S}/cand1/bin/{f}'), bc.read_binary(f'{S}/cand2/bin/{f}')
        eq += all(np.array_equal(np.asarray(A['mb_data'][v]), np.asarray(B['mb_data'][v])) for v in A['var_names'])
    md = lambda f: hashlib.md5(open(f, 'rb').read()).hexdigest()
    hs = {f: md(f'{S}/cand1/{f}') == md(f'{S}/cand2/{f}') for f in ('agcar3d.user.hst', 'agcar3d.hydro.hst')}
    print(f'\n== cand2 vs cand1: bin arrays identical {eq}/{len(fs)}; hst {hs}')

print('\n== E/J at the last dump (cycle 60): median (p5-p95) over cells in the band')
BE = [('1.0-1.05', 1.0, 1.05), ('1.05-1.5', 1.05, 1.5), ('1.5-2', 1.5, 2.0), ('2-3', 2.0, 3.01)]
for a in ARMS:
    m1 = last(a, 'm1'); fsb = last(a, 'm1_fs'); E = cells(m1, 'm1_e')
    sw = re.search(r'^implicit_thin_switch *= *(\S+)', open(f'{S}/{a}/deck.athinput').read(), re.M) if os.path.exists(f'{S}/{a}/deck.athinput') else None
    pairs = [('E/J_fs', cells(fsb, 'm1_jfs'))] + ([('E/J_tw', cells(last(a, 'm1_thsw'), 'm1_tj'))] if sw and sw.group(1) != 'none' else [])
    if len(pairs) == 1: print(f'  {a:6s} E/J_tw: n/a (implicit_thin_switch {sw.group(1) if sw else "?"}: m1_tj is not a twin J)')
    for nm, J in pairs:
        s = []
        for bn, lo, hi in BE:
            k = (r >= lo) & (r < hi); q = (E[..., k] / J[..., k]).ravel(); q = q[np.isfinite(q)]
            s.append(f'{bn}: {np.median(q):.4f} ({np.percentile(q,5):.4f}-{np.percentile(q,95):.4f})' if q.size else f'{bn}: n/a')
        print(f'  {a:6s} {nm} (cycle {m1["cycle"]}) ' + ' | '.join(s))
if 'cand1' in ARMS and 'ref' in ARMS:
    print('\n== cand1 vs ref at cycle 60: max over shells in the band of |<q>_cand1 - <q>_ref| / |<q>_ref| (shell means)')
    BB = [('<0.8', 0, 0.8), ('0.95-1.05', 0.95, 1.05), ('>1.05', 1.05, 1e9)]
    for kind, vs in (('hydro_w', ('dens', 'velx')), ('m1', ('m1_e', 'm1_f1'))):
        A, R = last('cand1', kind), last('ref', kind)
        for v in vs:
            ma, mr = smean(A, v), smean(R, v); rel = np.abs(ma - mr) / np.abs(mr)
            print(f'  {v:5s} ' + ' | '.join(f'{bn}: {rel[(r>=lo)&(r<hi)].max():.3e} (r/Rph {r[(r>=lo)&(r<hi)][rel[(r>=lo)&(r<hi)].argmax()]:.3f})' for bn, lo, hi in BB))
