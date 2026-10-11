# analyse.py <runs/JOBID>: TASK-2026-10-11-deltaai-icgen report numbers
import sys, os, re, glob, hashlib, warnings
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, '/u/jma20/ATHENAK/builds/wt/hes_gpu/vis/python')
import bin_convert as bc
S = sys.argv[1]; RPH = 2.7014e13; ROT = 10
ARMS = [a for a in ('g0', 'gha', 'gha2', 'ghafr', 'gra', 'gh') if os.path.exists(f'{S}/{a}/out.txt')]
gf = '/work/nvme/bivj/jma20/icgen_1011/agcar3d.x1grid.txt'
x1v = np.loadtxt(gf, max_rows=486)[3:483, 2]; r = x1v / RPH

def hst(a):
    f = f'{S}/{a}/agcar3d.user.hst'; hdr = [l for l in open(f) if l.startswith('#')][-1]
    cols = {m.group(2): int(m.group(1)) - 1 for m in re.finditer(r'\[(\d+)\]=(\S+)', hdr)}
    return np.loadtxt(f), cols
O = {a: open(f'{S}/{a}/out.txt', errors='replace').read() for a in ARMS}
print('== per run')
for a in ARMS:
    o = O[a]; el = [(int(m.group(2)), float(m.group(1))) for m in re.finditer(r'elapsed=(\S+) cycle=(\d+)', o)]
    last = el[-1][0]; st = [t1 - t0 for (c0, t0), (c1, t1) in zip(el[:-1], el[1:]) if c1 == c0 + 1 and c1 % ROT and c1 != last]
    ro = [t1 - t0 for (c0, t0), (c1, t1) in zip(el[:-1], el[1:]) if c1 == c0 + 1 and c1 % ROT == 0 and c1 != last]
    p = re.findall(r'solves=(\S+) Picard iterations mean=(\S+) max=(\S+) NON-CONVERGED=(\S+)', o)[-1]
    print(f'{a:6s} cycles {last} FATAL {o.count("FATAL")} NONCONV {float(p[3]):g} Picard {float(p[1]):.2f}/{float(p[2]):.0f} '
          f'(solves {float(p[0]):.0f})  s/cycle median {np.median(st):.3f} (rotation {np.mean(ro):.3f}, all {np.mean(st+ro):.3f})  setup elapsed@c{el[0][0]} {el[0][1]:.1f} s')
    for l in o.splitlines():
        if re.search(r'he_ic_rad = |ic_rad run_c0|ic_rad: done', l): print('    ' + l.strip()[:220])
H = {a: hst(a) for a in ARMS}
print('\n== Picard (hst) mean per 20-cycle window')
print('arm    ' + ' '.join(f'{lo:>3d}-{lo+19:<3d}' for lo in range(1, 121, 20)))
for a in ARMS:
    h, c = H[a]; P = h[:, c['Picard']]
    print(f'{a:6s} ' + ' '.join(f'{P[lo:lo+20].mean():7.2f}' for lo in range(1, 121, 20) if lo + 19 < len(P)))
print('\n== L_top / L_top(g0, cycle 0)')
L0 = H['g0'][0][0, H['g0'][1]['L_top']]
cyc = list(range(0, 11)) + [20, 60, 120]
print('cycle ' + ' '.join(f'{a:>9s}' for a in ARMS))
for c in cyc:
    print(f'{c:5d} ' + ' '.join(f'{H[a][0][c, H[a][1]["L_top"]] / L0:9.5f}' if c < len(H[a][0]) else '        -' for a in ARMS))

def dumps(a):
    out = {}
    for f in glob.glob(f'{S}/{a}/bin/agcar3d.hydro_w.*.bin'):
        b = bc.read_binary(f); out.setdefault(b['cycle'], b)
    return out
def smean(b, v):
    num = 0.; den = 0.
    for m in range(b['n_mbs']):
        g = b['mb_geometry'][m]; x = np.array(b['mb_data'][v][m]); n2 = x.shape[1]
        w = np.sin(g[2] + (np.arange(n2) + .5) * (g[3] - g[2]) / n2)[None, :, None]
        num = num + (x * w).sum(axis=(0, 1)); den = den + w.sum() * x.shape[0]
    return num / den
def cmax(b, v, k):
    return max(np.abs(np.array(b['mb_data'][v][m])[:, :, k]).max() for m in range(b['n_mbs']))
print('\n== photosphere 0.97-1.028 R_ph: max|v_r| km/s, max|<rho>/<rho>(c0)-1|; shell-mean v_r at 2.0 R_ph km/s')
k = (r >= 0.97) & (r <= 1.028); i2 = np.argmin(np.abs(r - 2.0))
for a in ARMS:
    D = dumps(a); r0 = smean(D[0], 'dens'); s = []
    for c in (10, 30, 60, 120):
        if c not in D: s.append(f'c{c}: -'); continue
        rm = smean(D[c], 'dens'); vr = smean(D[c], 'velx')
        s.append(f'c{c}: {cmax(D[c], "velx", k) / 1e5:.3f} / {np.max(np.abs(rm[k] / r0[k] - 1)):.4f} / {vr[i2] / 1e5:+.3f}')
    print(f'{a:6s} ' + ' | '.join(s))
print(f'  (shells in band: {k.sum()}, r/R_ph {r[k][0]:.4f}-{r[k][-1]:.4f}; 2.0 R_ph shell at {r[i2]:.4f})')

def md5(f): return hashlib.md5(open(f, 'rb').read()).hexdigest()
if 'gha' in ARMS and 'gha2' in ARMS:
    a, b = f'{S}/gha', f'{S}/gha2'
    fs = sorted(os.path.basename(x) for x in glob.glob(a + '/bin/*.bin'))
    eq = 0
    for f in fs:
        A, B = bc.read_binary(f'{a}/bin/{f}'), bc.read_binary(f'{b}/bin/{f}')
        eq += all(np.array_equal(np.asarray(A['mb_data'][v]), np.asarray(B['mb_data'][v])) for v in A['var_names'])
    hs = {f: md5(f'{a}/{f}') == md5(f'{b}/{f}') for f in ('agcar3d.user.hst', 'agcar3d.hydro.hst')}
    print(f'\n== gha vs gha2: bin arrays identical {eq}/{len(fs)}, bin files md5-identical '
          f'{sum(md5(f"{a}/bin/{f}") == md5(f"{b}/bin/{f}") for f in fs)}/{len(fs)}; hst {hs}')
