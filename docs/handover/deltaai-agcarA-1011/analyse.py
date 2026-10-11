# analyse.py <runs/JOBID>   (TASK-2026-10-11-deltaai-agcarA-face-icrad)
# Per run: rc/FATAL/NON-CONVERGED (from out.txt), Picard mean/max (user.hst col 22, cycles >= 1),
# s/cycle from the per-cycle 'elapsed=' lines (excluding rotation cycles and the last cycle).
# Run 1: max rel diff face vs cell of sin(theta)-weighted shell means (E, F1, dens, velx) in r/R_ph bands at cycle 50,
#        L_top; bitwise check face vs face2 (bin + hst).
# Run 2: L_top per cycle; dt vs dt/2 L_top difference over cycles 1-30 (dt/2 matched at cycle 2n); IC step cost.
import sys, os, re, glob, hashlib
import numpy as np
sys.path.insert(0, '/u/jma20/ATHENAK/builds/wt/hes_gpu/vis/python')
import bin_convert as bc

S = sys.argv[1]
RPH = 388.3 * 6.957e10
BANDS = [('<0.8', 0.0, 0.8), ('0.95-1.05', 0.95, 1.05), ('>1.05', 1.05, 1e9)]

def hst(d):
    f = os.path.join(d, 'agcar3d.user.hst')
    if not os.path.exists(f): return None
    return np.loadtxt(f)

def runinfo(n):
    d = os.path.join(S, n); o = open(os.path.join(d, 'out.txt'), errors='replace').read()
    el = [(int(m.group(2)), float(m.group(1))) for m in re.finditer(r'elapsed=(\S+) cycle=(\d+)', o)]
    nc = re.findall(r'NON-CONVERGED=(\S+)', o)
    rot = 20 if n.endswith('_h') else 10
    sc = []
    for (c0, t0), (c1, t1) in zip(el[:-1], el[1:]):
        if c1 == c0 + 1 and c1 % rot != 0 and c1 != el[-1][0]: sc.append(t1 - t0)
    h = hst(d); pic = h[1:, 21] if h is not None and len(h) > 1 else np.array([np.nan])
    ic = el[1][1] - el[0][1] if len(el) > 1 else np.nan   # cycle 1 (includes nothing of setup)
    return dict(fatal=o.count('FATAL'), nonconv=nc[-1] if nc else 'n/a', ncyc=el[-1][0] if el else -1,
                pmean=np.nanmean(pic), pmax=np.nanmax(pic), spc=np.median(sc) if sc else np.nan,
                spcmean=np.mean(sc) if sc else np.nan, t0=el[0][1] if el else np.nan, h=h)

def shellmeans(d, cyc):
    gf = os.path.join(d, 'agcar3d.x1grid.txt'); gf = gf if os.path.exists(gf) else '/work/nvme/bivj/jma20/agcarA_1011/agcar3d.x1grid.txt'
    x1v = np.loadtxt(gf, max_rows=486)[3:483, 2]
    out = {}
    for kind, vs in (('m1', ['m1_e', 'm1_f1']), ('hydro_w', ['dens', 'velx'])):
        f = glob.glob(os.path.join(d, 'bin', f'agcar3d.{kind}.*.bin'))
        f = [x for x in f if bc.read_binary(x)['cycle'] == cyc][0]
        b = bc.read_binary(f)
        for v in vs:
            num = 0.; den = 0.
            for m in range(b['n_mbs']):
                g = b['mb_geometry'][m]; a = np.array(b['mb_data'][v][m])   # (nx3, nx2, nx1)
                n2 = a.shape[1]; th = g[2] + (np.arange(n2) + 0.5) * (g[3] - g[2]) / n2
                w = np.sin(th)[None, :, None]
                num = num + (a * w).sum(axis=(0, 1)); den = den + w.sum() * a.shape[0]
            out[v] = num / den
    return x1v, out

def md5(f): return hashlib.md5(open(f, 'rb').read()).hexdigest()

R = {}
for n in sorted(os.listdir(S)):
    if os.path.exists(os.path.join(S, n, 'out.txt')): R[n] = runinfo(n)
print('run  ncyc  fatal  NON-CONV  Picard_mean  Picard_max  s/cycle(median, mean)  elapsed@cycle0')
for n, r in R.items():
    print(f"{n:13s} {r['ncyc']:4d} {r['fatal']:3d} {r['nonconv']:>14s} {r['pmean']:8.2f} {r['pmax']:6.0f}"
          f"  {r['spc']:.3f} {r['spcmean']:.3f}  {r['t0']:.1f}")

if all(k in R for k in ('A_cell', 'A_face')):
    print('\nRun 1: face vs cell at cycle 50, max |rel diff| of shell means per band')
    xc, mc = shellmeans(os.path.join(S, 'A_cell'), 50); xf, mf = shellmeans(os.path.join(S, 'A_face'), 50)
    r = xc / RPH
    for v in mc:
        s = []
        for bn, lo, hi in BANDS:
            k = (r >= lo) & (r < hi)
            rel = np.abs(mf[v][k] - mc[v][k]) / np.maximum(np.abs(mc[v][k]), 1e-300)
            s.append(f"{bn}: {rel.max():.3e} (r/Rph {r[k][rel.argmax()]:.3f})")
        print(f"  {v:6s} " + ' | '.join(s))
    hc, hf = R['A_cell']['h'], R['A_face']['h']; m = min(len(hc), len(hf))
    rl = np.abs(hf[:m, 5] - hc[:m, 5]) / np.abs(hc[:m, 5])
    print(f"  L_top: cell {hc[m-1,5]:.6e} face {hf[m-1,5]:.6e} at cycle {m-1}; max rel diff over cycles {rl.max():.3e} (cycle {rl.argmax()})")
    if 'A_face2' in R:
        d1, d2 = os.path.join(S, 'A_face'), os.path.join(S, 'A_face2')
        fs = sorted(os.path.basename(x) for x in glob.glob(d1 + '/bin/*.bin'))
        same = [md5(d1 + '/bin/' + f) == md5(d2 + '/bin/' + f) for f in fs if os.path.exists(d2 + '/bin/' + f)]
        hs = [md5(d1 + '/' + f) == md5(d2 + '/' + f) for f in ('agcar3d.user.hst', 'agcar3d.hydro.hst')]
        print(f"  face repeat bitwise: bin {sum(same)}/{len(fs)} identical, hst {hs}")

print('\nRun 2: L_top per cycle (erg/s)')
B = [n for n in ('B_def', 'B_fs', 'B_relax', 'B_fsrelax') if n in R]
if B:
    m = min(len(R[n]['h']) for n in B)
    print('cycle ' + ' '.join(f'{n:>14s}' for n in B))
    for c in range(m):
        print(f'{c:5d} ' + ' '.join(f"{R[n]['h'][c,5]:14.6e}" for n in B))
    print('\nRun 2: max |L_top(dt) - L_top(dt/2)| / L_top(dt/2) over cycles 1-30 (dt/2 at cycle 2c), and rms')
    for n in B:
        if n + '_h' not in R: continue
        h1, h2 = R[n]['h'], R[n + '_h']['h']
        c = np.arange(1, 31); c = c[(c < len(h1)) & (2 * c < len(h2))]
        assert np.allclose(h1[c, 0], h2[2 * c, 0], rtol=1e-12), 'time mismatch'
        d = np.abs(h1[c, 5] - h2[2 * c, 5]) / np.abs(h2[2 * c, 5])
        print(f"  {n:10s} max {d.max():.3e} (cycle {c[d.argmax()]})  rms {np.sqrt((d**2).mean()):.3e}  n={len(c)}")
    print('\nIC step cost: elapsed at cycle 0 (setup incl. IC radiation start), vs B_def')
    for n in B + [x + '_h' for x in B]:
        if n in R: print(f"  {n:12s} {R[n]['t0']:.2f} s  (+{R[n]['t0'] - R['B_def' if not n.endswith('_h') else 'B_def_h']['t0']:.2f})")
