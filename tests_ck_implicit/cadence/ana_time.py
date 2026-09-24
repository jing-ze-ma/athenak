# ms/cycle per arm from elapsed= diag lines (skip first 20 cycles), passes from verbose arms
import glob, re, os, statistics as st
import sys
P = sys.argv[1] if len(sys.argv) > 1 else '/viper/ptmp2/jinma/ckprof_0923/time'
res = {}
for d in sorted(glob.glob(P + '/*_r*')):
    arm = os.path.basename(d)
    L = [l for l in open(d + '/run.log') if l.startswith('elapsed=')]
    v = [(float(re.search(r'elapsed=(\S+)', l).group(1)), int(re.search(r'cycle=(\d+)', l).group(1)),
          float(re.search(r'dt=(\S+)', l).group(1))) for l in L]
    if len(v) < 3:
        print(arm, 'FAILED/short', len(v)); continue
    e0, c0, _ = v[1]; e1, c1, _ = v[-1]
    ms = 1e3 * (e1 - e0) / (c1 - c0)
    res.setdefault(arm.rsplit('_', 1)[0], []).append(ms)
    pl = [int(m) for m in re.findall(r'passes=(\d+) rank=0', open(d + '/run.log').read())]
    mem = ''
    if os.path.exists(d + '/mem.txt'):
        mem = ' '.join(re.findall(r'VRAM Total Used Memory \(B\): (\d+)', open(d + '/mem.txt').read()))
    print(f'{arm:10s} {ms:8.2f} ms/cyc  cycles {c1-c0}  dt {v[-1][2]:.3f}'
          + (f'  calls {len(pl)} mean passes {st.mean(pl):.2f} max {max(pl)}' if pl else '')
          + (f'  vramB {mem}' if mem else ''))
h = st.median(res['h']) if 'h' in res else None
for a, x in res.items():
    print(f'{a:6s} median {st.median(x):8.2f} ms  n={len(x)}' + (f'  x hydro {st.median(x)/h:.2f}' if h else ''))
