# usage: python3 tsum.py <tag> [c0=20] [c1=110]
# ms/cycle over cycles c0..c1 per run, then mean per arm; Picard/solve, inner/solve, NONCONV
import re, os, sys, collections
W = '/viper/ptmp2/jinma/launch_0923/runs/' + sys.argv[1]
c0 = int(sys.argv[2]) if len(sys.argv) > 2 else 20
c1 = int(sys.argv[3]) if len(sys.argv) > 3 else 110


def g(pat, L):
    x = [re.search(pat, l) for l in L]
    x = [m for m in x if m]
    return float(x[-1].group(1)) if x else float('nan')


arms = collections.OrderedDict()
for a in sorted(os.listdir(W)):
    p = f'{W}/{a}/run.log'
    if not os.path.exists(p):
        continue
    name, rep = a.rsplit('_', 1)
    L = open(p).read().splitlines()
    try:
        ec = {int(re.search(r'cycle=(\d+)', l).group(1)):
              float(re.search(r'elapsed=(\S+)', l).group(1))
              for l in L if l.startswith('elapsed=')}
        ms = 1e3*(ec[c1]-ec[c0])/(c1-c0)
    except Exception:
        ms = float('nan')
    n = g(r'solves=(\S+) Picard', L)
    pic = g(r'Picard iterations mean=(\S+)', L)
    its = g(r'inner iterations mean=\S+ max=\S+ total=(\S+)', L)/n if n == n else float('nan')
    nc = g(r'NON-CONVERGED=(\S+)', L)
    arms.setdefault(name, []).append((rep, ms, pic, its, nc))
print(f"{'arm':14s} {'ms/cyc runs':28s} {'mean':>6s} {'Picard':>7s} {'inner':>7s} {'NC':>4s}")
for k, v in arms.items():
    v.sort()
    ms = [x[1] for x in v]
    m = sum(ms)/len(ms)
    print(f"{k:14s} {', '.join('%.1f' % x for x in ms):28s} {m:6.1f} {v[0][2]:7.3f} "
          f"{v[0][3]:7.2f} {v[0][4]:4.0f}")
