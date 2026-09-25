# usage: python3 tsum.py <tag> [c0=10] [c1=50]: ms/cycle per run and per arm; NC; inner its
import re, os, sys, collections
W = '/viper/ptmp2/jinma/fast3/runs/' + sys.argv[1]
c0 = int(sys.argv[2]) if len(sys.argv) > 2 else 10
c1 = int(sys.argv[3]) if len(sys.argv) > 3 else 50
arms = collections.OrderedDict()
for a in sorted(os.listdir(W)):
    p = f'{W}/{a}/run.log'
    if not os.path.exists(p) or a.startswith('p'):
        continue
    name, rep = a.rsplit('_', 1)
    t = open(p).read()
    ec = {int(m.group(2)): float(m.group(1))
          for m in re.finditer(r'elapsed=(\S+) cycle=(\d+)', t)}
    ms = 1e3*(ec[c1]-ec[c0])/(c1-c0) if c0 in ec and c1 in ec else float('nan')
    nc = re.findall(r'NON-CONVERGED=(\S+)', t)
    it = re.findall(r'inner iterations mean=(\S+)', t)
    arms.setdefault(name, []).append((ms, nc[-1] if nc else '?', it[-1] if it else '?'))
for n, v in arms.items():
    ms = [x[0] for x in v]
    print(f'{n:10s} mean {sum(ms)/len(ms):7.2f}  reps ' + ' '.join(f'{x:6.2f}' for x in ms)
          + '  NC ' + ','.join(x[1] for x in v) + '  inner ' + ','.join(x[2] for x in v))
