# ms/cycle (cycles 20-120), Picard mean, NON-CONVERGED, fallbacks, rad_w builds per GPU run
import re, os
import sys
W = '/viper/ptmp2/jinma/h2fast_0924/' + (sys.argv[1] if len(sys.argv) > 1 else 'gpu')
def g(pat, L):
    x = [re.search(pat, l) for l in L]; x = [m for m in x if m]
    return x[-1].group(1) if x else 'nan'
for a in sorted(os.listdir(W)):
    if not os.path.isfile(f'{W}/{a}/run.log') or a.startswith('g0_'): continue
    L = open(f'{W}/{a}/run.log').read().splitlines()
    try:
        ec = {int(re.search(r'cycle=(\d+)', l).group(1)): float(re.search(r'elapsed=(\S+)', l).group(1))
              for l in L if l.startswith('elapsed=')}
        ms = 1e3*(ec[120]-ec[20])/100.0
    except Exception:
        ms = float('nan')
    fb = [re.findall(r'fallbacks=(\S+)', l) for l in L]
    fb = sorted(set(v for f in fb for v in f))
    pic = g(r'Picard iterations mean=(\S+)', L); nc = g(r'NON-CONVERGED=(\S+)', L)
    rw = g(r'### rad_w builds: (.*)', L)
    print(f"{a:14s} ms/cyc={ms:7.2f} Picard={pic} NC={nc} fallback-values={fb} radw=[{rw}]")
