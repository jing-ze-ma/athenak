# usage: python3 summarize.py [prefix]  -> ms/cycle (cycles 20-120), Picard/step, inner/step
import re, os, sys
W = '/viper/ptmp2/jinma/accel_0923/runs'
pre = sys.argv[1] if len(sys.argv) > 1 else ''
def g(pat, L):
    x = [re.search(pat, l) for l in L]; x = [m for m in x if m]
    return float(x[-1].group(1)) if x else float('nan')
print(f"{'arm':22s} {'ms/cyc':>8s} {'Picard':>7s} {'inner/step':>10s} {'NONCONV':>7s}")
for a in sorted(os.listdir(W)):
    if not a.startswith(pre) or a.startswith('LOG'): continue
    try:
        L = open(f'{W}/{a}/run.log').read().splitlines()
        ec = {int(re.search(r'cycle=(\d+)', l).group(1)): float(re.search(r'elapsed=(\S+)', l).group(1))
              for l in L if l.startswith('elapsed=')}
        ms = 1e3*(ec[120]-ec[20])/100.0
        n = g(r'solves=(\S+) Picard', L)
        pic = g(r'Picard iterations mean=(\S+)', L)
        its = g(r'inner iterations mean=\S+ max=\S+ total=(\S+)', L)/n
        nc = g(r'NON-CONVERGED=(\S+)', L)
        print(f"{a:22s} {ms:8.1f} {pic:7.3f} {its:10.2f} {nc:7.0f}")
    except Exception as e: print(a, 'ERR', e)
