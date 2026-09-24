#!/usr/bin/env python3
"""fbcmp.py <gpuf dir>: gas-Newton fallbacks per arm and distance (gates/cmp.py) to the
tight reference (_T) and to time2_lin_tol_fac = 1 (_L1) of the same closure."""
import os, re, subprocess, sys
G = sys.argv[1]
CMP = '/viper/ptmp2/jinma/wt_m1def2/tests_m1/gates/cmp.py'
for n in sorted(os.listdir(G)):
    lg = os.path.join(G, n, 'run.log')
    if not os.path.exists(lg):
        continue
    t = open(lg).read()
    fb = re.findall(r'newton fallbacks=(\S+) \((\S+) per', t)
    nc = re.findall(r'NON-CONVERGED\S*\s*=?\s*(\S+)', t)
    pic = re.findall(r'Picard\S* mean\S*=?\s*(\S+)', t)
    print(f"{n:7s} newton_fallbacks={fb[-1] if fb else '?'} NC={nc[-1] if nc else '?'}")
    for ref in ('T', 'L1'):
        r = n.split('_')[0] + '_' + ref
        if r == n or not os.path.isdir(os.path.join(G, r)):
            continue
        o = subprocess.run(['python3', CMP, os.path.join(G, r), os.path.join(G, n)],
                           capture_output=True, text=True).stdout.strip()
        print(f"    vs {r}: {o}")
