import numpy as np, re, sys
def hst(f):
    names = re.findall(r'\[\d+\]=(\S+)', open(f).readlines()[1]); d = np.loadtxt(f, ndmin=2)
    return {n: d[:, i] for i, n in enumerate(names)}
for a in sys.argv[1:]:
    try: u = hst(a+'/bsg3d.user.hst')
    except Exception as e: print(a, 'no hst', e); continue
    t=u['time']; E=u['Etot']
    r = (np.diff(E) - (u['L_in']-u['L_top']-u['E_hflx']+u['P_esrc'])[1:]*np.diff(t))/E[0]
    print(f"{a:8s} n {len(t)-1} dt mean {np.mean(u['dt'][1:]):.4g} min {u['dt'][1:].min():.4g}; E resid/step median {np.median(r):.3e} net {r.sum():.3e}")
