import glob, os, sys
import numpy as np
sys.path.insert(0, os.environ['ATHENA_VIS'])
import bin_convert
d = sys.argv[1]
out = []
for f in sorted(glob.glob(os.path.join(d, 'bin', '*.bin'))):
    b = bin_convert.read_binary(f)
    e = np.asarray(b['mb_data']['m1_e'][0])
    g = b['mb_geometry'][0]
    n3, n2, n1 = e.shape
    rf = np.linspace(g[0], g[1], n1 + 1); tf = np.linspace(g[2], g[3], n2 + 1)
    pf = np.linspace(g[4], g[5], n3 + 1)
    v = ((pf[1:] - pf[:-1])[:, None, None]*(np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :, None]
         *((rf[1:]**3 - rf[:-1]**3)/3.0)[None, None, :])
    out.append((b['cycle'], b['time'], np.sum(e*v)))
e0 = out[0][2]
prev = e0
for c, t, s in out:
    print('%5d %.5f %.17e rel %+.3e step %+.3e' % (c, t, s, s/e0 - 1, s/prev - 1))
    prev = s
