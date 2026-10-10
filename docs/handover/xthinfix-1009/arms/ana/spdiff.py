"""xthinfix-1009: sp wedge difference of two runs (mode all vs beam) by radius.
usage: REPO=<athenak> python3 spdiff.py <run A> <run B> [R_ph cm, default AG Car A 388.3 Rsun]
For the last bin dump of every output (same file name in both runs) and every variable: max over
(theta, phi) of |a - b| / max(|b|) per radial shell, reported in r/R_ph bands."""
import glob
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.environ['REPO'], 'vis/python'))
import bin_convert as bc  # noqa: E402

A, B = sys.argv[1], sys.argv[2]
RPH = float(sys.argv[3]) if len(sys.argv) > 3 else 388.3*6.957e10
BANDS = [(0, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 1.05), (1.05, 1.3), (1.3, 2.0), (2.0, 10.0)]
ids = {}
for f in sorted(glob.glob(B + '/bin/*.bin')):
    stem = os.path.basename(f).rsplit('.', 2)[0]
    ids[stem] = f            # last dump per output id
print('bands r/R_ph:', ' '.join('%g-%g' % b for b in BANDS))
for stem, fb in sorted(ids.items()):
    fa = os.path.join(A, 'bin', os.path.basename(fb))
    if not os.path.exists(fa):
        print(stem, 'missing in', A)
        continue
    a = bc.read_binary_as_athdf(fa)
    b = bc.read_binary_as_athdf(fb)
    r = np.asarray(b['x1v'])/RPH
    keys = [k for k in b.keys() if isinstance(b[k], np.ndarray) and b[k].ndim == 3]
    for k in keys:
        x, y = np.asarray(a[k], dtype=float), np.asarray(b[k], dtype=float)
        sc = np.maximum(np.abs(y).max(axis=(0, 1)), 1e-300)
        d = (np.abs(x - y).max(axis=(0, 1)))/sc
        out = []
        for lo, hi in BANDS:
            s = (r >= lo) & (r < hi)
            out.append('%.1e' % d[s].max() if s.any() else '   -   ')
        print('%-28s %-10s %s' % (os.path.basename(fb), k, ' '.join(out)))
