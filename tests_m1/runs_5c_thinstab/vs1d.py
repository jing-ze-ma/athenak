"""vs1d.py <2-D run> <1-D run>: max over x1 of |<E>_x2 - E_1D|/E_1D and the same for F1,
at matching dump cycles (every 100)."""
import sys
import glob
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_thinstab/vis/python')
import bin_convert as bc  # noqa: E402


def load(d):
    out = {}
    for fn in sorted(glob.glob(d + '/bin/*.bin')):
        f = bc.read_binary(fn)
        out[f['cycle']] = (np.asarray(f['mb_data']['m1_e'][0], float)[0].mean(axis=0),
                           np.asarray(f['mb_data']['m1_f1'][0], float)[0].mean(axis=0))
    return out


a, b = load(sys.argv[1]), load(sys.argv[2])
row = []
for c in sorted(set(a) & set(b)):
    if c % 100 == 0 or c == max(set(a) & set(b)):
        de = np.max(np.abs(a[c][0] - b[c][0]) / b[c][0])
        df = np.max(np.abs(a[c][1] - b[c][1])) / np.max(np.abs(b[c][1]))
        row.append('%d: dE %.1e dF %.1e' % (c, de, df))
print(sys.argv[1].split('/')[-1], 'vs', sys.argv[2].split('/')[-1], ' | '.join(row))
