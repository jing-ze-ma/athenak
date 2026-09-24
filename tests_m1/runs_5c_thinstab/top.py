"""top.py <rundir>...: per bin dump, the top-cell E (x2 mean), top-cell f = F1/(cE)
(c = 100) and the max over cells of |E - E_prev|/E between dumps."""
import sys
import glob
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_thinstab/vis/python')
import bin_convert as bc  # noqa: E402

for d in sys.argv[1:]:
    prev, out = None, []
    for fn in sorted(glob.glob(d + '/bin/*.bin')):
        f = bc.read_binary(fn)
        E = np.asarray(f['mb_data']['m1_e'][0], float)[0]
        F1 = np.asarray(f['mb_data']['m1_f1'][0], float)[0]
        ch = 0.0 if prev is None else np.max(np.abs(E - prev) / E)
        out.append('%d:E=%.4g f=%.3f d=%.1e' % (f['cycle'], E[:, -1].mean(),
                   (F1[:, -1] / (100 * E[:, -1])).mean(), ch))
        prev = E
    print(d.split('/')[-1], ' '.join(out[::2] + [out[-1]]))
