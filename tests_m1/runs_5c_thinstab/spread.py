"""spread.py <rundir>...: per bin dump, s = max over x1 of (max_x2 E - min_x2 E)/mean_x2 E
(also the same for F2/(c E)), and the per-step growth factor g fitted between the first
dump with s > 3 s0 floor-free and the first with s > 1e-2 (or the last dump)."""
import sys
import glob
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_thinstab/vis/python')
import bin_convert as bc  # noqa: E402


def series(d):
    out = []
    for fn in sorted(glob.glob(d + '/bin/*.bin')):
        f = bc.read_binary(fn)
        E = np.asarray(f['mb_data']['m1_e'][0], float)[0]      # [j, i]
        F2 = np.asarray(f['mb_data']['m1_f2'][0], float)[0]
        m = E.mean(axis=0)
        s = np.max((E.max(axis=0) - E.min(axis=0)) / m)
        sf = np.max(np.abs(F2)).item() / np.max(m)
        out.append((f['cycle'], s, sf))
    return out


if __name__ == '__main__':
    for d in sys.argv[1:]:
        s = series(d)
        cyc = np.array([x[0] for x in s])
        sp = np.array([x[1] for x in s])
        g = float('nan')
        ok = np.where((sp > 1e-5) & (sp < 1e-2))[0]
        if len(ok) >= 2:
            a, b = ok[0], ok[-1]
            if cyc[b] > cyc[a]:
                g = (sp[b] / sp[a]) ** (1.0 / (cyc[b] - cyc[a]))
        print('%-28s g/step=%7.3f  spread: %s' % (d.split('/')[-1], g,
              ' '.join('%d:%.1e' % (c, v) for c, v in zip(cyc, sp))))
