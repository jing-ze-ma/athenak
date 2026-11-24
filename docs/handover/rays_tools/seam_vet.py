# seam.py R lo hi f.m1.bin...: hole (|F|>cE) fraction vs distance to the cubed-sphere panel edge, d = 1 - max(|x2|,|x3|), and per panel
import sys, numpy as np
sys.path.insert(0, sys.argv[0].rsplit('/',1)[0]); import bin_convert as bc
R, lo, hi = map(float, sys.argv[1:4]); c = 2.99792458e10
for f in sys.argv[4:]:
    m = bc.read_binary(f); nb = m['n_mbs']
    bins = np.array([0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.01]); H = np.zeros(len(bins)-1); N = np.zeros(len(bins)-1); pan = {}
    for b in range(nb):
        x = np.array(m['mb_x1v'][b]); s = (x >= lo*R) & (x < hi*R)
        if not s.any(): continue
        E = np.array(m['mb_data']['vet_e'][b], dtype=np.float64)[:, :, s]
        Fm = np.sqrt(sum(np.array(m['mb_data'][k][b], dtype=np.float64)[:, :, s]**2 for k in ('vet_f1', 'vet_f2', 'vet_f3')))
        h = (Fm > c*E).mean(axis=2)
        x2 = np.array(m['mb_x2v'][b]); x3 = np.array(m['mb_x3v'][b])
        d = 1 - np.maximum(np.abs(x2)[None, :], np.abs(x3)[:, None])
        k = np.digitize(d, bins) - 1
        for i in range(len(N)): H[i] += h[k == i].sum(); N[i] += (k == i).sum()
        lg = tuple(m['mb_logical'][b]) if 'mb_logical' in m else b
        pan.setdefault(b % 6, []).append(h.mean())
    print(f.split('/')[-1], 'r/R %.1f-%.1f' % (lo, hi), ' '.join('d %.2f-%.2f: %.3f' % (bins[i], bins[i+1], H[i]/max(N[i], 1)) for i in range(len(N))))
