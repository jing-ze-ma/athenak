# compare spec_<tag>.npz and ana_<tag>.npz: new dir vs reference dir
import sys, numpy as np
new, ref = sys.argv[1], sys.argv[2]
for tag in sys.argv[3:]:
    a, b = np.load('%s/spec_%s.npz' % (new, tag)), np.load('%s/spec_%s.npz' % (ref, tag))
    s, r = a['spec'].astype(float), b['spec'].astype(float)
    d = np.abs(s - r); nz = np.abs(r) > 0
    rel = (d[nz]/np.abs(r[nz])).max() if nz.any() else 0
    W = np.abs(a['W'] - b['W']).max()/np.abs(b['W']).max()
    print('%s spec bitwise %s  max|d|/max|ref| %.2e  max rel (nonzero px) %.2e  n_diff_px %d / %d  W rel %.1e  ptan eq %s' % (
        tag, np.array_equal(a['spec'], b['spec']), d.max()/np.abs(r).max(), rel, (d > 0).sum(), d.size, W,
        np.array_equal(a['ptan'], b['ptan'])))
    try:
        x, y = np.load('%s/ana_%s.npz' % (new, tag), allow_pickle=True), np.load('%s/ana_%s.npz' % (ref, tag), allow_pickle=True)
    except FileNotFoundError:
        continue
    for k in ('full', 'rmcut'):
        kx, ky = x['kpm_' + k].item(), y['kpm_' + k].item()
        print('   %s Kp-marg median new %.4f ref %.4f | b_peak %.4f %.4f | kp0 peak %.4f %.4f | mapmax %s %s' % (
            k, kx['q'][1], ky['q'][1], kx['b_peak'], ky['b_peak'], x['kp0_' + k].item()['peak'], y['kp0_' + k].item()['peak'],
            np.round(x['mapmax_' + k], 3), np.round(y['mapmax_' + k], 3)))
    M = 'M_full'; print('   Kp-Vsys map max|d|/max %.2e' % (np.abs(x[M] - y[M]).max()/np.abs(y[M]).max()))
