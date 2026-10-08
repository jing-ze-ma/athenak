import numpy as np, premix_lib as p
out = np.load('out/res_1e8.npy', allow_pickle=True)
TN = np.array([o[0] for o in out])
KT = p.KT
iT = [list(KT['T']).index(t) for t in TN]
K = KT['K'][iT, 0]                       # (nT, nb, ng)
def st(d):
    d = np.abs(d[np.isfinite(d)])
    return np.median(d), np.percentile(d, 90), d.max()
gm = (TN >= 500) & (TN <= 3000)
for var in ('base', 'atomzero', 'co2x', 'ext'):
    for ci, cn in ((0, 'Gauss-value'), (1, 'subint-mean')):
        k = np.array([o[1][var][ci] for o in out])
        d = np.log10(np.maximum(k, 1e-300)/K)
        print(f'{var:9s} {cn:12s} GATE(500-3000,b2-b10) med/p90/max = %.3f %.3f %.2f' % st(d[gm][:, 2:11]),
              ' signed median %.3f' % np.median(d[gm][:, 2:11]))
k = np.array([o[1]['base'][0] for o in out]); d = np.log10(np.maximum(k, 1e-300)/K)
np.save('out/dlog_base.npy', d)
print('per band (base, Gauss-value, 500-3000 K): med p90 max signed-med')
for b in range(11):
    print(b, '%.2f-%.2f um' % (KT['wl'][b+1], KT['wl'][b]), '%.3f %.3f %.2f' % st(d[gm][:, b]), '%.3f' % np.median(d[gm][:, b]))
print('per T range (b2-b10)')
for lo, hi in ((100, 400), (500, 900), (1000, 1500), (1700, 2300), (2500, 2900), (3100, 4100), (4300, 6100)):
    s = (TN >= lo) & (TN <= hi)
    print(lo, hi, '%.3f %.3f %.2f' % st(d[s][:, 2:11]), '%.3f' % np.median(d[s][:, 2:11]))
print('per g (500-3000, b2-b10)')
for g in range(8):
    print(g, '%.3f %.3f %.2f' % st(d[gm][:, 2:11, g]), '%.3f' % np.median(d[gm][:, 2:11, g]))
print('signed dlog, base, Gauss value; rows T, cols band, g-mean over g0-3 | g4-7')
for i, t in enumerate(TN):
    print('%5d ' % t + ' '.join('%5.2f' % x for x in d[i, :, :4].mean(1)) + ' | ' + ' '.join('%5.2f' % x for x in d[i, :, 4:].mean(1)))
