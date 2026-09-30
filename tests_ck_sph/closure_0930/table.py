"""Compact VARIANTS.md table.  usage: python3 table.py SUF ..."""
import sys
sys.path.insert(0, '/viper/ptmp2/jinma/cksph_test_0930/ana')
from summ import exact, grad, newton   # noqa: E402
for suf in sys.argv[1:]:
    sh = [exact(c, suf)[0] for c in ['tr_x1.2', 'tr_x1.7', 'tr_x3', 'st_x1.2', 'st_x1.7']]
    row = [suf, ' / '.join(f'{v:.4f}' for v in sh)]
    for c in ['w_T2000', 'w_T3000', 'o_T2000', 'o_T3000']:
        tot, band, dep = exact(c, suf)
        b9 = {k: v for k, v in band.items() if k < 10}
        wb = max(b9, key=lambda k: abs(b9[k] - 1))
        row.append(f'{tot:.4f} (b{wb} {b9[wb]:.3f}; b10 {band[10]:.3f})')
    for c in ['w_T2000', 'w_T3000', 'o_T2000', 'o_T3000']:
        dep = exact(c, suf)[2]
        q = dep[dep[:, 3] < 40][:, 6]
        row.append(f'{q[:-1].min():.3f}-{q[:-1].max():.3f} [{q[-1]:.2f}]')
    n, mx, mn, rd, rv = grad(suf)
    row.append(f'{mn:+.1e} (max {mx:.1e})')
    n, p, nc, bm, bl = newton(suf)
    row.append(f'{p}/{nc}/{bm:.1e}')
    print('| ' + ' | '.join(row) + ' |')
