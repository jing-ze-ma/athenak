import sys, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhj_remap as R
out = {}
for a in ('L', 'H'):
    h = R.read_rst('%s/rst/dhj.00610.rst' % a); ng = h['ng']
    s = (slice(None), slice(ng, -ng), slice(ng, -ng), slice(ng, -ng))
    u = h['u']; rho = u[:, 0][s]; p = h['p'][s]
    vh = np.sqrt(u[:, 2][s]**2 + u[:, 3][s]**2)/rho   # covariant-ish, proxy only
    out[a] = dict(pb=p.mean((0, 1, 2))/1e6, pr=(p/rho).mean((0, 1, 2)), vh=np.sqrt((vh**2).mean((0, 1, 2))), vmax=vh.max((0, 1, 2)), E=u[:, 4][s].sum())
L, H = out['L'], out['H']
for i in range(0, len(L['pb']), 6):
    print('%3d p %.3g bar  d(p/rho)/ %+.2e  vh_rms L %.4g H %.4g  vh_max L %.4g H %.4g' % (i, L['pb'][i], H['pr'][i]/L['pr'][i]-1, L['vh'][i]/100, H['vh'][i]/100, L['vmax'][i]/100, H['vmax'][i]/100))
print('max |d(p/rho)| (p>1e-3 bar):', np.abs(H['pr']/L['pr']-1)[L['pb']>1e-3].max())
