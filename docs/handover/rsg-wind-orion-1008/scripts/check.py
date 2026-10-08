import os, sys
os.environ['CK_DATA'] = '/orion/ptmp/jinma/rsg_wind_1008/dace/data/'
sys.path.insert(0, '/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/scripts')
import numpy as np
import ck_lib as ck
D = ck.DATA + 'ck/'
fo = D + 'Premixed_1x_g8_11_hiT2.txt'
names = {'old': fo, 'lowP': D + 'Premixed_1x_g8_11_hiT2_lowP.txt',
         'A': D + 'Premixed_1x_g8_11_hiT2_dace_lowP_SiOfc.txt', 'B': D + 'Premixed_1x_g8_11_hiT2_dace_lowP_SiOtab.txt'}
kt = {k: ck.read_ktable(v) for k, v in names.items()}
o = kt['old']
lo = open(fo).read().split('\n')
for k in ('A', 'B'):
    t = kt[k]
    ln = open(names[k]).read().split('\n')
    print(k, 'header', ln[1], '| nT nP nb ng', t['K'].shape, '| P[0] %.3g P[17] %.3g P[18] %.3g' % (t['P'][0], t['P'][17], t['P'][18]),
          '| T,wl,wn,g identical:', all(np.array_equal(t[x], o[x]) for x in ('T', 'wl', 'wn', 'gx', 'gw')))
    nT, nP, nb = 93, 34, 11
    same = all(ln[10 + it*(nP+18)*nb + 18*nb + q] == lo[10 + it*nP*nb + q] for it in range(nT) for q in range(nP*nb))
    print('   rows P>=1e-8 textually identical to original:', same, '; K[:,18:] == K_orig:', np.array_equal(t['K'][:, 18:], o['K']),
          '; min K new rows %.3g; g-monotone new rows:' % t['K'][:, :18].min(), bool((np.diff(t['K'][:, :18], axis=-1) >= 0).all()))
ce = {'old': ck.read_ce(ck.DATA + 'CE_tables/FastChem_ck_1x_int_hiT2.txt'),
      'new': ck.read_ce(ck.DATA + 'CE_tables/FastChem_ck_1x_int_hiT2_lowP.txt')}
cia = ck.read_cia(); ray = ck.read_ray()
print('means(): kR, kP (cm^2/g)')
for T, P in ((1500., 1e-6), (1500., 1e-10), (2500., 1e-12), (3000., 3e-9), (1000., 1e-13)):
    r = []
    for k in ('old', 'lowP', 'A', 'B'):
        m = ck.means(kt[k], ce['old' if k == 'old' else 'new'], cia, ray, T, P)
        r.append('%s %.3e %.3e' % (k, m['kR'], m['kP']))
    print('  T %5d P %.0e : ' % (T, P) + ' | '.join(r))
# continuity: band-mean (sum gw k) log slope d log k / d log P just above / below 1e-8
print('continuity: dlogk/dlogP of band-mean k, above (1e-8..2.19e-8, table) vs below (4.68e-9..1e-8); median over bands, and max |diff|')
for k in ('A', 'B'):
    t = kt[k]; P = t['P']; K = (t['K']*t['gw']).sum(-1)
    up = np.log10(K[:, 19]/K[:, 18])/np.log10(P[19]/P[18]); dn = np.log10(K[:, 18]/K[:, 17])/np.log10(P[18]/P[17])
    for T in (1000., 1500., 2000., 2500., 3000., 3500.):
        i = int(np.argmin(abs(t['T']-T)))
        print('  %s T %5d up %s | dn %s' % (k, t['T'][i], ' '.join('%5.2f' % x for x in up[i]), ' '.join('%5.2f' % x for x in dn[i])))
# band-mean ratios at 1e-10 bar relative to 1e-8 (kline from means: bilinear log interp in T)
print('band-mean line k ratio k(1e-10)/k(1e-8) per band b0..b10 [lowP | A | B]')
for T in (1000., 1500., 2000., 2500., 3000., 3500.):
    rows = []
    for k in ('lowP', 'A', 'B'):
        a = ck.means(kt[k], ce['new'], cia, ray, T, 1e-10)['kline']; b = ck.means(kt[k], ce['new'], cia, ray, T, 1e-8)['kline']
        rows.append(a/b)
    print('T %d' % T)
    for k, r in zip(('lowP', 'A', 'B'), rows):
        print('   %-4s' % k, ' '.join('%8.2e' % x for x in r))
    print('   A/lowP', ' '.join('%8.2f' % x for x in rows[1]/rows[0]), '  B/A', ' '.join('%5.2f' % x for x in rows[2]/rows[1]))
