"""What heats the 1x equatorial knot (lon +100..+125, 1e-5..1e-4 bar)? Equatorial (|lat|<8) profiles vs lon of T, zonal
wind u, zonal Mach u/c_s (c_s from the EOS Gamma1), radial velocity v_r, mean molecular weight mu (recombination tracer),
mean of dumps 146-150 of w1x (rot ~300)."""
import sys, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts'); sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dhjcs, dclib
eos = dhjcs.EOS('/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt'); eg = dclib.EOS('/viper/ptmp2/jinma/deepconv_0925/dump/eos_grid.txt')
PB = [1e-5, 1e-4, 1e-3]; LEV = np.log10(PB) + 6
acc = None
for b in range(146, 151):
    d = dhjcs.load('/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.%05d.bin' % b, eos, 'cs')
    rho = d['rho']; T = d['T']
    mu = eg.at('mu', rho, T); chr_ = eg.at('chr', rho, T); cht = eg.at('cht', rho, T); cv = eg.at('cv', rho, T)
    g1 = chr_ + d['p']*cht**2/(rho*T*cv)
    cs = np.sqrt(g1*d['p']/rho)
    vr = np.asarray(d['raw']['mb_data']['velx'])
    Q = [T, d['u'], d['u']/cs, vr, mu, g1]
    L = [np.asarray(dhjcs.level_interp(q, d['lp'], LEV)) for q in Q]
    la, lo, M = dhjcs.latlon_bin(np.stack(L), d['lat'], d['lon'], 4, 4)
    acc = M if acc is None else acc + M
M = acc/5
eq = M[..., np.abs(la) < 8, :].mean(-2)          # (q, level, lon)
names = ['T[K]', 'u[km/s]', 'Mach_u', 'v_r[m/s]', 'mu', 'Gamma1']
sc = [1, 1e-5, 1, 1e-2, 1, 1]
for k, p in enumerate(PB):
    print('== %.0e bar, |lat|<8 deg' % p)
    print('   lon  ' + ''.join('%10s' % n for n in names))
    for i in np.where((lo >= 30) & (lo <= 178))[0][::2]:
        print('  %+5.0f ' % lo[i] + ''.join('%10.3g' % (eq[q, k, i]*sc[q]) for q in range(6)))
np.savez('knot_eq.npz', lat=la, lon=lo, M=M, names=names, pbar=PB)
