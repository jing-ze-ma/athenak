"""coverage.py (dhj_planets_1101): ck / CE / EOS table coverage of an IC column and of an expected
day/night range (the IC column with T x 1.3 (day) and x 0.6 (night) above 0.01 bar, the bounding
columns planet_setup.py reports, plus a hot dayside bound T_day,max = 2^(1/2) x Teq x 1.3 for
zero-redistribution substellar heating).
usage: coverage.py IC_PROFILE EOS_DUMP_DIR CK_TABLE CE_TABLE TEQ"""
import sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dclib  # noqa: E402

ic, eosd, ckf, cef, teq = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])


def grids(fn):
    with open(fn) as fh:
        fh.readline()
        fh.readline()
        T = np.array(fh.readline().split(), float)
        p = np.array(fh.readline().split(), float)
    return T, p


Tck, pck = grids(ckf)
Tce, pce = grids(cef)
eos = dclib.EOS(eosd + '/eos_grid.txt')
d = np.loadtxt(ic)
pbar, T = 10**d[:, 0]/1e6, d[:, 1]
print('IC %s: p %.2e..%.2e bar, T %.0f..%.0f K' % (ic.split('/')[-1], pbar.min(), pbar.max(),
                                                  T.min(), T.max()))
print('ck table T %.0f..%.0f K, p %.0e..%.0e bar; CE table T %.0f..%.0f K, p %.0e..%.0e bar'
      % (Tck.min(), Tck.max(), pck.min(), pck.max(), Tce.min(), Tce.max(), pce.min(), pce.max()))
print('code EOS table: log rho -19..0, log T 1.5..6 (eos_logd_min -19, eos_logt_min/max 1.5/6)')
for lab, fac in (('IC', 1.0), ('day (T x 1.3 above 0.01 bar)', 1.3),
                 ('night (T x 0.6 above 0.01 bar)', 0.6)):
    Tc = np.where(pbar < 0.01, T*fac, T)
    rho = np.ravel(eos.rho_pt(pbar*1e6, Tc))
    act = (pbar >= 1e-6) & (pbar <= 300)
    inck = (Tc >= Tck.min()) & (Tc <= Tck.max())
    inp = (pbar >= pck.min()) & (pbar <= pck.max())
    print('%-32s T %.0f..%.0f K (p >= 1e-6 bar: %.0f..%.0f); rho %.2e..%.2e; ck T in range: %s; '
          'rows above the ck p top (clamped, as W121): %d of %d (all at p < %.0e bar)'
          % (lab, Tc.min(), Tc.max(), Tc[act].min(), Tc[act].max(), rho.min(), rho.max(),
             bool(inck.all()), (~inp).sum(), len(pbar), pck.min()))
tsub = np.sqrt(2.0)*teq*1.3
print('hot-day bound sqrt2 Teq x 1.3 = %.0f K: inside the ck T range: %s; T > 6100 K uses the hiT2 '
      'extension (gen_hitemp.py, H bf/ff)' % (tsub, tsub <= Tck.max()))
