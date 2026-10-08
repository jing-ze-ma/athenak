"""write Premixed_1x_g8_11_hiT2_dace_lowP_SiO{fc,tab}.txt.  k_new = k_tab(T,1e-8) x R,
R = k_premix(T,P)/k_premix(T,1e-8) (1 where the denominator is <= 0/non-finite; floored at
1e-30, no cap); T > 6100 K: the 6100 K R; running max in g on the new rows."""
import sys
import numpy as np
sys.path.insert(0, '/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/scripts')
import ck_lib as ck  # noqa
SRC = '/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/data/ck/Premixed_1x_g8_11_hiT2.txt'
D = '/orion/ptmp/jinma/rsg_wind_1008/dace/'
z = np.load(D + 'out/premix_lowP.npz')
kt = ck.read_ktable(SRC)
K = kt['K']
nT, nP, nb, ng = K.shape
lines = open(SRC).read().split('\n')
assert len(lines) == 10 + nT*nP*nb + 1 and lines[-1] == ''
PN = z['P'][:18]
T = kt['T']
imap = [int(np.argmin(abs(z['T']-t))) if t <= 6100 else len(z['T'])-1 for t in T]
assert all(z['T'][i] == t for i, t in zip(imap, T) if t <= 6100)
for var, key, tag in (('SiOfc', 'KA', 'FastChem SiO'), ('SiOtab', 'KB', 'SiO x table-implied f(T)')):
    Kp = z[key][imap]                                  # (nT, 19, nb, ng)
    den = Kp[:, 18:19]
    with np.errstate(divide='ignore', invalid='ignore'):
        R = Kp[:, :18]/den
    bad = ~np.isfinite(R) | ~(den > 0)
    nbad = int(bad.sum())
    R[bad] = 1.
    nfl = int((R < 1e-30).sum())
    R = np.maximum(R, 1e-30)
    Knew = K[:, :1]*R
    Km = np.maximum.accumulate(Knew, axis=-1)
    nmono = int((Km != Knew).sum())
    Knew = Km
    TAG = (f' | dace_lowP ({tag}): P < 1e-8 bar = k(1e-8) x k_premix(P)/k_premix(1e-8), DACE 0.01 cm-1 '
           'line premix at 1e-8 bar cross sections, FastChem 4.0.3 eq-cond VMRs (rsg_wind_1008/dace)')
    fn = D + f'data/ck/Premixed_1x_g8_11_hiT2_dace_lowP_{var}.txt'
    with open(fn, 'w') as f:
        f.write(lines[0] + TAG + '\n')
        f.write(f'{nT} {nP + 18} {nb} {ng}\n')
        f.write(lines[2] + '\n')
        f.write(' '.join([repr(float(p)) for p in PN]) + ' ' + lines[3] + '\n')
        for q in range(4, 10):
            f.write(lines[q] + '\n')
        for it in range(nT):
            for ip in range(18):
                for b in range(nb):
                    f.write(' '.join(f'{v:.10e}' for v in Knew[it, ip, b]) + '\n')
            q0 = 10 + it*nP*nb
            for q in range(q0, q0 + nP*nb):
                f.write(lines[q] + '\n')
    print(var, 'entries', R.size, 'denominator guard', nbad, 'floored', nfl, 'g-runmax changed', nmono,
          'R>1:', int((R > 1).sum()), 'max R %.3g' % R.max())
