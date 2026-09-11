"""si/s05: radial profiles of the deepest drained columns, rot 1..5, and the same
columns in si/s01 for reference."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/analysis')
import bin_convert
import dhjcs
import csgeom
B = '/viper/u2/jinma/ATHENAK/bench/cs_ens/'
EOS = dhjcs.EOS(B+'analysis/eos_table.txt')
ROT = 3.05e5
NP, NC = 6, 32
RV, _ = csgeom.x1v(128)


def load(run, n):
    r = bin_convert.read_binary('%s%s/bin/dhj.hydro_w.%05d.bin' % (B, run, n))
    mb = r['mb_data']
    g = np.asarray(r['mb_geometry'])
    d = np.asarray(mb['dens']).astype(float)
    e = np.asarray(mb['eint']).astype(float)
    v = np.asarray(mb['velx']).astype(float)
    D = np.zeros((NP, NC, NC, 128))
    E = np.zeros_like(D)
    V = np.zeros_like(D)
    for m in range(24):
        p = m//4
        j0 = 0 if g[m, 2] < -0.5 else 16
        k0 = 0 if g[m, 4] < -0.5 else 16
        D[p, j0:j0+16, k0:k0+16] = np.transpose(d[m], (1, 0, 2))
        E[p, j0:j0+16, k0:k0+16] = np.transpose(e[m], (1, 0, 2))
        V[p, j0:j0+16, k0:k0+16] = np.transpose(v[m], (1, 0, 2))
    T, _ = EOS.invert(D, E)
    return r['time'], D, T, V


out = open('profiles_s05.txt', 'w')


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    out.write(s+'\n')


S = {n: load('si/s05', n) for n in range(1, 6)}
COLS = [(0, 18, 14, 'gid 1  deepest drained at rot 5'),
        (0, 14, 11, 'gid 0  2nd deepest'),
        (0, 21, 18, 'gid 3  coldest of the drained set')]
for (p, J, K, lab) in COLS:
    P('\n=== si/s05 column panel %d J %d K %d  (%s) ===' % (p, J, K, lab))
    P('%6s %11s | %s' % ('i', 'r[cm]', ' '.join('%25s' %
      ('rot %d: rho / T / v_r' % n) for n in range(1, 6))))
    for i in list(range(28, 80, 4))+list(range(80, 128, 8)):
        row = '%6d %11.4e |' % (i, RV[i])
        for n in range(1, 6):
            t, D, T, V = S[n]
            row += ' %8.2e %7.1f %9.2e' % (D[p, J, K, i], T[p, J, K, i], V[p, J, K, i])
        P(row)
    P('%6s %11s |' % ('', '') + ''.join(' rot%d minrho %.2e' %
      (n, S[n][1][p, J, K].min()) for n in range(1, 6)))
P('')
P('=== column-integrated mass of the drained columns, '
  'rot 1..5 (arbitrary weight = sum rho) ===')
P('%-22s %10s %10s %10s %10s %10s' % ('column', 'rot1', 'rot2', 'rot3', 'rot4', 'rot5'))
for (p, J, K, lab) in COLS:
    P('%-22s ' % ('p%d J%d K%d' % (p, J, K)) + ' '.join('%10.3e' %
      S[n][1][p, J, K].sum() for n in range(1, 6)))
out.close()
