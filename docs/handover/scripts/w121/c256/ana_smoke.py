"""smoke transient check: area-weighted isobar T at 1e-3/0.1/10 bar and rms v_r per band,
source rot-300 bin vs the smoke bins.  usage: ana_smoke.py RUNDIR"""
import sys, glob
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dhjcs  # noqa
eos = dhjcs.EOS('/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt')
PR = 1.101535e5
LEV = np.log10([1e-3, 0.1, 10.0]) + 6.0
BANDS = ((1e-4, 1e-2), (1e-2, 1.0), (1.0, 10.0), (10.0, 100.0), (100.0, 1e4))


def one(fn):
    raw = dhjcs.bin_convert.read_binary(fn)
    G = np.asarray(raw['mb_geometry'])
    D = raw['mb_data']
    nmb = len(G)
    sT = np.zeros(len(LEV)); sw = np.zeros(len(LEV))
    sv = np.zeros(len(BANDS)); nv = np.zeros(len(BANDS))
    for m in range(nmb):
        d = np.asarray(D['dens'][m], float)
        ei = np.asarray(D['eint'][m], float)
        vr = np.asarray(D['velx'][m], float)
        T, p = eos.invert(d, ei)
        nk, nj, ni = d.shape
        x = np.tan(np.pi/4*(G[m, 2] + (G[m, 3] - G[m, 2])*(np.arange(nj) + 0.5)/nj))
        y = np.tan(np.pi/4*(G[m, 4] + (G[m, 5] - G[m, 4])*(np.arange(nk) + 0.5)/nk))
        X, Y = np.meshgrid(x, y)
        w = (1 + X*X)*(1 + Y*Y)/(1 + X*X + Y*Y)**1.5*(2.0/nj)*(2.0/nk)
        Tl = dhjcs.level_interp(T[None], np.log10(p)[None], LEV)
        for i in range(len(LEV)):
            q = Tl[i][0]; ok = np.isfinite(q)
            sT[i] += (np.where(ok, q, 0)*w).sum(); sw[i] += (w*ok).sum()
        pb = p/1e6
        for b, (lo, hi) in enumerate(BANDS):
            s = (pb >= lo) & (pb < hi)
            sv[b] += (vr[s]**2).sum(); nv[b] += s.sum()
    return [raw['time']/PR] + list(sT/sw) + list(np.sqrt(sv/np.maximum(nv, 1))/1e5)


fs = ['/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.00150.bin'] + \
    sorted(glob.glob(sys.argv[1] + '/bin/*.bin'))
print('# file  rot  T(1e-3) T(0.1) T(10 bar) [K]  rms v_r [km/s] in bands',
      ' '.join('%g-%g' % b for b in BANDS))
for f in fs:
    r = one(f)
    print('%-26s %.5f ' % (f.split('/')[-1] if 'w1x' not in f else 'SOURCE C32x76', r[0])
          + ' '.join('%8.2f' % v for v in r[1:4]) + '  '
          + ' '.join('%7.4f' % v for v in r[4:]))
