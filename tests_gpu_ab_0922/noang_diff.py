"""Q6: is `noang` different from `sph` at all?  Compares the two hst time series at
IDENTICAL times (the only exactly comparable data) and the last dumps."""
import sys
import glob
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert                       # noqa: E402
import dhjcs                             # noqa: E402
import common as C                       # noqa: E402

NAMES = ['time', 'dt', 'mass', '1-mom', '2-mom', '3-mom', 'tot-E',
         '1-KE', '2-KE', '3-KE', '1-ME', '2-ME', '3-ME']


def main():
    a = C.hst(C.B + 'ck_sph_ab/sph')
    b = C.hst(C.B + 'cs_noang/noang')
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    same_t = np.allclose(a[:, 0], b[:, 0], rtol=0, atol=1e-3)
    fh = open(C.OUT + 'noang_diff.md', 'w')
    fh.write('# noang vs sph (noang_diff.py)\n\n')
    fh.write('hst rows compared: %d, times identical: %s\n\n' % (n, same_t))
    fh.write('| column | max |rel diff| | first row with rel diff > 1e-10 (t, rot since '
             'restart) |\n|---|---|---|\n')
    for i, nm in enumerate(NAMES):
        with np.errstate(invalid='ignore', divide='ignore'):
            r = np.abs(a[:, i] - b[:, i]) / np.maximum(np.abs(a[:, i]), 1e-300)
        k = np.where(r > 1e-10)[0]
        first = ('%.6e (%.4f)' % (a[k[0], 0], (a[k[0], 0] - 8.64e7) / C.PROT)
                 if len(k) else 'never')
        fh.write('| %s | %.3e | %s |\n' % (nm, np.nanmax(r), first))
    fh.write('\ndt: identical in every row: %s\n'
             % bool(np.all(a[:, 1] == b[:, 1])))

    # last dumps (different cycles, 20 s apart) -- structural comparison only
    fs = sorted(glob.glob(C.B + 'ck_sph_ab/sph/bin/*.bin'))[-1]
    fn = sorted(glob.glob(C.B + 'cs_noang/noang/bin/*.bin'))[-1]
    rs = bin_convert.read_binary(fs)
    rn = bin_convert.read_binary(fn)
    fh.write('\nlast dumps: sph t=%.6e cycle %d, noang t=%.6e cycle %d '
             '(NOT the same instant, so a cellwise diff mixes the physics difference '
             'with %0.f s of evolution)\n\n'
             % (rs['time'], rs['cycle'], rn['time'], rn['cycle'],
                abs(rs['time'] - rn['time'])))
    geo = dhjcs.cs_geometry(rs)
    fh.write('| variable | rms rel diff | max rel diff | location of max (lat, lon, i) '
             '|\n|---|---|---|---|\n')
    for v in rs['var_names']:
        x = np.asarray(rs['mb_data'][v], dtype=np.float64)
        y = np.asarray(rn['mb_data'][v], dtype=np.float64)
        sc = np.sqrt(np.mean(x**2))
        d = np.abs(x - y) / sc
        idx = np.unravel_index(np.argmax(d), d.shape)
        fh.write('| %s | %.3e | %.3e | %+.0f, %+.0f, %d |\n'
                 % (v, np.sqrt(np.mean(d**2)), d.max(),
                    np.degrees(geo['lat'][idx[0], idx[1], idx[2]]),
                    np.degrees(geo['lon'][idx[0], idx[1], idx[2]]), idx[3]))
    fh.close()
    print(open(C.OUT + 'noang_diff.md').read())


main()
