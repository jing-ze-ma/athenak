"""tr_fe.py: time-resolved Fe I (per exposure) from fe_hires outputs.

Input: fe_hires ana_<tag>.npz (key 'cc' = per-exposure, per-limb CCF with the pRT 3000 K Fe I
template, planet rest frame, instrument LSF of the planet JSON) -- or, if no ana file, the
spec_<tag>.npz of fe3d (old CPU or fehires_fast GPU path, same format), CCF'd here with ana.ccf.
Per exposure: Gaussian+offset fit (+-35 km/s, as ana.metrics 'wide') to
  both limbs (= what an observer measures), morning limb (theta with cos > 0... fe3d limb 0),
  evening limb (fe3d limb 1)
-> centre v (km/s, + = redshift), amplitude A (CCF units = template-weighted excess absorption,
reported in ppm of the stellar flux), FWHM. Mean and sd over the dumps.
  python tr_fe.py PLANET_JSON OUT_NPZ DIR tag1 [tag2 ...]
"""
import os, sys, warnings
import numpy as np
FH = os.environ.get('FEHIRES_CODE', '/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires/code')
sys.path.insert(0, FH)
import planet  # noqa: E402
import ana     # noqa: E402


def gfit(v, y, hw=35.0):
    from scipy.optimize import curve_fit
    m = abs(v) < hw
    vpk = v[m][np.argmax(y[m])]
    f = lambda x, A, mu, s, c: A*np.exp(-0.5*((x - mu)/s)**2) + c   # noqa: E731
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            p, _ = curve_fit(f, v[m], y[m], p0=[y[m].max() - np.median(y[m]), vpk, 4.0, np.median(y[m])],
                             maxfev=20000)
        except RuntimeError:
            return np.nan, np.nan, np.nan
    return p[0], p[1], 2.3548*abs(p[2])


def main():
    pj, outf, d = sys.argv[1], sys.argv[2], sys.argv[3]
    tags = sys.argv[4:]
    cfg = planet.load(pj)
    ana.configure(cfg)
    R = []
    for tag in tags:
        fa = os.path.join(d, 'ana_%s.npz' % tag)
        if os.path.exists(fa):
            z = np.load(fa, allow_pickle=True)
            cc, ph = z['cc'], z['phases']
        else:
            z = np.load(os.path.join(d, 'spec_%s.npz' % tag))
            ph = z['phases']
            cc = np.array([[ana.ccf(z['spec'][j, 0]), ana.ccf(z['spec'][j, 1])] for j in range(len(ph))])
        out = np.full((len(ph), 3, 3), np.nan)     # exposure, (both, mor, eve), (A, v, fwhm)
        for j in range(len(ph)):
            for k, y in enumerate((cc[j].sum(0), cc[j, 0], cc[j, 1])):
                if np.ptp(y) > 0:
                    out[j, k] = gfit(ana.VL, y)
        R.append(out)
        print(tag, 'both-limb v at T1+..T4-:', ' '.join('%+.1f' % x for x in out[1:-1:3, 0, 1]), flush=True)
    R = np.array(R)
    R[..., 0] *= 1e6                                 # amplitude -> ppm
    np.savez(outf, phases=ph, mean=np.nanmean(R, 0), sd=np.nanstd(R, 0), all=R, tags=tags,
             PH14=cfg['PH14'], PH23=cfg['PH23'], KP0=cfg['KP0'], LSF=cfg['LSF'], cc_mean=None)


if __name__ == '__main__':
    main()
