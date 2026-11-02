"""ana.py: CCF, Kp-Vsys (Seidel+2026) and morning/evening segments (Seidel+2025) from
OUTDIR/spec_<tag>.npz (fe3d.py).

Per exposure: excess absorption spectrum (planet frame) -> continuum removal (fe_core.hipass, as
the template) -> instrumental Gaussian FWHM = planet lsf_fwhm_kms (4.28 for WASP-121b) -> CCF
with the pRT 3000 K Fe I template (broadened by the same LSF):
CCF(v) = sum F(lam) T(lam/(1+v/c)) / sum T, lags +-60 km/s on the R = 1e6 pixel grid (0.2998 km/s).
Kp-Vsys: exposures shifted by (Kp0 - Kp) sin(2 pi phi) and summed (equal weights), as kp.py.
Kp-marg (Seidel+2026 Methods: "Given the uncertainty on the orbital velocity, we marginalise over
the orbital velocity as well as over the planetary wind velocity"; no Kp range/prior stated):
  (a) every Kp row in the prior window (planet key kp_marg_range, default Kp0 +- 50, flat prior) gets
      a Gaussian+offset fit (+-12 km/s around its maximum, as 'peak'); noise sigma = map std at
      |v| > 30 (ED Fig. 2 significance normalisation; or scaled so the map maximum = planet key
      kp_marg_snr); A = row maximum (|v| < 35) minus the row median at |v| > 30; row likelihood
      ln L = (A/sigma)^2/2 (matched-filter delta chi^2); the offset
      posterior is the L-weighted mixture of N(mu(Kp), sigma_mu(Kp)) -> mean, 16/50/84 %.
  (b) the map rows averaged over the same Kp window, then the usual peak/wide fits.
  python ana.py --planet planets/wasp121b.json --out OUTDIR <tag> [<tag> ...]
      -> OUTDIR/ana_<tag>.npz, summary printed and appended to OUTDIR/ana_summary.txt
"""
import argparse
import os
import sys
import warnings
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fe_core as F  # noqa: E402
import planet        # noqa: E402

NL = 200
VL = F.DV/1e5*np.arange(-NL, NL + 1)          # km/s
VSG = np.arange(-40, 40.001, 0.1)
_T = {}


def configure(cfg):
    """Kp0, RM mask, T2/T3, LSF from planet.load() (WASP-121b: Kp0 219.2, 0.017, 4.28 km/s)."""
    global KP0, KPG, PHRM, PH14, PH23, LSF, KPM_RANGE, KPM_SNR
    KP0 = cfg['KP0']
    KPM_RANGE = cfg.get('kp_marg_range') or [KP0 - 50.0, KP0 + 50.0]
    KPM_SNR = cfg.get('kp_marg_snr')
    KPG = np.arange(KP0 - 100, KP0 + 100.01, 1.0)
    PHRM = cfg['rm_mask_phase']                # Seidel+2025 RM exclusion |phi| < 0.017
    PH14 = cfg['PH14']                         # T1 / T4
    PH23 = cfg['PH23']                         # T2 / T3
    LSF = cfg['LSF']
    _T.clear()


configure(planet.load())


def template():
    if 'tb' not in _T:
        _T['tb'] = np.maximum(F.template_hp(LSF), 0)
        _T['ker'] = F.lsf_kernel(LSF)
    return _T['tb'], _T['ker']


def ccf(spec):
    tb, ker = template()
    f = np.convolve(F.hipass(spec.astype(np.float64)), ker, mode='same')
    n = len(f)
    out = np.empty(2*NL + 1)
    for i, k in enumerate(range(-NL, NL + 1)):
        if k >= 0:
            out[i] = np.dot(f[k:], tb[:n - k])/tb[:n - k].sum()
        else:
            out[i] = np.dot(f[:n + k], tb[-k:])/tb[-k:].sum()
    return out


def gfit(v, y, v0, hw, s0):
    from scipy.optimize import curve_fit
    m = abs(v - v0) < hw
    f = lambda x, A, mu, s, c: A*np.exp(-0.5*((x - mu)/s)**2) + c   # noqa: E731
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            p, _ = curve_fit(f, v[m], y[m], p0=[y[m].max() - y[m].min(), v0, s0, y[m].min()],
                             maxfev=20000)
        except RuntimeError:
            return np.nan, np.nan
    return p[1], 2.3548*abs(p[2])


def metrics(v, row):
    """peak-Gaussian (+-12 km/s around the maximum, as kp.py / Seidel 'highest peak'),
    single Gaussian over +-35 km/s, centroid over |v| < 25 after removing the |v| > 30 median."""
    vpk = v[np.argmax(np.where(abs(v) < 35, row, -np.inf))]
    pk, _ = gfit(v, row, vpk, 12.0, 3.0)
    wide, fw = gfit(v, row, 0.0, 35.0, 6.0)
    base = np.median(row[(abs(v) > 30) & (abs(v) < 40)])
    m = abs(v) < 25
    y = np.maximum(row[m] - base, 0)
    cen = (v[m]*y).sum()/y.sum()
    return dict(peak=pk, wide=wide, fwhm=fw, cen=cen, vmax=vpk)


def kpvsys(cc, phases, sel):
    sn = np.sin(2*np.pi*phases)
    M = np.zeros((len(KPG), len(VSG)))
    for i, kp in enumerate(KPG):
        for j in np.flatnonzero(sel):
            M[i] += np.interp(VSG - (KP0 - kp)*sn[j], VL, cc[j])
    return M


def gfit_cov(v, y, v0, hw, s0, sig):
    """Gaussian+offset fit with absolute noise sig: (A, mu, sigma_mu); nan on failure."""
    from scipy.optimize import curve_fit
    m = abs(v - v0) < hw
    f = lambda x, A, mu, s, c: A*np.exp(-0.5*((x - mu)/s)**2) + c   # noqa: E731
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            p, C = curve_fit(f, v[m], y[m], p0=[y[m].max() - y[m].min(), v0, s0, y[m].min()],
                             sigma=np.full(m.sum(), sig), absolute_sigma=True, maxfev=20000)
        except RuntimeError:
            return np.nan, np.nan, np.nan
    return p[0], p[1], np.sqrt(abs(C[1, 1]))


def kpmarg(M):
    """Kp-marginalised offset of a Kp-Vsys map M (KPG x VSG); see the module docstring."""
    w = (KPG >= KPM_RANGE[0] - 1e-9) & (KPG <= KPM_RANGE[1] + 1e-9)
    A = np.array([M[i][abs(VSG) < 35].max() - np.median(M[i][abs(VSG) > 30])
                  for i in np.flatnonzero(w)])
    sig = M[:, abs(VSG) > 30].std()
    if KPM_SNR:
        sig = A.max()/float(KPM_SNR)
    mu, smu = np.full(w.sum(), np.nan), np.full(w.sum(), np.nan)
    for n, i in enumerate(np.flatnonzero(w)):
        vpk = VSG[np.argmax(np.where(abs(VSG) < 35, M[i], -np.inf))]
        _, mu[n], smu[n] = gfit_cov(VSG, M[i], vpk, 12.0, 3.0, sig)
    ok = np.isfinite(mu) & np.isfinite(smu) & (smu > 0)
    lnl = 0.5*(A/sig)**2
    wt = np.where(ok, np.exp(lnl - np.nanmax(lnl[ok])), 0.0)
    wt /= wt.sum()
    vg = np.arange(-40, 40.0001, 0.005)
    pdf = sum(wt[n]*np.exp(-0.5*((vg - mu[n])/smu[n])**2)/smu[n] for n in np.flatnonzero(ok))
    cdf = np.cumsum(pdf)/pdf.sum()
    q = [float(np.interp(x, cdf, vg)) for x in (0.16, 0.5, 0.84)]
    row = M[w].mean(0)
    mb = metrics(VSG, row)
    return dict(mean=float((wt*np.where(ok, mu, 0)).sum()), q=q, sig=float(sig),
                smax=float(np.nanmax(A[ok])/sig), kpmean=float((wt*KPG[w]).sum()),
                flat=float(np.nanmean(mu[ok])), b_peak=mb['peak'], b_wide=mb['wide'],
                b_fwhm=mb['fwhm'], b_cen=mb['cen'], range=list(map(float, KPM_RANGE)),
                kp=KPG[w], mu=mu, smu=smu, wt=wt)


def analyse(tag, outdir):
    d = np.load(os.path.join(outdir, 'spec_%s.npz' % tag))
    ph = d['phases']
    spec = d['spec']                              # (nph, 2, NPIX)
    cc = np.array([[ccf(spec[j, 0]), ccf(spec[j, 1])] for j in range(len(ph))])
    cb = cc.sum(1)                                # both limbs
    R = {'tag': tag, 'rot': float(d['rot']), 'phases': ph, 'cc': cc}
    edge = np.isclose(abs(ph), PH14)              # exposures exactly at T1 / T4 (zero depth)
    sets = {'full': ~edge, 'rmcut': (~edge) & (abs(ph) > PHRM)}
    for sn_, sel in sets.items():
        M = kpvsys(cb, ph, sel)
        i0 = np.argmin(abs(KPG - KP0))
        R['M_' + sn_] = M
        R['kp0_' + sn_] = metrics(VSG, M[i0])
        im, jm = np.unravel_index(np.argmax(M), M.shape)
        R['mapmax_' + sn_] = (KPG[im], VSG[jm])
        R['kpm_' + sn_] = kpmarg(M)
    # Seidel+2025 segments: morning = ingress epoch, evening = egress epoch, |phi| > 0.017 (RM
    # mask) and "excluded true ingress and egress" (Table 3 text), i.e. PHRM < |phi| < T2/T3.
    # Their Table 3 gives the AVERAGE per-exposure shift and FWHM; also the stacked-CCF fit.
    for seg, sel in (('mor', (ph < -PHRM) & (ph > -PH23)), ('eve', (ph > PHRM) & (ph < PH23))):
        st = cb[sel].sum(0)
        R['seg_' + seg] = metrics(VL, st)
        R['seg_' + seg + '_n'] = int(sel.sum())
        R['seg_' + seg + '_ccf'] = st
        R['kpm_' + seg] = kpmarg(kpvsys(cb, ph, sel))
        per = np.array([[metrics(VL, cb[j])['wide'], metrics(VL, cb[j])['fwhm']]
                        for j in np.flatnonzero(sel)])
        R['seg_' + seg + '_perexp'] = float(np.nanmean(per[:, 0]))
        R['seg_' + seg + '_perexp_fwhm'] = float(np.nanmean(per[:, 1]))
        R['seg_' + seg + '_perexp_sd'] = float(np.nanstd(per[:, 0]))
        for li, ln in ((0, 'Lmor'), (1, 'Leve')):
            R['seg_%s_%s' % (seg, ln)] = metrics(VL, cc[sel, li].sum(0))['wide']
    # contribution: template-weighted absorption vs tangent log p (all phases, both limbs)
    W, pt = d['W'], d['ptan']
    edges = np.arange(-9, 1.01, 0.25)
    R['cf_edges'] = edges
    R['cf'] = np.histogram(pt.ravel(), edges, weights=W.ravel())[0]
    R['cf_lim'] = [np.histogram(pt[:, np.cos(d['th']) > 0 if k == 0 else np.cos(d['th']) <= 0].ravel(),
                                edges, weights=W[:, np.cos(d['th']) > 0 if k == 0 else
                                                 np.cos(d['th']) <= 0].ravel())[0] for k in (0, 1)]
    R['Wmid'] = W[len(ph)//2]
    R['ptmid'] = pt[len(ph)//2]
    R['th'] = d['th']
    R['bb'] = d['bb']
    c = np.cumsum(R['cf'])/R['cf'].sum()
    lc = 0.5*(edges[1:] + edges[:-1])
    R['cf_p16_50_84'] = [float(np.interp(q, c, lc)) for q in (0.16, 0.5, 0.84)]
    np.savez(os.path.join(outdir, 'ana_%s.npz' % tag), **{k: (np.array(v, dtype=object) if isinstance(v, dict)
                                                  else v) for k, v in R.items()})
    return R


def summary(R):
    s = ['%s rot %.2f' % (R['tag'], R['rot'])]
    for sn_ in ('full', 'rmcut'):
        k = R['kp0_' + sn_]
        s.append('  KpVsys %-5s at Kp0: peak %+6.2f  wide %+6.2f (FWHM %.1f)  centroid %+6.2f  '
                 'map max Kp %.1f v %+.2f' % (sn_, k['peak'], k['wide'], k['fwhm'], k['cen'],
                                              *R['mapmax_' + sn_]))
    for seg in ('mor', 'eve'):
        k = R['seg_' + seg]
        s.append('  seg %s (n=%d): mean per-exp %+6.2f (sd %.2f) FWHM %.1f | stacked: wide %+6.2f '
                 'FWHM %.1f peak %+6.2f centroid %+6.2f | limb rays mor %+6.2f eve %+6.2f' % (
                     seg, R['seg_%s_n' % seg], R['seg_%s_perexp' % seg],
                     R['seg_%s_perexp_sd' % seg], R['seg_%s_perexp_fwhm' % seg], k['wide'],
                     k['fwhm'], k['peak'], k['cen'], R['seg_%s_Lmor' % seg],
                     R['seg_%s_Leve' % seg]))
    for k in ('full', 'rmcut', 'mor', 'eve'):
        r = R['kpm_' + k]
        s.append('  Kp-marg %-5s Kp %.1f-%.1f: (a) posterior mean %+6.2f median %+6.2f [16 %+6.2f, '
                 '84 %+6.2f] <Kp> %.1f (S_max %.1f, flat-Kp mean %+6.2f) | (b) Kp-avg map: peak '
                 '%+6.2f wide %+6.2f (FWHM %.1f) centroid %+6.2f' % (
                     k, r['range'][0], r['range'][1], r['mean'], r['q'][1], r['q'][0], r['q'][2],
                     r['kpmean'], r['smax'], r['flat'], r['b_peak'], r['b_wide'], r['b_fwhm'],
                     r['b_cen']))
    s.append('  CF (template-weighted absorption) log p_tan 16/50/84 %%: %.2f %.2f %.2f' %
             tuple(R['cf_p16_50_84']))
    return '\n'.join(s)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--planet', help='planet JSON (default: WASP-121b values)')
    ap.add_argument('--out', required=True, help='directory holding spec_<tag>.npz')
    ap.add_argument('tags', nargs='+')
    a = ap.parse_args()
    configure(planet.load(a.planet))
    L = []
    for tag in a.tags:
        R = analyse(tag, a.out)
        L.append(summary(R))
        print(L[-1], flush=True)
    with open(os.path.join(a.out, 'ana_summary.txt'), 'a') as f:
        f.write('\n'.join(L) + '\n')
