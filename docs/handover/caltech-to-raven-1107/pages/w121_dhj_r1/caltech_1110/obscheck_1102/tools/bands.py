"""bands.py: em3d outputs -> planet/star flux ratio vs phase and wavelength, band phase curves,
amplitude, hotspot offset, day/night brightness temperatures.
  python bands.py SYSTEM_KEY OUT_NPZ em_1.npz [em_2.npz ...]   (dumps averaged)
Star: PHOENIX (pRT 3.4 table, Husser+2013), surface flux at T_eff, bin-averaged; R* from SYS.
Band ratios are photon-weighted (int F lam R dlam) with top-hat responses (BANDS).
Brightness temperature: Planck T that gives the band ratio with the OBSERVED Rp/Rs and the same
PHOENIX star (the observers' definition).
"""
import sys
import numpy as np

H, C, KB = 6.62607015e-27, 2.99792458e10, 1.380649e-16
RSUN = 6.957e10
SYS = {  # literature system values (the fe_hires planet JSONs; Teff/Rp/Rs below)
    'w121': dict(name='WASP-121b', Teff=6460.0, Rs=1.458, RpRs=0.1236),    # Delrez+2016 (as Mikal-Evans+23, Morello+23, Splinter+25)
    'w76': dict(name='WASP-76b', Teff=6329.0, Rs=1.756, RpRs=0.10852),     # Ehrenreich+2020
}
BANDS = {'IRAC1': (3.18, 3.94), 'IRAC2': (4.0, 5.0), 'G395H': (2.87, 5.18), 'G395H_NRS1': (2.87, 3.72),
         'G395H_NRS2': (3.82, 5.18), 'WFC3': (1.12, 1.65), 'SOSS': (0.85, 2.85), 'TESS': (0.6, 1.0)}


def star_flux(Teff, lam_um, edges_um):
    from petitRADTRANS.stellar_spectra.phoenix import phoenix_star_table
    s, _ = phoenix_star_table.compute_spectrum(Teff)
    lc = s[:, 0]
    Fl = s[:, 1]*C/lc**2                       # F_lambda [erg/s/cm^2/cm]
    o = np.argsort(lc)
    lc, Fl = lc[o]*1e4, Fl[o]
    out = np.empty(len(lam_um))
    for j, (a, b) in enumerate(edges_um):
        m = (lc >= a) & (lc < b)
        out[j] = Fl[m].mean() if m.sum() > 2 else np.interp(lam_um[j], lc, Fl)
    return out


def planck(lam_cm, T):
    return 2*H*C**2/lam_cm**5/np.expm1(H*C/(lam_cm*KB*T))


def tb(ratio, lam_um, w, Fs, RpRs):
    """band brightness temperature: solve sum w B(T) / sum w Fs/pi... = ratio/RpRs^2."""
    from scipy.optimize import brentq
    lc = lam_um*1e-4
    tgt = ratio/RpRs**2*(w*Fs).sum()/np.pi
    f = lambda T: (w*planck(lc, T)).sum() - tgt   # noqa: E731
    try:
        return brentq(f, 100, 8000)
    except ValueError:
        return np.nan


def main():
    key, outf = sys.argv[1], sys.argv[2]
    S = SYS[key]
    E = [np.load(f) for f in sys.argv[3:]]
    lam = E[0]['lam_um']
    ph = E[0]['phases']
    Fp = np.mean([e['Fp'] for e in E], 0)               # sum I dA [erg/s/cm/sr * cm^2]
    tab = np.load(str(E[0]['table']))
    edges = tab['lam_edges_um']
    Fs = star_flux(S['Teff'], lam, edges)
    Rs = S['Rs']*RSUN
    ratio = Fp/(Fs*Rs**2)                               # (phase, w)
    rall = np.array([e['Fp'] for e in E])/(Fs*Rs**2)
    dl = edges[:, 1] - edges[:, 0]
    R = dict(phases=ph, lam=lam, edges=edges, ratio=ratio, ratio_sd=rall.std(0), Fs=Fs)
    for bn, (a, b) in BANDS.items():
        w = ((lam >= a) & (lam < b))*lam*dl
        if w.sum() == 0:
            continue
        cs = (rall*Fs*w).sum(-1)/(Fs*w).sum()          # (dump, phase)
        c = cs.mean(0)
        # Fourier fit (2nd order) for offset/amplitude, as the observers
        A = np.stack([np.ones_like(ph), np.cos(2*np.pi*ph), np.sin(2*np.pi*ph), np.cos(4*np.pi*ph),
                      np.sin(4*np.pi*ph)], 1)
        co = np.linalg.lstsq(A, c, rcond=None)[0]
        pf = np.linspace(0, 1, 2001)
        Af = np.stack([np.ones_like(pf), np.cos(2*np.pi*pf), np.sin(2*np.pi*pf), np.cos(4*np.pi*pf),
                       np.sin(4*np.pi*pf)], 1)
        cf = Af@co
        pmax, pmin = pf[np.argmax(cf)], pf[np.argmin(cf)]
        day, night = np.interp(0.5, pf, cf), cf[0]
        R['pc_' + bn] = c
        R['pcsd_' + bn] = cs.std(0)
        R['st_' + bn] = dict(day=day, night=night, amp=cf.max() - cf.min(), max=cf.max(), min=cf.min(),
                             off_max_deg=(0.5 - pmax)*360, off_min_deg=(0.0 - (pmin if pmin < 0.5 else pmin - 1))*360,
                             Tday=tb(day, lam, w, Fs, S['RpRs']), Tnight=tb(max(night, 1e-12), lam, w, Fs, S['RpRs']))
        print('%-10s day %6.0f ppm night %6.0f ppm amp %6.0f ppm  peak offset %+5.1f deg (east +)  min offset %+5.1f'
              '  Tday %5.0f K  Tnight %5.0f K' % (bn, day*1e6, night*1e6, (cf.max() - cf.min())*1e6,
                                                   R['st_' + bn]['off_max_deg'], R['st_' + bn]['off_min_deg'],
                                                   R['st_' + bn]['Tday'], R['st_' + bn]['Tnight']))
    # spectral brightness temperature at day / night
    lc = lam*1e-4
    for nm, i in (('day', np.argmin(abs(ph - 0.5))), ('night', 0)):
        x = ratio[i]/S['RpRs']**2*Fs/np.pi
        R['Tb_' + nm] = H*C/(lc*KB)/np.log1p(2*H*C**2/lc**5/np.maximum(x, 1e-300))
    np.savez(outf, **{k: (np.array(v, dtype=object) if isinstance(v, dict) else v) for k, v in R.items()},
             sys=str(S))


if __name__ == '__main__':
    main()
