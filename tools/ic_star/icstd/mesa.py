"""MESA profile readers (comparison + optional L/M/Teff/R source), tools/ic_star.

format 'mesa': a MESA profile*.data (header names on line 2, values line 3, column names line 6;
  logR [Rsun], logT, logRho, mass [Msun], luminosity [Lsun]; r, m, L at the zone's outer face,
  T and rho at its centre -> centre radius ((r_out^3 + r_in^3)/2)^(1/3), as
  make_ic_giant_mesa.py).
format 'cgs6': r T rho m kappa F_r (cgs; L. Ma's MESA_profile_combined_cgs.txt).
Constants as make_ic_giant_mesa.py (MESA): Msun 1.9892e33, Rsun 6.957e10, Lsun 3.828e33.
"""
import numpy as np

MSUN_MESA = 1.9892e33
RSUN_MESA = 6.957e10
LSUN_MESA = 3.828e33


def read(path, fmt='mesa', r_max_rsun=None):
    m = _read(path, fmt)
    if r_max_rsun:                  # drop rows above r_max (e.g. a non-envelope top row)
        k = m['r'] <= r_max_rsun*RSUN_MESA
        for key, v in list(m.items()):
            if isinstance(v, np.ndarray) and v.size == k.size:
                m[key] = v[k]
        m['note'] += '; rows above %g Rsun dropped (spec r_max_rsun)' % r_max_rsun
    return m


def _read(path, fmt='mesa'):
    if fmt == 'cgs6':
        d = np.loadtxt(path)
        d = d[np.concatenate([[True], np.diff(d[:, 0]) > 0])]
        r, T, rho, m, kap, F = d.T
        L = 4*np.pi*r**2*F
        Lph = float(np.median(L))
        return dict(r=r, rc=r, T=T, rho=rho, m=m, L=L, kap=kap, L_ph=Lph, M_tot=float(m[-1]),
                    R_ph=None, Teff=None, fmt=fmt, path=path,
                    note='cgs6: L_ph = median 4 pi r^2 F_r (radiative profile); no photosphere'
                         ' header')
    lines = open(path).read().splitlines()
    hn, hv = lines[1].split(), lines[2].split()
    hdr = {}
    for k, v in zip(hn, hv):
        try:
            hdr[k] = float(v.replace('D', 'E'))
        except ValueError:
            pass
    names = lines[5].split()
    d = np.loadtxt(path, skiprows=6)
    p = {n: d[:, i] for i, n in enumerate(names)}
    rm = 10**p['logR']*RSUN_MESA
    o = np.argsort(rm)
    rm = rm[o]
    keep = np.concatenate([[True], np.diff(rm) > 0])
    rm = rm[keep]

    def sel(k):
        return p[k][o][keep]
    rcen = ((rm**3 + np.concatenate([[0.0], rm[:-1]])**3)/2.0)**(1.0/3.0)
    out = dict(r=rm, rc=rcen, T=10**sel('logT'), rho=10**sel('logRho'),
               m=sel('mass')*MSUN_MESA, L=sel('luminosity')*LSUN_MESA, fmt=fmt, path=path)
    out['M_tot'] = float(out['m'][-1])
    out['L_ph'] = hdr.get('photosphere_L', out['L'][-1]/LSUN_MESA)*LSUN_MESA
    out['R_ph'] = hdr.get('photosphere_r', rm[-1]/RSUN_MESA)*RSUN_MESA
    out['Teff'] = hdr.get('Teff')
    if 'logP' in p:
        out['P'] = 10**sel('logP')
    if 'opacity' in p:
        out['kap'] = sel('opacity')
    for k in ('x_mass_fraction_H', 'y_mass_fraction_He', 'z_mass_fraction_metals', 'h1',
              'he4'):
        if k in p:
            out[k] = sel(k)
    out['note'] = ('MESA header: Teff %s, photosphere_L %s Lsun, photosphere_r %s Rsun, '
                   'M %.6g Msun' % (hdr.get('Teff'), hdr.get('photosphere_L'),
                                   hdr.get('photosphere_r'), out['M_tot']/MSUN_MESA))
    return out


def at(m, key, r, centre=False):
    """log-interpolated profile quantity at radius r (T, rho at centres)"""
    x = m['rc'] if centre else m['r']
    return float(np.exp(np.interp(np.log(r), np.log(x), np.log(m[key]))))
