# wasp189b: copy of wasp121_0925/transit.py with a CHEOPS passband (0.35-1.10 um, top-hat) added first
"""Transit (chord) radius of a 1-D RCE column with the dhj model's own ck opacities.

Usage: python3 transit.py <rce.npz> [dT]      (env RCE_* as for rce_w121.py: RCE_MET, RCE_EOSDIR)
  rce.npz : rce_w121.py output (lpc = log10 p_c [barye], T, rf = face radii [cm]).
  dT      : optional uniform temperature offset [K] applied above 0.1 bar (limb sensitivity).

Method.  The column is made spherically symmetric (global mean, see README for the limb
caveat).  For every ck chain c = (band b, g-point) the chord optical depth at impact
parameter y is tau_c(y) = sum_i k_i L_i(y) with k = kappa rho of the cell (ck line +
continuum + Rayleigh + CIA + H-: rr_lib CK.kr, the code's own opacity) and L_i the exact
chord length through shell i.  Band transmission t_b(y) = sum_g w_g exp(-tau_bg(y))
(correlated-k), effective transit radius per band
    R_b^2 = r_0^2 + 2 int_{r_0}^{r_top} (1 - t_b(y)) y dy       (r_0 = the 250 bar anchor),
and a bandpass radius from the photon-weighted depths, R^2 = sum_b W_b R_b^2 / sum_b W_b,
W_b = stellar photons in (band b) cap (bandpass).  p_ref = p(R).  The tau = 0.56 chord
level of the band-mean transmission (Lecavelier des Etangs+2008) is printed alongside.
"""
import os
import sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/viper/ptmp2/jinma/wasp121_0925')
os.environ.setdefault('RCE_TSTAR', '6000.0')
import rce_w121 as R  # noqa: E402

H_, C_, K_ = 6.62607015e-27, 2.99792458e10, 1.380649e-16
PASSBANDS = {'CHEOPS 0.35-1.10um': (0.35, 1.10), 'NIRSpec NRS1 2.70-3.72um': (2.70, 3.72), 'TESS 0.6-1.0um': (0.60, 1.00), 'white 0.61-0.85um (band 8)': (0.61, 0.85),
             'z 0.82-0.92um': (0.82, 0.92), 'WFC3 G141 1.1-1.7um': (1.1, 1.7),
             'NIRSpec G395H 2.9-5.2um': (2.9, 5.2)}


def photons(l0, l1, teff, n=2000):
    """stellar photon count between l0 and l1 [um] (blackbody at teff; relative)."""
    lam = np.linspace(l0, l1, n)*1e-4
    x = H_*C_/(lam*K_*teff)
    return np.trapezoid(lam**-4/np.expm1(x), lam)


def chord_radii(rf, kr, gw):
    """rf (n+1,) faces, kr (n, 88).  Returns y (impact parameters), R_b (NB,), t_b(y)."""
    n = len(rf) - 1
    y = np.concatenate([rf[:-1] + 0.5*np.diff(rf)*f for f in np.linspace(0, 1, 5)[:-1]])
    y = np.sort(y)
    ro, ri = rf[1:][None, :], rf[:-1][None, :]
    Y = y[:, None]
    L = 2*(np.sqrt(np.maximum(ro**2 - Y**2, 0)) - np.sqrt(np.maximum(np.maximum(ri, Y)**2
                                                                      - Y**2, 0)))
    tau = L @ kr                                          # (ny, 88)
    t = np.exp(-tau).reshape(len(y), R.NB, R.NG) @ gw     # (ny, NB)
    r0 = rf[0]
    ys = np.concatenate([[r0], y, [rf[-1]]])
    ts = np.concatenate([np.zeros((1, R.NB)), t, t[-1:]], 0)
    Rb = np.sqrt(r0**2 + 2*np.trapezoid((1 - ts)*ys[:, None], ys, axis=0))
    return y, Rb, t, tau


def main(fn, dT=0.0, teff=6586.0, verbose=True):
    Z = np.load(fn)
    lpc, T, rf = Z['lpc'], Z['T'].copy(), Z['rf']
    pc = 10**lpc
    if dT != 0.0:
        T = np.where(pc < 0.1*R.BAR, T + dT, T)
        # re-integrate the hydrostatic radii with the same anchor as the RCE (r = AP at
        # 250 bar, g = G0 (AP/r)^2)
        col = R.Col(20)
        col.geom(T)
        rf = col.rf
    rho = R.rho_pt(pc, T)
    kr = R.ck.kr(T, pc/R.BAR, rho)
    y, Rb, t, tau = chord_radii(rf, kr, R.ck.gw)
    rc = 0.5*(rf[1:] + rf[:-1])
    p_at = lambda r: 10**np.interp(r, rc, lpc)/R.BAR          # noqa: E731
    wl = R.ck.wl
    out = {}
    for name, (l0, l1) in PASSBANDS.items():
        W = np.array([photons(max(l0, wl[b+1]), min(l1, wl[b]), teff)
                      if min(l1, wl[b]) > max(l0, wl[b+1]) else 0.0 for b in range(R.NB)])
        Rp = np.sqrt(np.sum(W*Rb**2)/W.sum())
        tb = t @ (W/W.sum())
        # tau_chord = 0.56 level of the bandpass-mean transmission
        y56 = np.interp(np.exp(-0.56), tb, y)
        out[name] = (Rp, p_at(Rp), y56, p_at(y56), W/W.sum())
        if verbose:
            print('%-28s R = %.5e cm  p_ref = %.3e bar | tau_ch=0.56: r = %.5e p = %.3e bar'
                  '  weights %s' % (name, Rp, p_at(Rp), y56, p_at(y56),
                                    ' '.join('%d:%.2f' % (b, w) for b, w in
                                             enumerate(W/W.sum()) if w > 0.005)))
    if verbose:
        print('per band: ' + ' '.join('%.2f-%.2fum %.2e bar;' % (wl[b+1], wl[b], p_at(Rb[b]))
                                      for b in range(R.NB)))
        print('H at p_ref(TESS) = %.3e cm' % (np.interp(out['TESS 0.6-1.0um'][0], rc,
                                                         pc/(rho*R.G0*(R.AP/rc)**2))))
    return out, rc, lpc, Rb


if __name__ == '__main__':
    dT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    teff = float(os.environ.get('TR_TEFF', '6628.0'))
    main(sys.argv[1], dT, teff)
