"""rw-drag-1112: multi-mode (acoustic + slow thermal-diffusion + left acoustic + fast) fit of the
x1-line Fourier amplitudes of (rho, v, T, E, F1) at the dumps.  At P 3.33 tau 10 the slow mode
(omega ~ -2.3i) decays slower than the acoustic one, so the 1-mode phase/log fit of rho is biased by
~1e-4 relative; the acoustic omega of this fit is the one to compare with the reference root."""
import glob
import os
import numpy as np


def _fh(x, a, k):
    return complex(((a - a.mean()) * np.exp(-1j * k * x)).mean())


def series(rdir, k):
    import common  # noqa: E402
    ts, Z = [], []
    for ph in sorted(glob.glob(os.path.join(rdir, "tab", "rw.hydro_w.*.tab"))):
        n = ph.split(".")[-2]
        h = common.load_tab(ph)
        m = common.load_tab(os.path.join(rdir, "tab", "rw.m1.%s.tab" % n))
        x = np.asarray(h.x1)
        d = h.var("dens")
        T = (2.0 / 3.0) * h.var("eint") / d      # gamma 5/3, P = rho T
        Z.append([_fh(x, d, k), _fh(x, h.var("velx"), k), _fh(x, T, k),
                  _fh(x, m.var("m1_e"), k), _fh(x, m.var("m1_f1"), k)])
        ts.append(h.time)
    return np.array(ts), np.array(Z)


def fit(ts, Z, wref, vref, guesses=None):
    """-> (omegas[4], amplitudes A[mode, var] in units of the ref eigvec, rms)."""
    from scipy.optimize import least_squares
    sc = np.abs(vref)
    Zs = Z / Z[0, 0] / sc
    g = guesses if guesses is not None else [wref, -2.38j, -np.conjugate(wref), -103j]

    def model(p):
        ws = p[0::2] + 1j * p[1::2]
        B = np.exp(-1j * np.outer(ts, ws))
        A, *_ = np.linalg.lstsq(B, Zs, rcond=None)
        return B, A, ws

    def res(p):
        B, A, _ = model(p)
        r = (B @ A - Zs).ravel()
        return np.concatenate([r.real, r.imag])
    p0 = np.ravel([[complex(w).real, complex(w).imag] for w in g])
    s = least_squares(res, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    B, A, ws = model(s.x)
    return ws, A, float(np.sqrt(np.mean(s.fun ** 2)))


def acoustic(rdir, bg, wref, vref):
    """dict for results.json: w4_re/w4_im, slow mode, its rho amplitude, acoustic eigvec/ref."""
    ts, Z = series(rdir, bg.k)
    ws, A, rms = fit(ts, Z, wref, vref)
    ev = (A[0] / A[0, 0]) * np.abs(vref) / (vref / vref[0])
    return {"w4_re": float(ws[0].real), "w4_im": float(ws[0].imag), "w4_rms": rms,
            "w4_slow": [float(ws[1].real), float(ws[1].imag)],
            "w4_slow_amp": float(abs(A[1, 0] / A[0, 0])),
            "w4_eig_ratio": [[float(z.real), float(z.imag)] for z in ev]}
