"""moveq-1011: dynamic diffusion (IDORT 4.2.1 / J21 4.4) diagnostics.
usage: ddts.py <rundir>...  per output: peak/analytic, centroid - v t, variance rate (D_eff/D),
L1 vs analytic, total E, gas rho range."""
import glob
import sys
import numpy as np

L, BG, V = 20.0, np.exp(-10.0), 1.0


def rd(f):
    t = float(open(f).readline().split('time=')[1].split()[0])
    hdr = [h for h in open(f).readlines()[1].lstrip('#').split() if h not in ('k', 'x3v', 'j', 'x2v')]
    a = np.loadtxt(f)
    return t, {h: a[:, n] for n, h in enumerate(hdr)}


def ana(d):
    c, ks = 1000.0, 4.0e4
    for ln in open(d + '/in.athinput'):
        if ln.strip().startswith('c_light'):
            c = float(ln.split('=')[1])
        if ln.strip().startswith('kappa_s'):
            ks = float(ln.split('=')[1])
    D = c/(3*ks)
    rows = []
    for f in sorted(glob.glob(d + '/tab/*.m1.*.tab')):
        t, m = rd(f)
        x, E = m['x1v'], m['m1_e']
        s = 1 + 160*D*t
        xx = (x - V*t + L/2) % L - L/2
        ex = np.exp(-40*xx**2/s)/np.sqrt(s) + BG
        w = np.maximum(E - BG, 0)
        xc = np.angle(np.sum(w*np.exp(2j*np.pi*x/L)))*L/(2*np.pi)
        xr = (x - xc + L/2) % L - L/2
        var = np.sum(w*xr**2)/np.sum(w)
        dx = L/len(x)
        rho = (np.nan, np.nan)
        fh = f.replace('.m1.', '.hydro_w.')
        try:
            _, h = rd(fh)
            rho = (h['dens'].min(), h['dens'].max())
        except Exception:
            pass
        rows.append((t, E.max()/ex.max(), ((xc - V*t + L/2) % L - L/2), var, s/80,
                     np.abs(E - ex).sum()/np.abs(ex - BG).sum(), E.sum()*dx, rho[0], rho[1]))
    return rows, D


if __name__ == '__main__':
    for d in sys.argv[1:]:
        rows, D = ana(d)
        print('== %s' % d)
        print('    t   peak/ex  xc-vt     var    var_ex  L1      sumE     rho_min  rho_max')
        for r in rows:
            print('%5.1f  %.4f  %+.4f  %.5f %.5f %.4f  %.5f  %.5f %.5f' % r)
        if len(rows) > 2:
            dv = (rows[-1][3] - rows[1][3])/(rows[-1][0] - rows[1][0])
            print('   D_eff/D (variance slope, t>=%g) = %.4f' % (rows[1][0], dv/(2*D)))
