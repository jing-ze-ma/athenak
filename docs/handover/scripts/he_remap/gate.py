# flake8: noqa
"""Gates for he_remap_rst.py.  usage:
  gate.py cmp A.rst B.rst              section-wise max |A-B|/max|A| (active angles, all r)
  gate.py cons SRC.rst OUT.rst         domain totals (mass, U = E - rho Phi, m_r, E_rad, F_r r^2)
  gate.py hse SRC.rst OUT.rst [tag]    shell-mean HSE residual before/after -> hse_<tag>.npz
"""
import he_remap_rst as H
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_mltpp_1002/remap128')
sys.path.insert(0, '/viper/ptmp2/jinma/he_mltpp_1002')

GM = 4.18143e26


def geom(r):
    p = [float(H.get_param(r.text, 'mesh', k) or 0.0) for k in H.SKEYS]
    f = H.faces_of(p, r.bind[1], r.ng, r.msize[0], r.msize[3])
    N2, N3 = r.mind[2], r.mind[3]
    th = np.linspace(r.msize[1], r.msize[4], N2 + 1)
    dA = (np.cos(th[:-1]) - np.cos(th[1:]))[None, :] * \
        np.full((N3, 1), (r.msize[5] - r.msize[2]) / N3)
    a = slice(r.ng, r.ng + r.bind[1])
    fa = f[r.ng:r.ng + r.bind[1] + 1]
    dV = dA[:, :, None] * (np.diff(fa**3) / 3.0)[None, None, :]
    xv = 0.75 * np.diff(fa**4) / np.diff(fa**3)
    return a, fa, xv, dA, dV


def load(fn):
    r = H.read_rst(fn)
    G = H.gather(r)
    return r, G


def cmp(fa, fb):
    ra, Ga = load(fa)
    rb, Gb = load(fb)
    for k in Ga:
        if k not in Gb:
            print('  %-5s only in A' % k)
            continue
        d = np.abs(Ga[k] - Gb[k]).max() / max(np.abs(Ga[k]).max(), 1e-300)
        print('  %-5s max|A-B|/max|A| = %.3e   bitwise equal: %s' %
              (k, d, np.array_equal(Ga[k], Gb[k])))


def totals(fn):
    r, G = load(fn)
    a, fa, xv, dA, dV = geom(r)
    u = G['hyd']
    rho = u[0][..., a]
    phi = GM / r.msize[0] - GM / xv
    U = u[4][..., a] - rho * phi[None, None, :]
    return dict(mass=(rho *
                      dV).sum(), U=(U *
                                    dV).sum(), E=(u[4][..., a] *
                                                  dV).sum(), m_r=(u[1][..., a] *
                                                                  dV).sum(), Lz_like=(u[3][..., a] *
                                                                                      xv *
                                                                                      dV).sum(), E_rad=(G['m1'][0][..., a] *
                                                                                                        dV).sum(), F_r=(G['m1'][1][..., a] *
                                                                                                                        dV).sum())


def cons(fs, fo):
    a, b = totals(fs), totals(fo)
    for k in a:
        print('  %-8s src %.12e out %.12e rel %.2e' %
              (k, a[k], b[k], (b[k] - a[k]) / abs(a[k])))


def hse(fs, fo, tag):
    import ana3d as A3
    out = {}
    for nm, fn in (('src', fs), ('out', fo)):
        r, G = load(fn)
        a, fa, xv, dA, dV = geom(r)
        u = G['hyd']
        rho = u[0][..., a]
        ke = 0.5 * (u[1][..., a]**2 + u[2][..., a]**2 + u[3][..., a]**2) / rho
        phi = GM / r.msize[0] - GM / xv
        eint = u[4][..., a] - ke - rho * phi[None, None, :]
        T = A3.temperature(rho, eint)
        lt = np.log10(np.clip(T, 10**3.8, 10**6.5))
        lr = np.clip(np.log10(rho), -16, -2)
        pg = rho * 10**A3.HC.SLP.ev(lt, lr)
        if nm == 'src' and 'wd' in G:
            pf = G['wd'][0][..., a]
            m = rho > 1e-12
            out['pg_py_vs_file'] = float(np.abs(pg[m] / pf[m] - 1).max())
        E = G['m1'][0][..., a]
        w = dA[:, :, None] / dA.sum()
        rm = (rho * w).sum(axis=(0, 1))
        P = ((pg + E / 3.0) * w).sum(axis=(0, 1))
        rf = 0.5 * (rm[1:] + rm[:-1])
        xf = 0.5 * (xv[1:] + xv[:-1])
        res = (np.diff(P) / np.diff(xv) + rf * GM / xf**2) / (rf * GM / xf**2)
        out[nm + '_x'] = xf
        out[nm + '_res'] = res
        out[nm + '_rho'] = rm
        out[nm + '_P'] = P
        out[nm + '_xv'] = xv
    np.savez('/viper/ptmp2/jinma/he_mltpp_1002/remap128/gate/hse_%s.npz' % tag, **out)
    R = 9.4868e10 / 0.4
    if 'pg_py_vs_file' in out:
        print(
            '  python EOS p_gas vs the file cache (rho > 1e-12): max rel %.2e' %
            out['pg_py_vs_file'])
    print('  shell-mean HSE residual (dP_tot/dr + <rho> g)/(<rho> g), P_tot = <p_g> + <E>/3:')
    print('   band R      src median |res|  src max   out median  out max   |out-src| median (interp)')
    for lo, hi in ((0.40, 0.45), (0.45, 0.55), (0.55, 0.60), (0.60, 0.65),
                   (0.65, 0.70), (0.70, 0.80), (0.80, 0.90)):
        ms = (out['src_x'] > lo * R) & (out['src_x'] < hi * R)
        mo = (out['out_x'] > lo * R) & (out['out_x'] < hi * R)
        di = np.abs(
            np.interp(
                out['src_x'][ms],
                out['out_x'],
                out['out_res']) -
            out['src_res'][ms])
        print(
            '   %.2f-%.2f  %.3e  %.3e  %.3e  %.3e  %.3e' %
            (lo, hi, np.median(
                np.abs(
                    out['src_res'][ms])), np.abs(
                out['src_res'][ms]).max(), np.median(
                    np.abs(
                        out['out_res'][mo])), np.abs(
                            out['out_res'][mo]).max(), np.median(di)))


if __name__ == '__main__':
    c = sys.argv[1]
    if c == 'cmp':
        cmp(sys.argv[2], sys.argv[3])
    elif c == 'cons':
        cons(sys.argv[2], sys.argv[3])
    elif c == 'hse':
        hse(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else 'x')
