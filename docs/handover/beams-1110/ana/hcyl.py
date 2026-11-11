"""suite-1011: homogeneous cylinder (Cartesian 2-D analogue of IDORT 4.1.2 / J21 4.6 homogeneous sphere).
Uniform cylinder r < R0 = 1 (rho kappa = 10, B = a T^4 = 1) in vacuum, box [-4,4]^2, x1 outflow (vacuum),
x2 periodic (period 8: an infinite row of cylinders).  Exact: I(n) = sum over images of
B (1 - exp(-tau_k)) exp(-tau of images nearer), tau_k = chi * chord_k / sin(theta).
Moments normalised so that an isotropic B gives E = 1 (E = <I>, F/c = <I n>, P = <I n n>).
usage: python hcyl.py exact            -> ana/hcyl_exact.npz (cut along y = 0, x in [0, 4])
       python hcyl.py <rundir> ...     -> metrics
"""
import sys
import os
import glob
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/blend_1022/wt/vis/python')
CHI, R0, L, PER = 10.0, 1.0, 4.0, 8.0
HERE = os.path.dirname(os.path.abspath(__file__))


def exact_point(x, y, nmu=64, nphi=4096, kmax=60):
    mu, wmu = np.polynomial.legendre.leggauss(nmu)       # mu = cos(theta) about z
    phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    st = np.sqrt(1 - mu**2)
    ux, uy = np.cos(phi), np.sin(phi)                       # in-plane propagation direction
    # trace backwards: points P - s u, s >= 0, until |x| > L
    with np.errstate(divide='ignore'):
        sx = np.where(ux > 0, (x + L) / np.where(ux > 0, ux, 1), np.where(ux < 0, (L - x) / np.where(ux < 0, -ux, 1), np.inf))
    sx = np.minimum(sx, 1e4)
    Iin = np.zeros(nphi)
    trans = np.ones((nmu, nphi))
    I = np.zeros((nmu, nphi))
    # images k: centre (0, PER k); chord of line P - s u with circle, s in [0, sx]
    segs = []
    for k in range(-kmax, kmax + 1):
        cx, cy = 0.0 - x, PER * k - y                        # centre relative to P
        b = cx * ux + cy * uy                                 # roots s = -b +- sqrt(b^2 - cc)
        cc = cx * cx + cy * cy - R0 * R0
        disc = b * b - cc
        ok = disc > 0
        sq = np.sqrt(np.where(ok, disc, 0))
        s1 = np.clip(-b - sq, 0, sx)
        s2 = np.clip(-b + sq, 0, sx)
        ch = np.where(ok, s2 - s1, 0.0)
        segs.append((s1, ch))
    s1 = np.array([s[0] for s in segs])
    ch = np.array([s[1] for s in segs])
    order = np.argsort(s1, axis=0)
    s1 = np.take_along_axis(s1, order, 0)
    ch = np.take_along_axis(ch, order, 0)
    for k in range(ch.shape[0]):
        tau = CHI * ch[k][None, :] / st[:, None]
        I += trans * (1 - np.exp(-tau))
        trans *= np.exp(-tau)
    w = (wmu[:, None] / 2) * (np.ones(nphi)[None, :] / nphi)
    nx = st[:, None] * ux[None, :]
    ny = st[:, None] * uy[None, :]
    E = np.sum(w * I)
    F = np.sum(w * I * nx)
    Pxx = np.sum(w * I * nx * nx)
    Pyy = np.sum(w * I * ny * ny)
    return E, F, Pxx, Pyy


def make_exact(nx=256):
    dx = 2 * L / nx
    xs = -L + (np.arange(nx) + 0.5) * dx
    xs = xs[xs > 0]
    out = np.array([exact_point(x, 0.0) for x in xs])
    np.savez(os.path.join(HERE, 'hcyl_exact_%d.npz' % nx), x=xs, E=out[:, 0], F=out[:, 1], Pxx=out[:, 2], Pyy=out[:, 3])
    print('wrote', nx, out[[0, len(xs) // 4, len(xs) // 2, -1]])


def lastbin(d, v):
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.%s.*.bin' % v)) + glob.glob(os.path.join(d, '*.%s.*.bin' % v)))
    return fs[-1] if fs else None


def metrics(d):
    import bin_convert as bc
    r = bc.read_binary_as_athdf(lastbin(d, 'm1'))
    E, F1 = r['m1_e'][0], r['m1_f1'][0]
    x1, x2 = r['x1v'], r['x2v']
    nx = len(x1)
    ex = np.load(os.path.join(HERE, 'hcyl_exact_%d.npz' % nx))
    j0 = np.argmin(np.abs(x2)); j1 = j0 - 1 if x2[j0] > 0 else j0 + 1     # average the two rows by y = 0
    m = x1 > 0
    Ec = 0.5 * (E[j0] + E[j1])[m]
    Fc = 0.5 * (F1[j0] + F1[j1])[m]
    rv = bc.read_binary_as_athdf(lastbin(d, 'm1_vet')) if lastbin(d, 'm1_vet') else None
    out = {}
    Fc = Fc / float(open(os.path.join(d, 'in.athinput')).read().split('c_light')[1].split('=')[1].split()[0])
    xin, xout = ex['x'] < R0, ex['x'] > R0 + 0.05
    l1 = lambda a, b, s: np.sum(np.abs(a[s] - b[s])) / np.sum(np.abs(b[s]))
    out['E_L1_in'] = l1(Ec, ex['E'], xin)
    out['E_L1_out'] = l1(Ec, ex['E'], xout)
    out['F_L1_out'] = l1(Fc, ex['F'], xout)
    # r F at the box edge region (flux conservation r F = const for a cylinder)
    sel = (ex['x'] > 2.0) & (ex['x'] < 3.5)
    out['rF_code/rF_ex(2-3.5)'] = np.mean(Fc[sel] * ex['x'][sel]) / np.mean(ex['F'][sel] * ex['x'][sel])
    out['E(0)/ex'] = Ec[0] / ex['E'][0]
    out['E(R0)/ex'] = np.interp(R0, ex['x'], Ec) / np.interp(R0, ex['x'], ex['E'])
    if rv is not None and 'm1_d11' in rv:
        D = 0.5 * (rv['m1_d11'][0][j0] + rv['m1_d11'][0][j1])[m]
        fex = ex['Pxx'] / ex['E']
        out['fxx(0.5) code/ex'] = '%.4f/%.4f' % (np.interp(0.5, ex['x'], D), np.interp(0.5, ex['x'], fex))
        out['fxx(3) code/ex'] = '%.4f/%.4f' % (np.interp(3.0, ex['x'], D), np.interp(3.0, ex['x'], fex))
    return out, dict(x=ex['x'], E=Ec, F=Fc, Eex=ex['E'], Fex=ex['F'])


if __name__ == '__main__':
    if sys.argv[1] == 'exact':
        for n in sys.argv[2:] or ['256']:
            make_exact(int(n))
    else:
        for d in sys.argv[1:]:
            o, _ = metrics(d)
            print(os.path.basename(d.rstrip('/')), ' '.join('%s=%s' % (k, ('%.4g' % v) if isinstance(v, float) else v) for k, v in o.items()))
