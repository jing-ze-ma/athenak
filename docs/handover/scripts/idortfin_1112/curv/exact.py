"""tanacc-1107 exact moments for a uniformly bright sphere (r = 1, F_in = 1, I_s = 4 in E units,
c = 1) in vacuum plus a uniform spherical clump (centre xc, radius a, extinction chi, source B):
  I(n) = I_s [ray hits star] exp(-tau(n)) + B (1 - exp(-tau(n))),
tau(n) = chi x (length of the backward ray x - s n, 0 < s < s_star, inside the clump).
Moments J = <I>, H = <I n>, K = <I n n> (<.> = int dOmega/4pi), by product Gauss-Legendre (mu) x
uniform (phi) quadrature on the star cap and on the clump cap (full sphere inside the clump)."""
import numpy as np

IS = 4.0


def _frame(a):
    a = a/np.linalg.norm(a, axis=-1, keepdims=True)
    t = np.where(np.abs(a[..., 2:3]) < 0.9, np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]))
    e1 = np.cross(a, t)
    e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(a, e1)
    return a, e1, e2


def _cap_dirs(axis, cmin, nmu, nph):
    """directions (P, Q, 3) and weights (P, Q) (sum = cap solid angle / 4 pi) on caps
    {n . axis >= cmin} (cmin per point)."""
    x, w = np.polynomial.legendre.leggauss(nmu)
    ph = (np.arange(nph) + 0.5)*2*np.pi/nph
    mu = cmin[:, None] + (1.0 - cmin[:, None])*0.5*(x[None, :] + 1.0)        # (P, nmu)
    wm = (1.0 - cmin[:, None])*0.5*w[None, :]                                # (P, nmu)
    s = np.sqrt(np.clip(1.0 - mu**2, 0, None))
    a, e1, e2 = _frame(axis)
    cp, sp = np.cos(ph), np.sin(ph)
    n = (mu[:, :, None, None]*a[:, None, None, :]
         + s[:, :, None, None]*(cp[None, None, :, None]*e1[:, None, None, :]
                                + sp[None, None, :, None]*e2[:, None, None, :]))
    wt = np.repeat(wm[:, :, None], nph, axis=2)/(2.0*nph)                    # dmu dphi/4pi
    P = axis.shape[0]
    return n.reshape(P, -1, 3), wt.reshape(P, -1)


def _intensity(x, n, clump, star_mask_only=False):
    """x (P,3), n (P,Q,3): the exact I(n)."""
    # backward ray y(s) = x - s n; star hit: |x - s n|^2 = 1, s > 0
    xn = np.einsum('pi,pqi->pq', x, n)
    r2 = np.sum(x*x, axis=1)[:, None]
    disc = xn**2 - (r2 - 1.0)
    hit = (disc >= 0) & (xn > 0)
    sstar = np.where(hit, xn - np.sqrt(np.clip(disc, 0, None)), np.inf)
    ib = np.where(hit, IS, 0.0)
    if clump is None:
        return ib
    xc, a, chi, B = clump
    d = x - xc[None, :]
    dn = np.einsum('pi,pqi->pq', d, n)
    d2 = np.sum(d*d, axis=1)[:, None]
    dc = dn**2 - (d2 - a*a)
    ok = dc > 0
    sq = np.sqrt(np.clip(dc, 0, None))
    s1 = np.clip(dn - sq, 0, None)
    s2 = np.minimum(dn + sq, sstar)
    L = np.where(ok, np.clip(s2 - s1, 0, None), 0.0)
    tau = chi*L
    et = np.exp(-tau)
    return ib*et + B*(1.0 - et)


def moments(x, clump=None, nmu=96, nph=160, chunk=64):
    """x (N,3) points (r > 1).  clump = (xc(3), a, chi, B) or None.  Returns J (N), H (N,3),
    K (N,3,3)."""
    N = x.shape[0]
    J = np.zeros(N)
    H = np.zeros((N, 3))
    K = np.zeros((N, 3, 3))
    for c0 in range(0, N, chunk):
        xs = x[c0:c0 + chunk]
        r = np.linalg.norm(xs, axis=1)
        caps = [(xs/r[:, None], np.sqrt(np.clip(1 - 1/r**2, 0, None)))]   # star cap
        if clump is not None and clump[3] > 0:
            xc, a = clump[0], clump[1]
            d = xs - xc[None, :]
            dd = np.linalg.norm(d, axis=1)
            d = np.where(dd[:, None] > 1e-12, d, np.array([1.0, 0.0, 0.0]))
            dd = np.maximum(dd, 1e-12)
            cm = np.where(dd > a*1.0000001, np.sqrt(np.clip(1 - (a/dd)**2, 0, None)), -1.0)
            caps.append((d/dd[:, None], cm))
        tot_I = []
        for ci, (ax, cmin) in enumerate(caps):
            n, w = _cap_dirs(ax, cmin, nmu, nph)
            I = _intensity(xs, n, clump)
            if ci == 0:
                if clump is not None and clump[3] > 0:
                    I = I - _emis_only(xs, n, clump)       # clump emission taken on its cap
            else:
                I = _emis_only(xs, n, clump)
                # remove the star term double count: only the emission on the clump cap
            wi = w*I
            J[c0:c0 + chunk] += wi.sum(1)
            H[c0:c0 + chunk] += np.einsum('pq,pqi->pi', wi, n)
            K[c0:c0 + chunk] += np.einsum('pq,pqi,pqj->pij', wi, n, n)
    return J, H, K


def _emis_only(x, n, clump):
    xc, a, chi, B = clump
    xn = np.einsum('pi,pqi->pq', x, n)
    r2 = np.sum(x*x, axis=1)[:, None]
    disc = xn**2 - (r2 - 1.0)
    hit = (disc >= 0) & (xn > 0)
    sstar = np.where(hit, xn - np.sqrt(np.clip(disc, 0, None)), np.inf)
    d = x - xc[None, :]
    dn = np.einsum('pi,pqi->pq', d, n)
    d2 = np.sum(d*d, axis=1)[:, None]
    dc = dn**2 - (d2 - a*a)
    ok = dc > 0
    sq = np.sqrt(np.clip(dc, 0, None))
    s1 = np.clip(dn - sq, 0, None)
    s2 = np.minimum(dn + sq, sstar)
    L = np.where(ok, np.clip(s2 - s1, 0, None), 0.0)
    return B*(1.0 - np.exp(-chi*L))


if __name__ == '__main__':
    # checks: star only, E = 2 (1 - sqrt(1 - 1/r^2)), F r^2 = 1, K_rr/J exact
    rr = np.array([1.05, 2.0, 5.0])
    x = np.stack([rr, 0*rr, 0*rr], 1)
    J, H, K = moments(x)
    print('star J', J, 'exact', 2*(1 - np.sqrt(1 - 1/rr**2)), 'H r^2', H[:, 0]*rr**2)
    mu0 = np.sqrt(1 - 1/rr**2)
    print('D_rr', K[:, 0, 0]/J, 'exact', (1 + mu0 + mu0**2)/3)
    # clump checks: emission-only point far away: J -> B Omega_c/4pi (1 - e^-tau) approx
    cl = (np.array([3.0, 0, 0]), 0.4, 1e3, 0.0)
    for nm in (48, 96, 192):
        J2, H2, K2 = moments(np.array([[5.0, 0.0, 0.0], [5.0, 0.3, 0.0], [5.0, 0.6, 0.0]]),
                             cl, nmu=nm, nph=int(nm*5/3))
        print('shadow nmu', nm, 'J', J2, 'Kyy-Kzz', K2[:, 1, 1] - K2[:, 2, 2])
    cl = (np.array([3.0, 0, 0]), 0.4, 50.0, 0.1)
    J3, H3, K3 = moments(np.array([[3.0, 0.0, 0.0], [3.3, 0.2, 0.0], [5.0, 0.0, 0.0]]), cl)
    print('hot J (inside, edge, far)', J3)
