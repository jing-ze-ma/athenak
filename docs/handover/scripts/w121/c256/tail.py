# ------------------------------------------------------------------ horizontal (dhj-remap-h)
# z components of PanelFrame (a, b, n) in src/coordinates/cubed_sphere.hpp, per panel
PANEL_Z = ((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
           (0.0, 1.0, 0.0), (0.0, 0.0, -1.0))


def corner_omega(X, Y):
    return np.arctan(X*Y/np.sqrt(1.0 + X*X + Y*Y))


def panel_cols(k0, k1, j0, j1, n2, n3, panel):
    """exact solid angle (GnomonicSolidAngle), cos_cell and sin^2(colatitude) of the panel
    columns [k0,k1) x [j0,j1) of a panel with n2 x n3 cells; arrays (k1-k0, j1-j0)"""
    xe = np.tan(0.25*np.pi*(-1.0 + 2.0*np.arange(j0, j1 + 1)/n2))
    ye = np.tan(0.25*np.pi*(-1.0 + 2.0*np.arange(k0, k1 + 1)/n3))
    om = (corner_omega(xe[None, 1:], ye[1:, None]) - corner_omega(xe[None, :-1], ye[1:, None])
          - corner_omega(xe[None, 1:], ye[:-1, None])
          + corner_omega(xe[None, :-1], ye[:-1, None]))
    X = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(j0, j1) + 0.5)/n2))[None, :]
    Y = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(k0, k1) + 0.5)/n3))[:, None]
    C, D = np.sqrt(1 + X*X), np.sqrt(1 + Y*Y)
    cosc = -X*Y/(C*D)
    a2, b2, n2z = PANEL_Z[panel]
    cz = (a2*X + b2*Y + n2z)/np.sqrt(1.0 + X*X + Y*Y)
    sin2 = 1.0 - cz*cz
    return om, cosc + 0*om, sin2 + 0*om


def panel_arrays(h, q, n2, n3):
    """per-panel array (6, n3+2, n2+2, O1) of a restart field q (nmb, O3, O2, O1): the
    active cells plus a one-cell face-ghost ring from the restart's own ghosts (at a panel
    seam these are the code's resampled cross-panel ghosts, vectors already in THIS panel's
    basis). Diagonal corners stay NaN: the dimension-by-dimension slopes never read them."""
    ng, nx2, nx3 = h['ng'], h['nx2'], h['nx3']
    P = np.full((6, n3 + 2, n2 + 2, q.shape[-1]), np.nan)
    a2, a3 = slice(ng, ng + nx2), slice(ng, ng + nx3)
    for m in range(h['nmb']):
        _, lx2, lx3, _, pn = h['lloc'][m]
        j0, k0 = lx2*nx2, lx3*nx3
        sj, sk = slice(1 + j0, 1 + j0 + nx2), slice(1 + k0, 1 + k0 + nx3)
        P[pn, sk, sj] = q[m, a3, a2]
        if j0 == 0:
            P[pn, sk, 0] = q[m, a3, ng - 1]
        if j0 + nx2 == n2:
            P[pn, sk, n2 + 1] = q[m, a3, ng + nx2]
        if k0 == 0:
            P[pn, 0, sj] = q[m, ng - 1, a2]
        if k0 + nx3 == n3:
            P[pn, n3 + 1, sj] = q[m, ng + nx3, a2]
    return P


def mc(dm, dp):
    """monotonized-central limited slope (per parent cell width)"""
    return np.where(dm*dp > 0.0,
                    np.sign(dm)*np.minimum(np.minimum(2*np.abs(dm), 2*np.abs(dp)),
                                           0.5*np.abs(dm + dp)), 0.0)


def slopes(P, pn, k0, k1, j0, j1):
    """cell value and MC slopes in j (x2) and k (x3) of panel parents [k0,k1) x [j0,j1)"""
    c = P[pn, 1 + k0:1 + k1, 1 + j0:1 + j1]
    sj = mc(c - P[pn, 1 + k0:1 + k1, j0:j1], P[pn, 1 + k0:1 + k1, 2 + j0:2 + j1] - c)
    sk = mc(c - P[pn, k0:k1, 1 + j0:1 + j1], P[pn, 2 + k0:2 + k1, 1 + j0:1 + j1] - c)
    return c, sj, sk


class HProlong:
    """conservative limited piecewise-linear prolongation of one panel patch by an
    integer factor f per direction.  In each parent the linear profile
    q + s_j (d_a - d_j0) + s_k (d_b - d_k0) is evaluated at the child offsets d_a =
    (a + 1/2)/f - 1/2 (parent units, the equiangular chart is uniform in xi, eta), with
    d_j0, d_k0 the EXACT-solid-angle weighted means of the child offsets, so sum_children
    q_c dOmega_c = q dOmega_parent to round-off (the volume is dOmega x the radial shell,
    so this holds per radial cell)."""

    def __init__(self, f, pk0, pk1, pj0, pj1, n2, n3, pn):
        self.f = f
        self.shape = (pk1 - pk0, pj1 - pj0)
        w, cosc, sin2 = panel_cols(pk0*f, pk1*f, pj0*f, pj1*f, n2*f, n3*f, pn)
        kp, jp = self.shape
        w4 = w.reshape(kp, f, jp, f)
        d = (np.arange(f) + 0.5)/f - 0.5
        W = w4.sum(axis=(1, 3))
        dj0 = (w4*d[None, None, None, :]).sum(axis=(1, 3))/W
        dk0 = (w4*d[None, :, None, None]).sum(axis=(1, 3))/W
        self.oj = d[None, None, None, :] - dj0[:, None, :, None]     # (kp, f, jp, f)
        self.ok = d[None, :, None, None] - dk0[:, None, :, None]
        self.w4 = w4
        self.cosc, self.sin2 = cosc, sin2          # of the children (kp f, jp f)

    def children(self, c, sj, sk, positive=False, floor=0.1):
        """children (kp f, jp f, O1) of parents c (kp, jp, O1) with slopes sj, sk.
        positive: both slopes scaled down in any parent whose smallest child would drop
        below floor x parent (conservation kept; returns the number of parents scaled)."""
        f = self.f
        kp, jp = self.shape
        oj = self.oj[..., None]
        ok = self.ok[..., None]
        sj4 = sj[:, None, :, None, :]
        sk4 = sk[:, None, :, None, :]
        nsc = 0
        if positive:
            dmin = (sj4*oj + sk4*ok).min(axis=(1, 3))       # (kp, jp, O1)
            lim = (floor - 1.0)*c
            bad = dmin < lim
            nsc = int(bad.sum())
            if nsc:
                th = np.where(bad, lim/np.where(bad, dmin, -1.0), 1.0)
                sj4 = sj4*th[:, None, :, None, :]
                sk4 = sk4*th[:, None, :, None, :]
        q = c[:, None, :, None, :] + sj4*oj + sk4*ok
        return q.reshape(kp*f, jp*f, -1), nsc

    def parent_sum(self, q):
        """sum_children q dOmega per parent: (kp f, jp f, O1) -> (kp, jp, O1)"""
        kp, jp = self.shape
        f = self.f
        return (q.reshape(kp, f, jp, f, -1)*self.w4[..., None]).sum(axis=(1, 3))


def zorder(lloc_list, level):
    """sort (lx1, lx2, lx3, level, panel) tuples as CreateZOrderedLLList numbers them:
    panel by panel, then the Morton key at the root level (x1 lowest bit)"""
    def key(ll):
        c = 0
        for b in range(level - 1, -1, -1):
            c = c*8 + (((ll[2] >> b) & 1) << 2) + (((ll[1] >> b) & 1) << 1) + ((ll[0] >> b) & 1)
        return (ll[4], c)
    return sorted(lloc_list, key=key)


def root_level(nb):
    lev = 0
    while (1 << lev) < max(nb):
        lev += 1
    return lev


def selftest():
    """gate (c): a smooth analytic field on one panel, exact cell averages (Gauss
    quadrature in dOmega) on coarse and fine grids; the prolonged coarse averages against
    the exact fine ones, L1 and Linf, for n = 8, 16, 32 per panel and f = 2, 8"""
    xg, wg = np.polynomial.legendre.leggauss(6)

    def fun(X, Y):
        return 2.0 + np.sin(1.3*np.arctan(X) + 0.4)*np.cos(0.9*np.arctan(Y) - 0.2) + 0.3*X*Y

    def avgs(n, lo, hi):
        """exact averages of cells lo..hi-1 (x and y) incl. a ring; returns (hi-lo)^2"""
        e = 0.25*np.pi*(-1.0 + 2.0*np.arange(lo, hi + 1)/n)
        A = np.empty((hi - lo, hi - lo))
        xi = 0.5*(e[1:] + e[:-1])[:, None] + 0.5*np.diff(e)[:, None]*xg[None, :]
        wx = 0.5*np.diff(e)[:, None]*wg[None, :]
        X = np.tan(xi)
        for kk in range(hi - lo):
            Y = X[kk][:, None, None]
            Xx = X[None, :, :]
            dens = (1 + Xx*Xx)*(1 + Y*Y)/(1 + Xx*Xx + Y*Y)**1.5
            ww = wx[kk][:, None, None]*wx[None, :, :]*dens
            A[kk] = (fun(Xx, Y)*ww).sum(axis=(0, 2))/ww.sum(axis=(0, 2))
        return A
    out = []
    for f in (2, 8):
        prev = None
        for n in (8, 16, 32):
            Pc = avgs(n, -1, n + 1)[None, :, :, None]           # ring = extended chart
            hp = HProlong(f, 0, n, 0, n, n, n, 0)
            c, sj, sk = slopes(Pc, 0, 0, n, 0, n)
            q, _ = hp.children(c, sj, sk)
            ex = avgs(n*f, 0, n*f)
            err = np.abs(q[..., 0] - ex)
            cons = np.abs(hp.parent_sum(q)[..., 0]/hp.w4.sum(axis=(1, 3)) - c[..., 0]).max()
            l1 = err.mean()
            out.append('  f %d n %3d: L1 %.3e  Linf %.3e  ratio L1 %s  cons %.1e'
                       % (f, n, l1, err.max(), '%.2f' % (prev/l1) if prev else '  - ', cons))
            prev = l1
    print('gate (c) smooth field, prolonged vs exact fine averages (order 2 = ratio 4):')
    print('\n'.join(out))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument('rst', nargs='?')
    ap_.add_argument('outdir', nargs='?')
    ap_.add_argument('--grid', default='')
    ap_.add_argument('--method', default='avg', choices=('avg', 'point'))
    ap_.add_argument('--interp-identity', action='store_true',
                     help='run the interpolation even when the grid is unchanged')
    ap_.add_argument('--henergy', default='total', choices=('total', 'eint'),
                     help='horizontal prolongation: conserve the total energy per parent '
                     '(e + KE + rho Phi, the etotgrav convention; e rescaled) or e only')
    ap_.add_argument('--hse-sample', type=int, default=16,
                     help='HSE statistics from every n-th new MeshBlock')
    ap_.add_argument('--check-only', action='store_true')
    ap_.add_argument('--selftest', action='store_true')
    a = ap_.parse_args()
    if a.selftest:
        selftest()
        return

    h = read_rst(a.rst)
    par = h['par']
    ng, nmb = h['ng'], h['nmb']
    nx1o, r0, r1, co = grid_of(par)
    if nx1o != h['nx1'] or int(par['meshblock']['nx1']) != nx1o:
        raise SystemExit('one MeshBlock in r is required')
    if 'cubed_sphere' in par['mesh'] and \
            par['mesh']['cubed_sphere'].lower() not in ('true', '1'):
        raise SystemExit('cubed-sphere restarts only')
    # target grid
    ov = dict(kv.split('=', 1) for kv in a.grid.split())
    parn = {b: dict(v) for b, v in par.items()}
    for k, v in ov.items():
        blk, key = k.split('/')
        parn.setdefault(blk, {})[key] = v
    nx1n, r0n, r1n, cn = grid_of(parn)
    if abs(r0n/r0 - 1) > 1e-12 or abs(r1n/r1 - 1) > 1e-12:
        raise SystemExit('x1min/x1max must not change (%g %g -> %g %g)'
                         % (r0, r1, r0n, r1n))
    # horizontal: per-panel cells and MeshBlock size, old -> new
    n2o, n3o = int(par['mesh']['nx2']), int(par['mesh']['nx3'])
    n2n, n3n = int(parn['mesh']['nx2']), int(parn['mesh']['nx3'])
    b2n, b3n = int(parn['meshblock']['nx2']), int(parn['meshblock']['nx3'])
    if n2n % n2o or n3n % n3o or n2n//n2o != n3n//n3o:
        raise SystemExit('horizontal: nx2, nx3 must grow by the same integer factor')
    if int(parn['meshblock']['nx1']) != nx1n:
        raise SystemExit('one MeshBlock in r is required (meshblock/nx1 = mesh/nx1)')
    if n2n % b2n or n3n % b3n:
        raise SystemExit('meshblock nx2/nx3 must divide mesh nx2/nx3')
    fac = n2n//n2o
    same_layout = (fac == 1 and b2n == h['nx2'] and b3n == h['nx3'])
    identity = (nx1n == nx1o and np.array_equal(cn, co) and same_layout)

    g = float(par['problem']['grav'])
    apl = float(par['problem']['ap'])
    etg = par['hydro'].get('etotgrav', 'false').lower() in ('true', '1')
    rotp = par['problem'].get('rot_potential', 'false').lower() in ('true', '1')
    om2 = float(par['problem'].get('omega', '0'))**2 if rotp else 0.0
    if par['problem'].get('stellar_tide', 'false').lower() in ('true', '1') and etg:
        raise SystemExit('stellar_tide with etotgrav: tidal potential not handled')

    eo = edges(nx1o, r0, r1, co, ng)            # nx1o + 2 ng + 1 edges
    ro = centroids(eo)
    en = edges(nx1n, r0, r1, cn, 0)
    rn = centroids(en)
    dvo = (eo[1:]**3 - eo[:-1]**3)/3.0
    dvn = (en[1:]**3 - en[:-1]**3)/3.0
    s3 = slice(ng, -ng)
    act = slice(ng, ng + nx1o)
    phig_o = g*apl*(1.0 - apl/ro)
    phig_n = g*apl*(1.0 - apl/rn)
    if h['eint'] is None:
        raise SystemExit('restart without the EINTRST1 slab (not a general-EOS run)')

    # ---------------------------------------------------------- old state checks
    u = h['u'][:, :, s3, s3, :]                 # (nmb, 5, nx3, nx2, O1)
    rho, m1, m2, m3, E = (u[:, n] for n in range(5))
    ei = h['eint'][:, s3, s3, :]
    geo = [panel_cols(ll[2]*h['nx3'], (ll[2] + 1)*h['nx3'], ll[1]*h['nx2'],
                      (ll[1] + 1)*h['nx2'], n2o, n3o, ll[4]) for ll in h['lloc']]
    dom = np.stack([gg[0] for gg in geo])
    cosc = np.stack([gg[1] for gg in geo])
    Bana = -0.5*om2*np.stack([gg[2] for gg in geo])
    ke = kinetic(rho, m1, m2, m3, cosc[..., None])
    if etg:
        phi_c = (E - ke - ei)/rho
        w = ro[act]**2
        B = ((phi_c[..., act] - phig_o[act])*w).sum(-1)/(w*w).sum()
        resid = phi_c[..., act] - phig_o[act] - B[..., None]*w
        print('Phi fit: max |Phi_code - (g ap (1-ap/r) + B r^2)| / |Phi| = %.2e; '
              'max |B_fit - B_analytic| / (Omega^2/2) = %.2e (rot_potential %s)'
              % (np.abs(resid/phi_c[..., act]).max(),
                 np.abs(B - Bana).max()/max(0.5*om2, 1e-300), rotp))
    else:
        B = np.zeros(rho.shape[:-1])
    if not same_layout:
        B = Bana                                  # the analytic B of TotPotAt
    phio = phig_o[None, None, None, :] + B[..., None]*ro**2
    rel_c2p = np.abs(E - ke - ei - (rho*phio if etg else 0.0))[..., act] / ei[..., act]
    print('state check: max |E - KE - e - rho Phi| / e (active) = %.2e' % rel_c2p.max())

    # ---------------------------------------------------------- budgets
    KEYS = ('mass', 'mom1', 'mom2', 'mom3', 'eint', 'ekin', 'egrav', 'etot')

    def budget(r_, a1, a2, a3, e_, c_, B_, dv, rc, pg, dom_):
        k_ = kinetic(r_, a1, a2, a3, c_[..., None])
        eg = r_*(pg + B_[..., None]*rc**2) if etg else 0.0*r_
        cols = dict(mass=r_, mom1=a1, mom2=a2, mom3=a3, eint=e_, ekin=k_, egrav=eg,
                    etot=e_ + k_ + eg)
        return {k: ((v*dv).sum(-1)*dom_).sum() for k, v in cols.items()}
    bo = budget(rho[..., act], m1[..., act], m2[..., act], m3[..., act], ei[..., act],
                cosc, B, dvo[act], ro[act], phig_o[act], dom)
    bh = dict.fromkeys(KEYS, 0.0)
    bn = dict.fromkeys(KEYS, 0.0)
    colmax = dict(mass=0.0, eint=0.0, mom1=0.0)

    # ---------------------------------------------------------- new layout
    nb2, nb3 = n2n//b2n, n3n//b3n
    if same_layout:
        lloc_new = [tuple(int(x) for x in ll) for ll in h['lloc']]
    else:
        lev = root_level((1, nb2, nb3))
        lloc_new = zorder([(0, j, k, lev, p) for p in range(6) for k in range(nb3)
                           for j in range(nb2)], lev)
        # the ordering rule reproduces the restart's own list
        lev_o = root_level((1, n2o//h['nx2'], n3o//h['nx3']))
        chk = zorder([tuple(int(x) for x in ll) for ll in h['lloc']], lev_o)
        if lev_o != h['lloc'][0, 3] or \
                chk != [tuple(int(x) for x in ll) for ll in h['lloc']]:
            raise SystemExit('Z-order / root level rule does not reproduce the restart')
    nmbn = len(lloc_new)
    print('horizontal: C%d -> C%d (factor %d), MeshBlocks %d x %d x %d -> %d x %d x %d '
          '(%d -> %d); henergy %s'
          % (n2o, n2n, fac, nx1o, h['nx2'], h['nx3'], nx1n, b2n, b3n, nmb, nmbn,
             a.henergy))

    if fac > 1 or not same_layout:
        names = ('rho', 'm1', 'm2', 'm3', 'ei', 'p', 'E')
        full = dict(rho=h['u'][:, 0], m1=h['u'][:, 1], m2=h['u'][:, 2], m3=h['u'][:, 3],
                    ei=h['eint'], p=h['p'], E=h['u'][:, 4])
        PA = {k: panel_arrays(h, full[k], n2o, n3o) for k in names}
        nsc_tot = dict(rho=0, ei=0, p=0)
        fstat = [np.inf, -np.inf, 0]

    method = a.method
    if identity and not a.interp_identity:
        method = 'copy'
    dxo, dxn = np.diff(eo), np.diff(en)

    def remap_ln(q):
        L = np.log(q)
        if method == 'avg':
            L = L - avg_corr(L, ro, dxo)
        Ln = PchipInterpolator(ro, L, axis=-1)(rn)
        if method == 'avg':
            Ln = Ln + avg_corr(Ln, rn, dxn)
        return np.exp(Ln)

    fmr, fer = [np.inf, -np.inf], [np.inf, -np.inf]
    dfl = float(par['hydro'].get('dfloor', '0'))
    ncl = 0
    hse_o, hse_n = [], []
    rmin, emin, pmin, pmax = np.inf, np.inf, np.inf, -np.inf

    def hse(rr, pp, r, B_):
        rf = 0.5*(r[1:] + r[:-1])
        geff = g*apl**2/rf**2 + 2.0*B_[..., None]*rf
        dpdr = np.diff(pp, axis=-1)/np.diff(r)
        rhof = 0.5*(rr[..., 1:] + rr[..., :-1])
        return dpdr/(-rhof*geff) - 1.0, np.sqrt(pp[..., 1:]*pp[..., :-1])
    for m in range(nmb):
        if m % a.hse_sample == 0:
            hse_o.append(hse(rho[m][..., act], h['p'][m, s3, s3, act], ro[act], B[m]))

    fo = None
    if not a.check_only:
        os.makedirs(a.outdir, exist_ok=True)
        fn = os.path.join(a.outdir, 'remap.dat')
        fo = open(fn, 'wb')
        fo.write(b'RADREMPE')
        fo.write(struct.pack('<4i', nmbn, nx1n, b2n, b3n))

    for mn, ll in enumerate(lloc_new):
        _, lx2, lx3, _, pn = ll
        # ---- horizontal: this block's columns on the OLD radial grid (O1, ghosts incl.)
        if same_layout and fac == 1:
            c_rho, c_m1, c_m2, c_m3 = rho[mn], m1[mn], m2[mn], m3[mn]
            c_ei, c_p = ei[mn], h['p'][mn, s3, s3, :]
            c_cos, c_B, c_dom = cosc[mn], B[mn], dom[mn]
        else:
            J0, K0 = lx2*b2n, lx3*b3n
            pj0, pj1 = J0//fac, -(-(J0 + b2n)//fac)
            pk0, pk1 = K0//fac, -(-(K0 + b3n)//fac)
            hp = HProlong(fac, pk0, pk1, pj0, pj1, n2o, n3o, pn)
            sl = (slice(K0 - pk0*fac, K0 - pk0*fac + b3n),
                  slice(J0 - pj0*fac, J0 - pj0*fac + b2n))
            ch = {}
            for k in ('rho', 'm1', 'm2', 'm3', 'ei', 'p'):
                c, sj, sk = slopes(PA[k], pn, pk0, pk1, pj0, pj1)
                ch[k], nsc = hp.children(c, sj, sk, positive=(k in nsc_tot))
                if k in nsc_tot:
                    nsc_tot[k] += nsc
            B_f = -0.5*om2*hp.sin2
            if a.henergy == 'total' and fac > 1:
                # e x f per parent so that sum_children (e + KE + rho Phi) dOmega is the
                # parent's E dOmega (active radial cells; radial ghosts keep f = 1)
                kef = kinetic(ch['rho'], ch['m1'], ch['m2'], ch['m3'], hp.cosc[..., None])
                rpf = ch['rho']*(phig_o + B_f[..., None]*ro**2) if etg else 0.0*kef
                Ep = PA['E'][pn, 1 + pk0:1 + pk1, 1 + pj0:1 + pj1]
                Wp = hp.w4.sum(axis=(1, 3))[..., None]
                ff = (Ep*Wp - hp.parent_sum(kef + rpf))/hp.parent_sum(ch['ei'])
                ff[..., :ng] = 1.0
                ff[..., ng + nx1o:] = 1.0
                fstat[0] = min(fstat[0], ff.min())
                fstat[1] = max(fstat[1], ff.max())
                bad = ff < 0.5
                fstat[2] += int(bad.sum())
                ff = np.where(bad, 0.5, ff)
                kp, jp = hp.shape
                ch['ei'] = (ch['ei'].reshape(kp, fac, jp, fac, -1) *
                            ff[:, None, :, None, :]).reshape(kp*fac, jp*fac, -1)
            c_rho, c_m1, c_m2, c_m3 = (ch[k][sl] for k in ('rho', 'm1', 'm2', 'm3'))
            c_ei, c_p = ch['ei'][sl], ch['p'][sl]
            c_cos, c_B = hp.cosc[sl], B_f[sl]
            c_dom = panel_cols(K0, K0 + b3n, J0, J0 + b2n, n2n, n3n, pn)[0]
        bb = budget(c_rho[..., act], c_m1[..., act], c_m2[..., act], c_m3[..., act],
                    c_ei[..., act], c_cos, c_B, dvo[act], ro[act], phig_o[act], c_dom)
        for k in KEYS:
            bh[k] += bb[k]
        # ---- radial, per column (the dhj-remap method, unchanged)
        if method == 'copy':
            rn_, m1n, m2n, m3n = (q[..., act] for q in (c_rho, c_m1, c_m2, c_m3))
            ein, pn_ = c_ei[..., act], c_p[..., act]
        else:
            rr = remap_ln(c_rho)
            ee = remap_ln(c_ei)
            pp = remap_ln(c_p)
            f_m = (c_rho[..., act]*dvo[act]).sum(-1)/(rr*dvn).sum(-1)
            f_e = (c_ei[..., act]*dvo[act]).sum(-1)/(ee*dvn).sum(-1)
            fmr = [min(fmr[0], f_m.min()), max(fmr[1], f_m.max())]
            fer = [min(fer[0], f_e.min()), max(fer[1], f_e.max())]
            rn_ = f_m[..., None]*rr
            ein = f_e[..., None]*ee
            pn_ = f_e[..., None]*pp
            mom = []
            for q in (c_m1, c_m2, c_m3):
                un = PchipInterpolator(ro, q/c_rho, axis=-1)(rn)
                du = ((q[..., act]*dvo[act]).sum(-1) - (rn_*un*dvn).sum(-1)) / \
                    (rn_*dvn).sum(-1)
                mom.append(rn_*(un + du[..., None]))
            m1n, m2n, m3n = mom
        nc = int((rn_ < dfl).sum())
        if nc:
            s_ = np.maximum(dfl/rn_, 1.0)
            rn_, m1n, m2n, m3n = rn_*s_, m1n*s_, m2n*s_, m3n*s_
            ncl += nc
        bb = budget(rn_, m1n, m2n, m3n, ein, c_cos, c_B, dvn, rn, phig_n, c_dom)
        for k in KEYS:
            bn[k] += bb[k]
        for k, q_o, q_n in (('mass', c_rho, rn_), ('eint', c_ei, ein), ('mom1', c_m1, m1n)):
            co_ = (q_o[..., act]*dvo[act]).sum(-1)
            cn_ = (q_n*dvn).sum(-1)
            ref = np.abs(co_) if k != 'mom1' else (c_rho[..., act]*dvo[act]).sum(-1)*1e5
            colmax[k] = max(colmax[k], np.abs((cn_ - co_)/ref).max())
        if mn % a.hse_sample == 0:
            hse_n.append(hse(rn_, pn_, rn, c_B))
        rmin, emin = min(rmin, rn_.min()), min(emin, ein.min())
        pmin, pmax = min(pmin, pn_.min()), max(pmax, pn_.max())
        if not (np.isfinite(rn_).all() and np.isfinite(ein).all() and
                np.isfinite(m1n).all() and np.isfinite(m2n).all() and
                np.isfinite(m3n).all()):
            raise SystemExit('non-finite value in new MeshBlock %d %s' % (mn, ll))
        if fo is not None:
            fo.write(np.array(ll, '<i4').tobytes())
            rec = np.stack([rn_, m1n/rn_, m2n/rn_, m3n/rn_, ein])
            fo.write(np.ascontiguousarray(rec, '<f8').tobytes())
    if fo is not None:
        fo.close()

    # ---------------------------------------------------------- report
    print('rst %s: t = %.6e s, cycle %d, sections %s'
          % (a.rst, h['t'], h['ncycle'], h['sections']))
    print('grid nx1 %d -> %d (%s); method %s; dr min/max new %.3e / %.3e cm'
          % (nx1o, nx1n, 'identity' if identity else 'remap', method,
             np.diff(en).min(), np.diff(en).max()))
    if fac > 1 or not same_layout:
        print('horizontal: parents with slopes scaled for positivity: rho %d, e %d, p %d; '
              'energy factor f (e x f per parent) %.6f..%.6f, clipped at 0.5: %d'
              % (nsc_tot['rho'], nsc_tot['ei'], nsc_tot['p'],
                 fstat[0] if fstat[0] < np.inf else 1.0,
                 fstat[1] if fstat[1] > -np.inf else 1.0, fstat[2]))
    if method != 'copy':
        print('column factors: f_m %.6f..%.6f, f_e %.6f..%.6f (1 = the interpolation '
              'alone conserved it)' % (fmr[0], fmr[1], fer[0], fer[1]))
    print('cells raised to dfloor = %g: %d' % (dfl, ncl))
    print('%-7s %16s %16s %16s %10s %10s' % ('budget', 'old', 'horizontal', 'new',
                                             'h/old-1', 'new/h-1'))
    for k in KEYS:
        print('%-7s %16.9e %16.9e %16.9e %+10.2e %+10.2e'
              % (k, bo[k], bh[k], bn[k], bh[k]/bo[k] - 1 if bo[k] else 0.0,
                 bn[k]/bh[k] - 1 if bh[k] else 0.0))
    print('  (momenta: sums of the stored covariant components; relative to |mom| they '
          'are near-cancelling sums: vs mass x 1 km/s: h %.1e %.1e %.1e)'
          % tuple(abs(bh[k] - bo[k])/(abs(bo['mass'])*1e5) for k in
                  ('mom1', 'mom2', 'mom3')))
    print('  radial step, max |column change| (mass, eint rel; mom1 vs column mass x 1 '
          'km/s): %.2e %.2e %.2e' % (colmax['mass'], colmax['eint'], colmax['mom1']))
    print('HSE residual |dp/dr / (-rho g_eff) - 1| at faces, median / p90 (old -> new, '
          'every %d-th MeshBlock)' % a.hse_sample)
    ho = np.concatenate([x[0].ravel() for x in hse_o])
    po = np.concatenate([x[1].ravel() for x in hse_o])
    hn = np.concatenate([x[0].ravel() for x in hse_n])
    pnf = np.concatenate([x[1].ravel() for x in hse_n])
    for lo, hi in ((1e-9, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 1.0), (1.0, 10.0),
                   (10.0, 100.0), (100.0, 1e4)):
        x = np.abs(ho[(po >= lo*BAR) & (po < hi*BAR)])
        y = np.abs(hn[(pnf >= lo*BAR) & (pnf < hi*BAR)])
        if x.size and y.size:
            print('  %8.0e-%-8.0e bar  %8.2e / %8.2e  ->  %8.2e / %8.2e' % (
                lo, hi, np.median(x), np.percentile(x, 90), np.median(y),
                np.percentile(y, 90)))
    print('rho min %.3e -> %.3e; e min %.3e -> %.3e; p range new %.3e..%.3e bar'
          % (rho[..., act].min(), rmin, ei[..., act].min(), emin, pmin/BAR, pmax/BAR))
    if a.check_only:
        return
    txt = h['text']
    txt = txt[:txt.find('<par_end>')]
    for k, v in ov.items():
        blk, key = k.split('/')
        txt = set_key(txt, blk, key, v, 'remap target grid')
    txt = set_key(txt, 'time', 'start_time', repr(h['t']),
                  'remap: t of %s' % os.path.basename(a.rst))
    txt = set_key(txt, 'problem', 'remap_file', os.path.abspath(fn), 'remap')
    with open(os.path.join(a.outdir, 'remap.athinput'), 'w') as fo:
        fo.write(txt)
    print('wrote', fn, 'and', os.path.join(a.outdir, 'remap.athinput'))


if __name__ == '__main__':
    sys.exit(main())
