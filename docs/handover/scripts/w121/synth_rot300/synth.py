"""synth.py: WASP-121b w1x rot-300 synthetic observables (synth_rot300).

  python synth.py phase  [olr files...]   phase curves, T_b, energy budget   -> phase.txt, PNGs
  python synth.py winds                   limb line-of-sight velocities      -> winds.txt

Geometry: dhjcs panel frames (substellar = -x, lon = atan2(-Y,-X), east = +lon = the
rotation direction).  Column solid angle from the gnomonic area weight (as ana300.py).
Phase-curve projection: FLUX-BASED, Lambertian: each column emits its net top flux F
isotropically (I = F/pi), so the disk flux is r_top^2 sum F mu dOmega / pi (no limb
darkening beyond the cos weighting; no angle information survives in rt_Fb).
"""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhjcs   # noqa: E402

O = '/viper/ptmp2/jinma/w121prod_0929/synth_rot300/'
BIN = '/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.%05d.bin'
# constants (cgs)
h, c, kB, sig = 6.62607015e-27, 2.99792458e10, 1.380649e-16, 5.670374419e-5
RJ, RSUN = 7.1492e9, 6.957e10
RP = 1.742*RJ              # Sing+2024 (the input's R_p)
RS = 1.461*RSUN            # Sing+2024 / input
TS = 6628.0                # Sing+2024 T_eff; stellar spectrum = blackbody
AOR = 3.7844               # a/R*, Sing+2024
PROT = 1.27492504*86400.0
OMEGA = 5.704026e-05       # input problem/omega
RTOP = 1.631037e10         # x1max = the face ie+1 where rt_Fb is taken
R0, R1 = 1.126575e10, 1.631037e10
CSTR = [-0.379629, -1.797216, 7.873416, -79.711522, 347.885804, -753.943230,
        784.221097, -314.370119]
INSTR = {'NRS1': (2.70, 3.72), 'NRS2': (3.82, 5.15), 'SOSS1': (0.85, 2.85),
         'IRAC1': (3.13, 3.96), 'IRAC2': (3.92, 5.06), 'TESS': (0.6, 1.0),
         'WFC3': (1.12, 1.64)}


def planck_int(l0, l1, T, n=400):
    """int_{l0}^{l1} B_lambda(T) dlambda, lambda in um; erg/s/cm2/sr"""
    lam = np.linspace(l0, l1, n)*1e-4
    T = np.atleast_1d(np.asarray(T, float))
    x = h*c/(lam[None, :]*kB*T[:, None])
    B = 2*h*c*c/lam[None, :]**5/np.expm1(np.minimum(x, 700))
    return np.trapezoid(B, lam, axis=1)


def tb_invert(target, l0, l1):
    """T with planck_int(l0,l1,T) = target (vectorised bisection in log T)"""
    target = np.atleast_1d(target)
    lo, hi = np.full(target.shape, 100.0), np.full(target.shape, 2e4)
    for _ in range(60):
        mid = np.sqrt(lo*hi)
        f = planck_int(l0, l1, mid) - target
        lo = np.where(f < 0, mid, lo)
        hi = np.where(f >= 0, mid, hi)
    return np.sqrt(lo*hi)


def geometry(raw):
    g = np.asarray(raw['mb_geometry'])
    nmb, nj, nk = g.shape[0], raw['nx2_mb'], raw['nx3_mb']
    per = nmb//6
    X = np.zeros((nmb, nk, nj, 3))
    EXI = np.zeros_like(X)
    EET = np.zeros_like(X)
    w = np.zeros((nmb, nk, nj))
    for m in range(nmb):
        p = m//per
        xi = np.pi/4*(g[m, 2] + (g[m, 3] - g[m, 2])*(np.arange(nj) + 0.5)/nj)
        eta = np.pi/4*(g[m, 4] + (g[m, 5] - g[m, 4])*(np.arange(nk) + 0.5)/nk)
        XI, ETA = np.meshgrid(xi, eta)
        X[m] = dhjcs.panel_map(p, XI, ETA)
        d = 1e-6
        e = dhjcs.panel_map(p, XI + d, ETA) - dhjcs.panel_map(p, XI - d, ETA)
        EXI[m] = e/np.linalg.norm(e, axis=-1)[..., None]
        e = dhjcs.panel_map(p, XI, ETA + d) - dhjcs.panel_map(p, XI, ETA - d)
        EET[m] = e/np.linalg.norm(e, axis=-1)[..., None]
        x, y = np.tan(XI), np.tan(ETA)
        w[m] = (1 + x*x)*(1 + y*y)/(1 + x*x + y*y)**1.5
    dOm = 4*np.pi*w/w.sum()
    return X, EXI, EET, dOm


def read_olr(fn):
    edges = None
    with open(fn) as f:
        for ln in f:
            if ln.startswith('# band edges'):
                edges = np.array(ln.split(':')[1].split(), float)
            if ln.startswith('# t ='):
                t = float(ln.split()[3])
    d = np.loadtxt(fn)
    return t, edges, d


def phase():
    raw = dhjcs.bin_convert.read_binary(BIN % 150)
    X, _, _, dOm = geometry(raw)
    nmb, nk, nj = dOm.shape
    files = sys.argv[2:]
    phi = np.linspace(0, 1, 721)
    lo = np.pi - 2*np.pi*phi                           # sub-observer longitude
    obs = np.stack([-np.cos(lo), -np.sin(lo), 0*lo], -1)
    mu = np.einsum('mkjx,px->pmkj', X, obs)
    mu = np.maximum(mu, 0.0)
    out = []
    curves = {}
    for fn in files:
        t, edges, d = read_olr(fn)
        gid, kk, jj = d[:, 0].astype(int), d[:, 1].astype(int), d[:, 2].astype(int)
        nb = d.shape[1] - 5
        F = np.zeros((nb, nmb, nk, nj))
        mu0 = np.zeros((nmb, nk, nj))
        fsw = np.zeros((nmb, nk, nj))
        for b in range(nb):
            F[b, gid, kk, jj] = d[:, 5 + b]
        mu0[gid, kk, jj] = d[:, 3]
        fsw[gid, kk, jj] = d[:, 4]
        # consistency: stellar mu0 of the dump vs -X_x of the geometry (maps gid -> bin order)
        chk = np.corrcoef(mu0.ravel(), -X[..., 0].ravel())[0, 1]
        lam = np.sort(edges)                              # ascending um
        blo, bhi = (edges[1:], edges[:-1]) if edges[0] > edges[-1] else (edges[:-1], edges[1:])
        Lb = RTOP**2*np.einsum('bmkj,mkj->b', F, dOm)
        Lsw = RTOP**2*(fsw*dOm).sum()
        # disk flux x d^2, per band, per phase: r_top^2 sum F mu dOm / pi
        P = RTOP**2*np.einsum('bmkj,pmkj,mkj->pb', F, mu, dOm)/np.pi
        Pbol = P.sum(1)
        star_b = np.array([np.pi*planck_int(a, b_, TS)[0] for a, b_ in zip(blo, bhi)])*RS**2
        fpfs = {'bol': Pbol/(RS**2*sig*TS**4)}
        for b in range(nb):
            fpfs['ck%d_%.2f-%.2f' % (b, blo[b], bhi[b])] = P[:, b]/star_b[b]
        # planet brightness T per ck band (disk flux / (pi R_p^2) = pi int B)
        Tbb = np.zeros_like(P)
        for b in range(nb):
            Tbb[:, b] = tb_invert(P[:, b]/(np.pi*RP**2)/np.pi, blo[b], bhi[b])
        # instrument bands: piecewise-Planck within each ck band at its own T_b
        for nm, (l0, l1) in INSTR.items():
            Pi = np.zeros(len(phi))
            for b in range(nb):
                a, z = max(l0, blo[b]), min(l1, bhi[b])
                if z <= a:
                    continue
                Pi += P[:, b]*planck_int(a, z, Tbb[:, b])/planck_int(blo[b], bhi[b], Tbb[:, b])
            fpfs[nm] = Pi/(np.pi*planck_int(l0, l1, TS)[0]*RS**2)
        curves[fn] = (t, fpfs)
        out.append('# %s  t=%.6e s (rot %.2f)  gid-map corr(mu0,-X)=%.5f  L_IR=%.4e  '
                   'L_sw_abs=%.4e erg/s' % (fn.split('/')[-1], t, t/PROT, chk, Lb.sum(), Lsw))
        out.append('#   band L_IR fractions: ' + ' '.join('%.2f-%.2f:%.3f' % (blo[b], bhi[b],
                   Lb[b]/Lb.sum()) for b in range(nb)))
        # day / night effective T from the disk flux (Lambertian: P = R^2 sigma T^4)
        i0, i5 = 0, np.argmin(abs(phi - 0.5))
        for R, nmR in ((RP, 'R_p'), (RTOP, 'r_top')):
            Td = (Pbol[i5]/(R**2*sig))**0.25
            Tn = (Pbol[i0]/(R**2*sig))**0.25
            out.append('#   T_day,eff %.1f  T_night,eff %.1f  (disk flux normalised at %s)'
                       % (Td, Tn, nmR))
        # hemispheric (unweighted) means of the bolometric OLR at the top, scaled to R_p
        Ft = F.sum(0)*RTOP**2/RP**2
        day = (-X[..., 0]) > 0
        for nmh, msk in (('day', day), ('night', ~day)):
            out.append('#   hemisphere-mean OLR (at R_p) %s: T = %.1f K' % (
                nmh, ((Ft*dOm*msk).sum()/(dOm*msk).sum()/sig)**0.25))
        for nm, y in fpfs.items():
            # peak from a 2-harmonic fit
            A = np.stack([np.ones_like(phi), np.cos(2*np.pi*phi), np.sin(2*np.pi*phi),
                          np.cos(4*np.pi*phi), np.sin(4*np.pi*phi)], 1)
            cfit = np.linalg.lstsq(A, y, rcond=None)[0]
            yf = A @ cfit
            pk = phi[np.argmax(yf)]
            mn = phi[np.argmin(yf)]
            row = (nm, y[i5]*1e6, y[i0]*1e6, (y.max() - y.min())*1e6, (0.5 - pk)*360.0,
                   (mn - 0.0 if mn < 0.5 else mn - 1.0)*360.0)
            tb = ''
            if nm in INSTR:
                l0, l1 = INSTR[nm]
                sb = planck_int(l0, l1, TS)[0]
                tb = '  Tb_day %.0f  Tb_night %.0f' % tuple(
                    tb_invert(y[[i5, i0]]*sb/(RP/RS)**2, l0, l1))
            elif nm.startswith('ck'):
                b = int(nm[2:].split('_')[0])
                sb = planck_int(blo[b], bhi[b], TS)[0]
                tb = '  Tb_day %.0f  Tb_night %.0f' % tuple(
                    tb_invert(y[[i5, i0]]*sb/(RP/RS)**2, blo[b], bhi[b]))
            out.append('%-18s day %8.1f ppm  night %8.1f ppm  ampl %8.1f  peak offset %+6.2f deg'
                       ' (+ = east)  min at %+6.1f deg%s' % (row + (tb,)))
    open(O + 'phase.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
    np.savez(O + 'phase_curves.npz', phi=phi,
             **{('%s__%s' % (fn.split('_')[-1][:5], k)): v for fn, (t, f) in curves.items()
                for k, v in f.items()})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for i, (nm, obsd, obsn) in enumerate((('NRS1', 3924, 136), ('NRS2', 4924, 630),
                                          ('bol', None, None))):
        for fn, (t, f) in curves.items():
            ax[i].plot(phi, f[nm]*1e6, lw=1.2, label='rot %.1f' % (t/PROT))
        if obsd:
            ax[i].plot([0.5], [obsd], 'k*', ms=12, label='ME23 day')
            ax[i].plot([0.0, 1.0], [obsn, obsn], 'kx', ms=10, label='ME23 night')
        ax[i].set_title(nm + (' (Fp/F*, BB %d K star)' % TS))
        ax[i].set_xlabel('orbital phase')
        ax[i].set_ylabel('Fp/F* [ppm]')
        ax[i].axvline(0.5, color='0.7', lw=0.5)
        ax[i].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(O + 'phase_curves.png', dpi=120)
    print('wrote', O + 'phase_curves.png')


def rgrid():
    """cell-centre radii of the stretched grid (StretchRPoly of the uniform centres)"""
    xi = (np.arange(76) + 0.5)/76
    u = xi.copy()
    xk = xi.copy()
    for ck in CSTR:
        u += ck*xk*(1 - xi)
        xk *= xi
    return R0 + (R1 - R0)*u


def winds():
    from scipy.spatial import cKDTree
    raw = dhjcs.bin_convert.read_binary(BIN % 150)
    X, EXI, EET, dOm = geometry(raw)
    mb = raw['mb_data']
    rho = np.asarray(mb['dens'], float)
    eos = dhjcs.EOS(sys.argv[2] if len(sys.argv) > 2 else
                    '/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt')
    T, p = eos.invert(rho, np.asarray(mb['eint'], float))
    vr, v2, v3 = (np.asarray(mb[k], float) for k in ('velx', 'vely', 'velz'))
    rc = rgrid()
    # Cartesian velocity in the rotating frame + solid-body rotation (inertial, as observed)
    V = (vr[..., None]*X[:, :, :, None, :] + v2[..., None]*EXI[:, :, :, None, :]
         + v3[..., None]*EET[:, :, :, None, :])
    Vrot = OMEGA*rc[None, None, None, :, None]*np.stack(
        [-X[..., 1], X[..., 0], 0*X[..., 0]], -1)[:, :, :, None, :]
    ncol = X[..., 0].size
    Xf = X.reshape(ncol, 3)
    tree = cKDTree(Xf)
    rho_f = rho.reshape(ncol, -1)
    lp_f = np.log10(p.reshape(ncol, -1)/1e6)          # log10 bar
    T_f = T.reshape(ncol, -1)
    Vx_f = V[..., 0].reshape(ncol, -1)
    Vy_f = V[..., 1].reshape(ncol, -1)
    Vrx_f = Vrot[..., 0].reshape(ncol, -1)
    Vry_f = Vrot[..., 1].reshape(ncol, -1)
    windows = {'Fe (1e-5..1e-3 bar)': (-5, -3), 'Na (1e-6..1e-4 bar)': (-6, -4),
               'CO/H2O (1e-4..1e-2 bar)': (-4, -2), 'Halpha proxy (1e-8..1e-6, sponge)': (-8, -6)}
    out = ['# limb LOS velocity [km/s], observer on +x at mid-transit rotated by alpha;'
           ' RV = -(v_obs), negative = blueshift; rotation included',
           '# weights: ray = rays with tangent-point p in the window, rho ds along the ray; win = all'
           ' rays, rho ds of samples inside the window; rays uniform in position angle and b']
    th = np.deg2rad(np.arange(0, 360, 2.0) + 1.0)
    bb = np.linspace(rc[0], rc[-1], 120)
    s = np.linspace(-R1, R1, 801)
    ds = s[1] - s[0]
    res = {}
    for alpha in (0.0, -8.5, 8.5):
        a = np.deg2rad(alpha)
        # observer direction in the planet frame: at phase phi the sub-observer longitude is
        # 180 - 360 phi, i.e. the +x direction rotated by -alpha about z (alpha = 360 phi)
        ob = np.array([np.cos(a), -np.sin(a), 0.0])
        e1 = np.array([np.sin(a), np.cos(a), 0.0])     # sky axis in the orbital plane
        e2 = np.array([0.0, 0.0, 1.0])                 # sky axis, planet north
        acc = {}
        for it, t_ in enumerate(th):
            for b_ in bb:
                P0 = b_*(np.cos(t_)*e1 + np.sin(t_)*e2)
                pts = P0[None, :] + s[:, None]*ob[None, :]
                r = np.linalg.norm(pts, axis=1)
                ok = (r >= rc[0]) & (r <= rc[-1])
                if ok.sum() < 3:
                    continue
                pts, r = pts[ok], r[ok]
                _, col = tree.query(pts/r[:, None])
                fi = np.interp(r, rc, np.arange(76))
                i0 = np.clip(fi.astype(int), 0, 74)
                f = fi - i0

                def ip(A):
                    return A[col, i0]*(1 - f) + A[col, i0 + 1]*f
                rr = np.exp(np.log(rho_f[col, i0])*(1 - f) + np.log(rho_f[col, i0 + 1])*f)
                lpp = ip(lp_f)
                vobs = (ip(Vx_f) + ip(Vrx_f))*ob[0] + (ip(Vy_f) + ip(Vry_f))*ob[1]
                vwind = ip(Vx_f)*ob[0] + ip(Vy_f)*ob[1]
                lat = np.degrees(np.arcsin(np.sin(t_)*b_/b_))
                limb = 'morning' if np.cos(t_) > 0 else 'evening'   # +e1 = +y = lon -90 = west
                lptan = lpp.max()                    # tangent (closest-approach) level
                for wn, (l0, l1) in windows.items():
                    m = (lpp >= l0) & (lpp <= l1)
                    rsel = (lptan >= l0) & (lptan <= l1)
                    for sub in ('all', 'eq30'):
                        if sub == 'eq30' and abs(np.sin(t_)) > 0.5:
                            continue
                        # 'ray': rays whose tangent point lies in the window, rho ds along
                        # the whole ray; 'win': every ray, only samples inside the window
                        for wt, sel in (('ray', np.ones_like(m) if rsel else None),
                                        ('win', m if m.any() else None)):
                            if sel is None or b_ < rc[0]:
                                continue
                            k = (alpha, wn, limb, sub, wt)
                            q = acc.setdefault(k, np.zeros(4))
                            ww = rr[sel]*ds*b_
                            q += [ww.sum(), (ww*vobs[sel]).sum(), (ww*vwind[sel]).sum(),
                                  (ww*vobs[sel]**2).sum()]
                _ = lat
        res[alpha] = acc
    for alpha, acc in res.items():
        out.append('## alpha = %+.1f deg (planet rotation angle from mid-transit)' % alpha)
        for k in sorted(acc):
            W, s1, s2, s3 = acc[k]
            mean = s1/W
            sd = np.sqrt(max(s3/W - mean**2, 0))
            out.append('%-36s %-8s %-5s w=%-5s RV %+7.2f  (wind only %+7.2f, rotation %+6.2f)'
                       '  spread %5.2f' % (k[1], k[2], k[3], k[4], -mean/1e5, -s2/W/1e5,
                                            -(mean - s2/W)/1e5, sd/1e5))
        # both limbs together (full-transit proxy)
        for wn in windows:
            for wt in ('ray', 'win'):
                ks = [k for k in acc if k[1] == wn and k[3] == 'all' and k[4] == wt]
                W = sum(acc[k][0] for k in ks)
                if W > 0:
                    out.append('%-36s BOTH     all   w=%-5s RV %+7.2f' % (
                        wn, wt, -sum(acc[k][1] for k in ks)/W/1e5))
    # equatorial jet maxima at the limbs vs pressure (for the Na jet comparison)
    lat = np.degrees(np.arcsin(X[..., 2]))
    lon = np.degrees(np.arctan2(-X[..., 1], -X[..., 0]))
    u = v2*np.einsum('mkjx,mkjx->mkj', EXI, np.stack([-X[..., 1], X[..., 0], 0*X[..., 0]],
                     -1)/np.maximum(np.hypot(X[..., 0], X[..., 1]), 1e-30)[..., None])[..., None] \
        + v3*np.einsum('mkjx,mkjx->mkj', EET, np.stack([-X[..., 1], X[..., 0], 0*X[..., 0]],
                       -1)/np.maximum(np.hypot(X[..., 0], X[..., 1]), 1e-30)[..., None])[..., None]
    for lv in (-2, -3, -4, -5, -6, -7):
        ul = dhjcs.level_interp(u, np.log10(p/1e6), [lv])[0]
        for lnm, l0 in (('evening(+90)', 90.0), ('morning(-90)', -90.0)):
            m = (abs(lat) < 30) & (abs(lon - l0) < 15)
            out.append('# zonal u at 10^%d bar, |lat|<30, lon %s+-15: mean %+6.2f km/s, max %+6.2f'
                       % (lv, lnm, np.nanmean(ul[m])/1e5, np.nanmax(ul[m])/1e5))
    open(O + 'winds.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    {'phase': phase, 'winds': winds}[sys.argv[1]]()
