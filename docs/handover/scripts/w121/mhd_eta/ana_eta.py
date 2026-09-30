"""ana_eta.py: w121_mhd_0929 max_eta scan analysis (WASP-121b 1x, rot 300 -> resistive MHD).

Usage:  python3 ana_eta.py ARMDIR [ARMDIR ...]   (writes ARMDIR/ana_eta.txt, prints it)
        env ANA_BINS=all  -> every bin (default: first and last bin of the arm)

Per arm (inputs: ARMDIR/mhd.athinput, bin/dhj.mhd_w_bcc.*.bin, dhj.mhd.hst, run*.log):
 * cost:   wall per rotation (log 'elapsed=' lines), mean dt; dt limiter = Ohmic if the run dt is
           within 3 % of the Ohmic estimate cfl/6 * min(dx_d^2)/eta (resistivity.cpp, cfl from input)
 * eta:    x_e/T from the run's composition (coarse.bin = the EOS model on the run's own table grid,
           xetab.cpp of rmdx_0928), eta = ResistivityEOS capped (arm) and uncapped; per pressure band:
           f(eta at the cap), median eta, t_drag = eta/v_A^2 (capped / uncapped) vs t_adv = pi r/|v_h|
 * flow:   max |lat|<6 zonal-mean u and its p; equatorial (|lat|<20) u at 0.1/1 bar; day-night T at
           0.1 bar (dayside |lon|<90 minus nightside, area weighted) -- as ana_rot300/ana300.py
 * energy: resolved Ohmic heating Q = sum eta |curl B|^2 dV (B Heaviside-Lorentz; curl by a
           2nd-order Cartesian chain rule on the cell centres; grid-scale J under-estimated),
           Q/L_abs with L_abs = 4 pi R^2 sigma Teq^4; Lorentz work W = sum v.(J x B) dV;
           numerical eta_num = (W - Q - dME/dt)/sum J^2 dV (dME/dt from the hst around the
           bin; ignores the boundary Poynting flux: an ESTIMATE); div B max |divB| dr/|B|
"""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bin_convert as bc  # noqa: E402

SIGMA = 5.670374419e-5
PROT = 110153.5
T0 = 33046050.0
BANDS = [(1e-9, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 3e-2), (3e-2, 0.1), (0.1, 0.3),
         (0.3, 1.0), (1.0, 10.0), (10.0, 100.0), (100.0, 1e4)]
FRAMES = {0: ([0, 1, 0], [0, 0, 1], [1, 0, 0]), 1: ([-1, 0, 0], [0, 0, 1], [0, 1, 0]),
          2: ([0, -1, 0], [0, 0, 1], [-1, 0, 0]), 3: ([0, 1, 0], [-1, 0, 0], [0, 0, 1]),
          4: ([1, 0, 0], [0, 0, 1], [0, -1, 0]), 5: ([0, 1, 0], [1, 0, 0], [0, 0, -1])}


def panel_map(p, xi, eta):
    a, b, n = [np.array(v, float) for v in FRAMES[p]]
    x, y = np.tan(xi), np.tan(eta)
    dl = np.sqrt(1 + x * x + y * y)
    return (a * x[..., None] + b * y[..., None] + n) / dl[..., None]


class Tab:
    """log e, log p, mu, log x_e on a (log rho, log T) grid (xetab.cpp format)."""
    def __init__(self, fn):
        with open(fn, 'rb') as f:
            nx, ny = np.fromfile(f, np.int32, 2)
            self.x0, self.dx, self.y0, self.dy = np.fromfile(f, np.float64, 4)
            d = np.fromfile(f, np.float64).reshape(ny, nx, 4)
        self.nx, self.ny = nx, ny
        self.le, self.lp, self.lxe = d[..., 0].copy(), d[..., 1].copy(), d[..., 3].copy()

    def _x(self, rho):
        x = (np.log10(rho) - self.x0) / self.dx
        ix = np.clip(np.floor(x).astype(int), 0, self.nx - 2)
        return ix, np.clip(x - ix, 0.0, 1.0)

    def invert(self, rho, eint):
        ix, fx = self._x(rho)
        tgt = np.log10(eint / rho)

        def le_at(j):
            return self.le[j, ix] * (1 - fx) + self.le[j, ix + 1] * fx
        lo = np.zeros(rho.size, int)
        hi = np.full(rho.size, self.ny - 1)
        for _ in range(14):
            mid = (lo + hi) // 2
            up = le_at(mid) >= tgt
            hi = np.where(up, mid, hi)
            lo = np.where(up, lo, mid)
        e0, e1 = le_at(lo), le_at(hi)
        f = np.clip((tgt - e0) / np.where(e1 > e0, e1 - e0, 1.0), 0.0, 1.0)
        lp = (self.lp[lo, ix] * (1 - fx) + self.lp[lo, ix + 1] * fx) * (1 - f) + \
             (self.lp[hi, ix] * (1 - fx) + self.lp[hi, ix + 1] * fx) * f
        return 10 ** (self.y0 + (lo + f) * self.dy), 10 ** lp * rho

    def xe(self, rho, T):
        ix, fx = self._x(rho)
        y = (np.log10(T) - self.y0) / self.dy
        iy = np.clip(np.floor(y).astype(int), 0, self.ny - 2)
        fy = np.clip(y - iy, 0.0, 1.0)
        L = self.lxe
        return 10 ** ((L[iy, ix] * (1 - fx) + L[iy, ix + 1] * fx) * (1 - fy)
                      + (L[iy + 1, ix] * (1 - fx) + L[iy + 1, ix + 1] * fx) * fy)


def eta_eos(xe, T, cap):
    """Resistivity::ResistivityEOS (src/diffusion/resistivity.hpp)."""
    xu = np.maximum(np.maximum(xe, 230.0 * np.sqrt(T) / cap), 1e-300)
    return np.minimum(230.0 * np.sqrt(T) / xu + 5.2e11 * 20.0 / (T * np.sqrt(T)), cap)


def read_input(fn):
    blk, d = None, {}
    for line in open(fn):
        t = line.split('#')[0].strip()
        if t.startswith('<'):
            blk = t.strip('<>')
        elif '=' in t:
            k, v = [s.strip() for s in t.split('=', 1)]
            d[blk + '/' + k] = v
    return d


def redges(inp):
    n = int(inp['mesh/nx1'])
    a, b = float(inp['mesh/x1min']), float(inp['mesh/x1max'])
    c = [float(inp['mesh/f_stretch_r_c%d' % k]) for k in range(1, 9)
         if 'mesh/f_stretch_r_c%d' % k in inp]
    x = np.linspace(0.0, 1.0, n + 1)
    u, xk = x.copy(), x.copy()
    for ck in c:
        u += ck * xk * (1.0 - x)
        xk = xk * x
    return a + (b - a) * u


def geometry(fd, re):
    """per cell: unit rhat/exi/eeta (Cartesian), position, lat, lon, area weight, dx_d."""
    g = np.asarray(fd['mb_geometry'])
    nmb, nj, nk = fd['n_mbs'], fd['nx2_mb'], fd['nx3_mb']
    per = nmb // 6
    rc = 0.5 * (re[1:] + re[:-1])
    dr = np.diff(re)
    G = {k: np.zeros((nmb, nk, nj, 3)) for k in ('X', 'exi', 'eeta')}
    G.update({k: np.zeros((nmb, nk, nj)) for k in ('lat', 'lon', 'w', 'axi', 'aeta', 'cos')})
    for m in range(nmb):
        p = m // per
        xe_ = np.pi / 4 * (g[m, 2] + (g[m, 3] - g[m, 2]) * np.arange(nj + 1) / nj)
        ee_ = np.pi / 4 * (g[m, 4] + (g[m, 5] - g[m, 4]) * np.arange(nk + 1) / nk)
        XC, EC = np.meshgrid(0.5 * (xe_[1:] + xe_[:-1]), 0.5 * (ee_[1:] + ee_[:-1]))
        X = panel_map(p, XC, EC)
        h = 1e-6
        exi = panel_map(p, XC + h, EC) - panel_map(p, XC - h, EC)
        eet = panel_map(p, XC, EC + h) - panel_map(p, XC, EC - h)
        exi /= np.linalg.norm(exi, axis=-1)[..., None]
        eet /= np.linalg.norm(eet, axis=-1)[..., None]
        G['X'][m], G['exi'][m], G['eeta'][m] = X, exi, eet
        G['cos'][m] = (exi * eet).sum(-1)
        G['lat'][m] = np.arcsin(X[..., 2])
        G['lon'][m] = np.arctan2(-X[..., 1], -X[..., 0])     # 0 = substellar (-x), east +
        tx, ty = np.tan(XC), np.tan(EC)
        G['w'][m] = (1 + tx * tx) * (1 + ty * ty) / (1 + tx * tx + ty * ty) ** 1.5 \
            * (xe_[1] - xe_[0]) * (ee_[1] - ee_[0])
        XL, _ = np.meshgrid(xe_[:-1], EC[:, 0])
        XR, _ = np.meshgrid(xe_[1:], EC[:, 0])
        G['axi'][m] = np.arccos(np.clip((panel_map(p, XL, EC) * panel_map(p, XR, EC)).sum(-1),
                                        -1, 1))
        _, EL = np.meshgrid(XC[0], ee_[:-1])
        _, ER = np.meshgrid(XC[0], ee_[1:])
        G['aeta'][m] = np.arccos(np.clip((panel_map(p, XC, EL) * panel_map(p, XC, ER)).sum(-1),
                                         -1, 1))
    G['w'] *= 4 * np.pi / G['w'].sum()                          # solid angle per column
    G['rc'], G['dr'] = rc, dr
    G['dV'] = G['w'][..., None] * (rc ** 2 * dr)[None, None, None, :]
    return G


def curl_cart(Bc, P):
    """Bc, P: (nmb,nk,nj,ni,3) -> curl B (Cartesian) by d/dx = J^-T d/d(index)."""
    dP = np.stack([np.gradient(P, axis=a) for a in (1, 2, 3)], -1)   # (...,3 xyz,3 idx)
    Minv = np.linalg.inv(dP)                                           # (...,3 idx, 3 xyz)
    grad = np.zeros(Bc.shape + (3,))                                   # dB_a/dx_b
    for a in range(3):
        dB = np.stack([np.gradient(Bc[..., a], axis=ax) for ax in (1, 2, 3)], -1)
        grad[..., a, :] = np.einsum('...i,...ib->...b', dB, Minv)
    return np.stack([grad[..., 2, 1] - grad[..., 1, 2], grad[..., 0, 2] - grad[..., 2, 0],
                     grad[..., 1, 0] - grad[..., 0, 1]], -1)


def lev(q, lp, L):
    """q, log p (..., ni) with p decreasing in i -> q at log p = L (nan outside)."""
    idx = (lp > L).sum(-1)
    ok = (idx > 0) & (idx < q.shape[-1])
    i1 = np.clip(idx, 1, q.shape[-1] - 1)
    a = np.take_along_axis(lp, (i1 - 1)[..., None], -1)[..., 0]
    b = np.take_along_axis(lp, i1[..., None], -1)[..., 0]
    qa = np.take_along_axis(q, (i1 - 1)[..., None], -1)[..., 0]
    qb = np.take_along_axis(q, i1[..., None], -1)[..., 0]
    f = np.clip((a - L) / np.where(a > b, a - b, 1.0), 0, 1)
    return np.where(ok, qa + f * (qb - qa), np.nan)


def hst(fn):
    if not os.path.exists(fn):
        return None
    with open(fn) as f:
        f.readline()
        names = [c.split('=')[1] for c in f.readline().split()[1:]]
    d = np.loadtxt(fn, ndmin=2)
    return {n: d[:, i] for i, n in enumerate(names) if i < d.shape[1]}


def cost(arm):
    rows = []
    for fn in sorted(glob.glob(arm + '/run*.log')):
        seg = []
        for line in open(fn):
            m = re.match(r'elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)', line)
            if m:
                seg.append([float(m.group(1)), int(m.group(2)), float(m.group(3)),
                            float(m.group(4))])
        if len(seg) > 2:
            rows.append(np.array(seg))
    if not rows:
        return None
    wall = sum(s[-1, 0] - s[0, 0] for s in rows)
    trot = sum(s[-1, 2] - s[0, 2] for s in rows) / PROT
    cyc = sum(s[-1, 1] - s[0, 1] for s in rows)
    allr = np.vstack(rows)
    return dict(wall_per_rot=wall / max(trot, 1e-30), rot=trot, cycles=cyc,
                ms_per_cycle=1e3 * wall / max(cyc, 1), dt_mean=np.mean(allr[:, 3]),
                dt_min=np.min(allr[:, 3]), rot_end=(allr[-1, 2] - T0) / PROT)


def analyse_bin(fn, inp, tab, G, cap, H):
    fd = bc.read_binary(fn)
    md = fd['mb_data']
    f64 = lambda k: np.asarray(md[k], dtype=np.float64)  # noqa: E731
    rho = f64('dens')
    sh = rho.shape
    T, p = tab.invert(rho.ravel(), f64('eint').ravel())
    T, p = T.reshape(sh), p.reshape(sh)
    xe = tab.xe(rho.ravel(), T.ravel()).reshape(sh)
    eta = eta_eos(xe, T, cap)
    eta_u = eta_eos(xe, T, 1e30)
    atcap = eta_u >= cap * (1 - 1e-9)
    v1, v2, v3 = f64('velx'), f64('vely'), f64('velz')
    b1, b2, b3 = f64('bcc1'), f64('bcc2'), f64('bcc3')
    X, exi, eet = [G[k][:, :, :, None, :] for k in ('X', 'exi', 'eeta')]
    c = G['cos'][..., None, None]
    s = np.sqrt(1 - c * c)
    V = v1[..., None] * X + v2[..., None] * exi + v3[..., None] * eet
    Bc = b1[..., None] * X + b2[..., None] * exi + b3[..., None] * (eet - c * exi) / s
    P = X * G['rc'][None, None, None, :, None]
    J = curl_cart(Bc, P)
    J2 = (J * J).sum(-1)
    B2 = (Bc * Bc).sum(-1)
    dV = G['dV']
    Q = np.sum(eta * J2 * dV)
    W = np.sum((V * np.cross(J, Bc)).sum(-1) * dV)
    J2V = np.sum(J2 * dV)
    vh = np.sqrt(np.maximum(((V - (V * X).sum(-1)[..., None] * X) ** 2).sum(-1), 1e-30))
    rc = G['rc'][None, None, None, :]
    t_adv = np.pi * rc / vh
    vA2 = np.maximum(B2 / rho, 1e-300)
    tdr, tdr_u = eta / vA2, eta_u / vA2
    # Ohmic dt estimate, as Resistivity::NewTimeStepGeneralResist x cfl
    cfl = float(inp.get('time/cfl_number', 0.3))
    dxi = rc * G['axi'][..., None]
    det = rc * G["aeta"][..., None]
    dmin = np.minimum(np.broadcast_to(G['dr'], sh), np.minimum(dxi, det))
    tdo = cfl / 6.0 * dmin ** 2 / eta
    io = np.unravel_index(np.argmin(tdo), sh)
    pbar = p / 1e6
    out = {'time': fd['time'], 'rot': (fd['time'] - T0) / PROT + 300.0, 'Q': Q, 'W': W,
           'J2V': J2V, 'EM': 0.5 * np.sum(B2 * dV), 'EK': 0.5 * np.sum(rho * (V * V).sum(-1) * dV),
           'dt_ohm': tdo[io], 'dt_ohm_p': pbar[io], 'dt_ohm_lat': np.degrees(G['lat'][io[:3]]),
           'dt_ohm_lon': np.degrees(G['lon'][io[:3]]), 'Bmax': np.sqrt(B2.max())}
    rows = []
    for lo, hi in BANDS:
        k = (pbar >= lo) & (pbar < hi)
        if not k.any():
            continue
        wk = dV * k
        kc = k & atcap                       # the cells where the cap acts
        lr = np.log10(tdr / t_adv)
        lru = np.log10(tdr_u / t_adv)
        med = lambda q, msk: np.median(q[msk]) if msk.any() else np.nan  # noqa: E731
        rows.append((lo, hi, k.sum(), np.sum(atcap * wk) / wk.sum(),
                     np.median(np.log10(eta[k])), np.median(np.log10(eta_u[k])),
                     med(lr, kc), med(lru, kc),
                     np.sum((tdr < t_adv) * kc * dV) / max(np.sum(kc * dV), 1e-300),
                     med(lr, k), np.median(np.sqrt(B2[k])), np.sum(eta * J2 * wk)))
    out['bands'] = rows
    # cap sweep on THIS state (B, v, T fixed): Ohmic dt and the volume fraction per band
    # where the cap creates drag the physical eta would not (t_drag(cap) < t_adv <= t_drag(eta))
    sweep = []
    for cs_ in (1e11, 1e12, 1e13, 1e14, 1e15):
        es = np.minimum(eta_u, cs_)
        sp = (es / vA2 < t_adv) & (tdr_u >= t_adv)
        fr = []
        for lo, hi in BANDS:
            k = (pbar >= lo) & (pbar < hi)
            fr.append(np.sum(sp * k * dV) / max(np.sum(k * dV), 1e-300))
        sweep.append((cs_, np.min(cfl / 6.0 * dmin ** 2 / es), fr))
    out['sweep'] = sweep
    # circulation, as ana300.py
    lat, lon = G['lat'], G['lon']
    zon = np.stack([-G['X'][..., 1], G['X'][..., 0], np.zeros_like(lat)], -1)
    zon /= np.maximum(np.linalg.norm(zon, axis=-1), 1e-30)[..., None]
    u = (V * zon[:, :, :, None, :]).sum(-1)
    lp = np.log10(pbar)
    w = G['w']
    eq6, eq20 = np.abs(lat) < np.radians(6), np.abs(lat) < np.radians(20)
    PF = np.logspace(-5, 2.5, 61)
    ujet = []
    for L in np.log10(PF):
        ul = lev(u, lp, L)
        ok = np.isfinite(ul) & eq6
        ujet.append(np.sum(np.where(ok, ul, 0) * w) / max(np.sum(w * ok), 1e-30))
    ujet = np.array(ujet) / 100.0
    out['ujet'], out['pjet'] = np.nanmax(ujet), PF[np.nanargmax(ujet)]
    for P_ in (1e-3, 0.1, 1.0):
        ul = lev(u, lp, np.log10(P_))
        ok = np.isfinite(ul) & eq20
        out['ueq_%g' % P_] = np.sum(np.where(ok, ul, 0) * w) / np.sum(w * ok) / 100.0
    T1 = lev(T, lp, -1.0)
    day = np.abs(lon) < np.pi / 2
    ok = np.isfinite(T1)
    out['Tday'] = np.sum(np.where(ok & day, T1, 0) * w) / np.sum(w * (ok & day))
    out['Tnight'] = np.sum(np.where(ok & ~day, T1, 0) * w) / np.sum(w * (ok & ~day))
    # dME/dt from the hst around this bin (+-0.25 rot; hst cadence is 0.1 rot -- DeltaAI local fix, was 0.1)
    out['dEMdt'] = np.nan
    if H is not None and '1-ME' in H:
        me = H['1-ME'] + H['2-ME'] + H['3-ME']
        sel = np.abs(H['time'] - fd['time']) < 0.25 * PROT
        if sel.sum() >= 3:
            out['dEMdt'] = np.polyfit(H['time'][sel], me[sel], 1)[0]
    return out


def main():
    tab = Tab(os.path.join(HERE, 'coarse.bin'))
    for arm in sys.argv[1:]:
        arm = arm.rstrip('/')
        inp = read_input(arm + '/mhd.athinput')
        cap = float(inp['mhd/max_eta'])
        teq = float(inp['problem/Teq'])
        ap = float(inp['problem/ap'])
        Labs = 4 * np.pi * ap ** 2 * SIGMA * teq ** 4
        re = redges(inp)
        H = hst(arm + '/dhj.mhd.hst')
        bins = sorted(glob.glob(arm + '/bin/dhj.mhd_w_bcc.*.bin'))
        if os.environ.get('ANA_BINS', '') != 'all' and len(bins) > 2:
            bins = [bins[0], bins[-1]]
        L = ['# arm %s: bbot %s (code), max_eta %.1e, use_rkg_sts %s, rsolver %s' % (
            arm, inp['problem/bbot'], cap, inp.get('mhd/use_rkg_sts', 'false'),
            inp['mhd/rsolver'])]
        C = cost(arm)
        if C:
            L.append('# cost: %.2f rot run (to rot %.2f), wall %.1f min/rot, %.2f ms/cycle, '
                     'dt mean %.3f s min %.3f s' % (C['rot'], 300 + C['rot_end'],
                                                    C['wall_per_rot'] / 60, C['ms_per_cycle'],
                                                    C['dt_mean'], C['dt_min']))
        G = None
        for fn in bins:
            if G is None:
                G = geometry(bc.read_binary(fn), re)
            o = analyse_bin(fn, inp, tab, G, cap, H)
            dtrun = np.nan
            if H is not None:
                dtrun = np.interp(o['time'], H['time'], H['dt'])
            lim = 'Ohmic' if abs(dtrun / o['dt_ohm'] - 1) < 0.03 else 'CFL/other'
            L.append('## %s rot %.3f: dt run %.3f s, Ohmic estimate %.3f s at p %.3g bar lat %+.0f '
                     'lon %+.0f -> limiter %s' % (os.path.basename(fn), o['rot'], dtrun,
                                                  o['dt_ohm'], o['dt_ohm_p'], o['dt_ohm_lat'],
                                                  o['dt_ohm_lon'], lim))
            L.append('   jet max %.0f m/s at %.2g bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = %.0f/%.0f/%.0f '
                     'm/s; T(0.1 bar) day %.0f night %.0f contrast %.0f K'
                     % (o['ujet'], o['pjet'], o['ueq_0.001'], o['ueq_0.1'], o['ueq_1'],
                        o['Tday'], o['Tnight'], o['Tday'] - o['Tnight']))
            eta_num = (o['W'] - o['Q'] - o['dEMdt']) / o['J2V']
            L.append('   E_mag %.3e E_kin %.3e erg, |B|max %.3g; Q_ohm %.3e erg/s = %.2e L_abs; '
                     'W_lorentz %.3e; dEM/dt %.3e; eta_num(est) %.2e cm2/s'
                     % (o['EM'], o['EK'], o['Bmax'], o['Q'], o['Q'] / Labs, o['W'], o['dEMdt'],
                        eta_num))
            L.append('   p band [bar]  ncell f(cap) medlog eta cap/uncap | AT-CAP cells: medlog '
                     't_drag/t_adv cap/uncap, f(t_drag<t_adv) | all: medlog t_drag/t_adv, '
                     'med|B|, Q_band/L')
            for r in o['bands']:
                L.append('   %7.0e-%7.0e %6d %5.3f %6.2f %6.2f | %7.2f %7.2f %6.3f | %7.2f '
                         '%8.3g %8.2e' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7],
                                          r[8], r[9], r[10], r[11] / Labs))
        if G is not None:
            L.append('# cap sweep on the last bin (state fixed): cap, Ohmic dt [s], then per '
                     'band the volume fraction where the cap creates spurious drag')
            L.append('#   bands: ' + ' '.join('%.0e-%.0e' % b for b in BANDS))
            for cs_, dto, fr in o['sweep']:
                L.append('   %.0e %8.3f  ' % (cs_, dto) + ' '.join('%.3f' % x for x in fr))
        # div B
        db = sorted(glob.glob(arm + '/bin/dhj.mhd_divb.*.bin'))
        if db and G is not None:
            fd = bc.read_binary(db[-1])
            dv = np.abs(np.asarray(fd['mb_data']['divb'], float))
            fw = bc.read_binary(bins[-1])['mb_data']
            B = np.sqrt(sum(np.asarray(fw[k], float) ** 2 for k in ('bcc1', 'bcc2', 'bcc3')))
            if dv.shape == B.shape:
                L.append('# div B (%s): max |divB| dr/|B| = %.2e' % (
                    os.path.basename(db[-1]),
                    np.max(dv * G['dr'][None, None, None, :] / np.maximum(B, 1e-30))))
        txt = '\n'.join(L)
        open(arm + '/ana_eta.txt', 'w').write(txt + '\n')
        print(txt)


if __name__ == '__main__':
    main()
