#!/usr/bin/env python3
"""zzzzd5 item 2 (dump part) and item 4 (READ ONLY).  groups.py DUMP...  bsg sp wedge, hydro_w float32 primitives (post-merge state).
merge_lat_rc 2.43495e12 (f = 2), 3.82635e12 (f = 1), 5.5656e12 (f = 2): groups 2x2 lateral, block aligned."""
import sys, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
RC = (2.43495e12, 3.82635e12, 5.5656e12); F = (2, 1, 2)
def fac(r):
    return 2 if r < RC[0] else (1 if r < RC[1] else (2 if r < RC[2] else 1))
for f in sys.argv[1:]:
    h = bc.read_binary(f); out = []
    Lz = {k: 0.0 for k in ('all', 'm1', 'f1', 'm2')}; Lg = dict(Lz); M = dict(Lz); tot_g = 0
    ng = dict(inner=0, outer=0); c = {k: np.zeros(6, int) for k in ('inner', 'outer')}; sp = {'inner': [], 'outer': []}; spread = []
    for b in range(h['n_mbs']):
        rr = np.asarray(h['mb_x1v'][b], float); rf = np.asarray(h['mb_x1f'][b], float); th = np.asarray(h['mb_x2v'][b], float); thf = np.asarray(h['mb_x2f'][b], float); ph = np.asarray(h['mb_x3v'][b], float); phf = np.asarray(h['mb_x3f'][b], float)
        D = {k: np.asarray(h['mb_data'][k][b], float) for k in ('dens', 'velx', 'vely', 'velz', 'eint')}
        dV = np.diff(phf)[:, None, None] * (np.cos(thf[:-1]) - np.cos(thf[1:]))[None, :, None] * ((rf[1:] ** 3 - rf[:-1] ** 3) / 3.)[None, None, :]
        rho = D['dens']; sint = np.sin(th)[None, :, None]; lzc = rho * rr[None, None, :] * sint * D['velz'] * dV
        fr = np.array([fac(x) for x in rr])
        for nm, q in (('all', np.ones(len(rr), bool)), ('m1', rr < RC[0]), ('f1', (rr >= RC[0]) & (rr < RC[1])), ('m2', (rr >= RC[1]) & (rr < RC[2]))):
            Lz[nm] += lzc[:, :, q].sum(); Lg[nm] += abs(lzc[:, :, q]).sum(); M[nm] += (rho * dV)[:, :, q].sum()
        # Cartesian velocity (physical r, theta, phi components) and positions
        st, ct, sph, cph = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
        er = np.stack([st[None, :] * cph[:, None], st[None, :] * sph[:, None], ct[None, :] * np.ones_like(ph)[:, None]], -1)      # (k, j, 3)
        et = np.stack([ct[None, :] * cph[:, None], ct[None, :] * sph[:, None], -st[None, :] * np.ones_like(ph)[:, None]], -1)
        ep = np.stack([-sph[:, None] * np.ones_like(th)[None, :], cph[:, None] * np.ones_like(th)[None, :], np.zeros((len(ph), len(th)))], -1)
        nk, nj = len(ph), len(th)
        for i in np.where(fr == 2)[0]:
            w = 'inner' if rr[i] < RC[0] else 'outer'
            v = D['velx'][:, :, i, None] * er + D['vely'][:, :, i, None] * et + D['velz'][:, :, i, None] * ep
            x = rr[i] * er; dv = dV[:, :, i]; rh = rho[:, :, i]; ei = D['eint'][:, :, i]
            g = lambda q: q.reshape(nk // 2, 2, nj // 2, 2, *q.shape[2:]).swapaxes(1, 2).reshape(nk // 2, nj // 2, 4, *q.shape[2:])
            V, X, dw, R_, E_ = g(v), g(x), g(dv), g(rh), g(ei)
            vs = dw.sum(2); xc = (dw[..., None] * X).sum(2) / vs[..., None]
            vb = (dw[..., None] * V).sum(2) / vs[..., None]; d = X - xc[:, :, None, :]
            Lm = (dw[..., None] * np.cross(d, V)).sum(2)                                # per unit rho
            d2 = (d ** 2).sum(-1); Jt = (dw[..., None, None] * (d2[..., None, None] * np.eye(3) - d[..., :, None] * d[..., None, :])).sum(2)
            tr = np.trace(Jt, axis1=-2, axis2=-1); det = np.linalg.det(Jt)
            om = np.zeros_like(Lm); ok = det > 1e-12 * tr ** 3
            om[ok] = np.linalg.solve(Jt[ok], Lm[ok][..., None])[..., 0]
            rbar = (dw * R_).sum(2) / vs; Er = 0.5 * rbar * (om * Lm).sum(-1); Eth = (dw * E_).sum(2)
            rig = np.cross(om[:, :, None, :], d); res = V - vb[:, :, None, :] - rig
            Eres = 0.5 * rbar * (dw[..., None] * res ** 2).sum((2, 3))
            dvm = np.sqrt(((V - vb[:, :, None, :]) ** 2).sum(-1)).max(2); vm = np.sqrt((vb ** 2).sum(-1))
            uni = (dvm <= 3e-7 * np.maximum(vm, 1e-30))                                 # cell velocities equal to float32 round-off = zero rotation
            c[w][0] += ok.size; c[w][1] += int((Er > 0.5 * Eth).sum()); c[w][2] += int((Er > Eth).sum()); c[w][3] += int((Er > 0.5 * (Eth + Er + Eres)).sum()); c[w][4] += int(uni.sum()); c[w][5] += int((~ok).sum())
            spread.append(float((abs(R_ - R_[:, :, :1]).max(2) / R_[:, :, 0]).max()))
            m = Er > 0.5 * Eth
            if m.any():
                for q in np.argwhere(m)[:200000]:
                    sp[w].append((rr[i], rbar[tuple(q)], Eth[tuple(q)] / vs[tuple(q)], Er[tuple(q)] / Eth[tuple(q)]))
    print('FILE %s t %.5e cycle %d' % (f.split('jingze/')[1], h['time'], h['cycle']))
    print('  wedge M all %.4e; r < 2.435e12 %.4e; 2.435-3.826e12 %.4e; 3.826-5.566e12 %.4e (wedge, theta 60-120, phi 0-60 deg, not scaled)' % tuple(M[k] for k in ('all', 'm1', 'f1', 'm2')))
    for nm in ('all', 'm1', 'f1', 'm2'):
        print('  L_z net (sum rho r sin(theta) v_phi dV) %s %.4e, gross (sum of |.|) %.4e' % (nm, Lz[nm], Lg[nm]))
    print('  max relative spread of rho inside a merged group over the dump: %.2e (post-merge state if ~0)' % max(spread))
    for w in ('inner', 'outer'):
        print('  groups %s: total %d; Erot > 0.5 Eth: %d; Erot > Eth: %d; code test Erot > 0.5 (Eth + Erot + Eres): %d; all cells same velocity to float32 (zero rotation): %d; singular J: %d' % (w, c[w][0], c[w][1], c[w][2], c[w][3], c[w][4], c[w][5]))
        if sp[w]:
            q = np.array(sp[w]); print('    Erot > 0.5 Eth set: r %.4e-%.4e cm, rho %.3e-%.3e, e_int density %.3e-%.3e, Erot/Eth max %.3f' % (q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max(), q[:, 2].min(), q[:, 2].max(), q[:, 3].max()))
