#!/usr/bin/env python3
"""zzzzd8 (READ ONLY): KE of CSHERO arms by component, radius bin and octahedral-symmetric part.  ke.py ARMDIR_Bp ARMDIR_L
Cubed sphere, 6 panels N x N equiangular.  hydro_w: v_r, v_2, v_3.  Two conventions for the lateral part are tested against hydro.hst:
 A: components on unit tangents (contravariant), KE_lat = v2^2 + v3^2 + 2 v2 v3 c, c = e1.e2 = -sin xi sin eta;
 B: covariant on the unit tangents: KE_lat = (v2^2 + v3^2 - 2 v2 v3 c) / (1 - c^2).
Symmetric part: v_sym(x) = (1/48) sum_R R^T v(R x) over the 48 signed permutation matrices (Cartesian vectors), rho_sym likewise."""
import os, re, sys, itertools, importlib.util, numpy as np
H = '/resnick/groups/carnegie_poc/jingze/bsglz_1124/'
sp = importlib.util.spec_from_file_location('imp', H + 'phase.py'); imp = importlib.util.module_from_spec(sp); sp.loader.exec_module(imp)
RS = 6.957e10; EDGES = [0, 100, 200, 400, 600, 800, 850, 900, 1000, 2000]
FR = {0: ((0, 1, 0), (0, 0, 1), (1, 0, 0)), 1: ((-1, 0, 0), (0, 0, 1), (0, 1, 0)), 2: ((0, -1, 0), (0, 0, 1), (-1, 0, 0)),
      3: ((0, 1, 0), (-1, 0, 0), (0, 0, 1)), 4: ((1, 0, 0), (0, 0, 1), (0, -1, 0)), 5: ((0, 1, 0), (1, 0, 0), (0, 0, -1))}
FA, FB, FN = (np.array([FR[p][q] for p in range(6)], float) for q in range(3))
def geom(N):
    e = (-1 + (np.arange(N) + 0.5) * 2 / N) * np.pi / 4; XI, ETA = np.meshgrid(e, e, indexing='xy')      # [k(eta), j(xi)]
    X, Y = np.tan(XI), np.tan(ETA); dl = np.sqrt(1 + X ** 2 + Y ** 2)
    Q = (FA[:, None, None, :] * X[None, :, :, None] + FB[:, None, None, :] * Y[None, :, :, None] + FN[:, None, None, :]) / dl[None, :, :, None]
    sx = 1 + X ** 2; sy = 1 + Y ** 2
    # unit tangents d Q / d xi, d Q / d eta (normalised)
    e1 = FA[:, None, None, :] * sx[None, :, :, None]; e2 = FB[:, None, None, :] * sy[None, :, :, None]
    def tang(d, dv):
        g = d / dl[None, :, :, None] - Q * (dv / dl ** 2)[None, :, :, None]
        return g / np.sqrt((g ** 2).sum(-1))[..., None]
    t1 = tang(e1, X * sx); t2 = tang(e2, Y * sy)
    return Q, t1, t2, (t1 * t2).sum(-1)[0]
def imgmap(N, Q):
    P = []
    flat = Q.reshape(-1, 3)
    for perm in itertools.permutations(range(3)):
        for s in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3)); R[range(3), perm] = s
            q = flat @ R.T                                         # R x for every cell centre x
            d = q @ FN.T; p = np.argmax(d, 1); qa = (q * FA[p]).sum(1) / (q * FN[p]).sum(1); qb = (q * FB[p]).sum(1) / (q * FN[p]).sum(1)
            j = np.rint((np.arctan(qa) / (np.pi / 4) + 1) * 0.5 * N - 0.5).astype(int); k = np.rint((np.arctan(qb) / (np.pi / 4) + 1) * 0.5 * N - 0.5).astype(int)
            idx = p * N * N + k * N + j
            assert len(np.unique(idx)) == len(idx), 'not a bijection'
            assert abs(flat[idx] - q).max() < 1e-10, 'cell centres not mapped'
            P.append((R, idx))
    assert len(P) == 48
    return P
def load(f):
    d = imp.read_bin(f); rho, rf, cs = imp.assemble(d, 'dens'); v1, v2, v3 = (imp.assemble(d, k)[0] for k in ('velx', 'vely', 'velz'))
    return d['time'], d['cycle'], rf, rho, v1, v2, v3
def rhst(f):
    n = None; rows = []
    for ln in open(f):
        if ln.startswith('#'):
            if '[1]=' in ln: n = [m.split('=', 1)[1] for m in re.findall(r'\[\d+\]=\S+', ln)]
        elif n and len(ln.split()) == len(n): rows.append([float(x) for x in ln.split()])
    return n, np.array(rows)
arms = {'Bp': sys.argv[1], 'L': sys.argv[2]}; geo = None; RES = {}
for a, d in arms.items():
    nh, hy = rhst(d + '/run_cs/rsg.hydro.hst'); nu, us = rhst(d + '/run_cs/rsg.user.hst')
    for c in (1, 2, 3, 4):
        f = d + '/run_cs/bin/rsg.hydro_w.%05d.bin' % c
        if not os.path.exists(f): continue
        t, cyc, rf, rho, v1, v2, v3 = load(f); N = rho.shape[1]; nr = rho.shape[3]
        if geo is None:
            geo = geom(N); IM = imgmap(N, geo[0]); w = imp.solid_angle(N); print('N %d, 48 images verified as bijections of the cell centres; max |e1.e2 + sin xi sin eta| checked via tangents' % N)
        Q, t1, t2, cg = geo; W = np.broadcast_to(w, (6, N, N)); dV = (rf[1:] ** 3 - rf[:-1] ** 3) / 3.0; r = 0.5 * (rf[1:] + rf[:-1]) / RS
        m = rho * W[..., None] * dV
        A = v2 ** 2 + v3 ** 2 + 2 * v2 * v3 * cg[None, :, :, None]; B = (v2 ** 2 + v3 ** 2 - 2 * v2 * v3 * cg[None, :, :, None]) / (1 - cg[None, :, :, None] ** 2)
        Kr = 0.5 * m * v1 ** 2
        k = int(np.argmin(abs(hy[:, 0] - t))); hk = [hy[k, nh.index(x)] for x in ('1-KE', '2-KE', '3-KE')]
        print('%s c%d t %.1f: dump KE_r %.6e, lat A %.6e, lat B %.6e | hydro.hst 1-KE %.6e, 2+3-KE %.6e, sum %.6e | dump sums A %.6e B %.6e' % (a, cyc, t, Kr.sum(), (0.5 * m * A).sum(), (0.5 * m * B).sum(), hk[0], hk[1] + hk[2], sum(hk), (Kr + 0.5 * m * A).sum(), (Kr + 0.5 * m * B).sum()))
        # Cartesian vectors with convention A (contravariant components on unit tangents); B would need raising first
        cov = abs((0.5 * m * B).sum() - (hk[1] + hk[2])) < abs((0.5 * m * A).sum() - (hk[1] + hk[2]))
        if cov:
            s2 = 1 - cg[None, :, :, None] ** 2; u2 = (v2 - cg[None, :, :, None] * v3) / s2; u3 = (v3 - cg[None, :, :, None] * v2) / s2
        else: u2, u3 = v2, v3
        V = v1[..., None] * Q[:, :, :, None, :] + u2[..., None] * t1[:, :, :, None, :] + u3[..., None] * t2[:, :, :, None, :]
        Vf = V.reshape(6 * N * N, nr, 3); rf_ = rho.reshape(6 * N * N, nr); Vs = np.zeros_like(Vf); rs_ = np.zeros_like(rf_)
        for R, idx in IM:
            Vs += Vf[idx] @ R; rs_ += rf_[idx]                    # R^T v(R x): row vector v R = (R^T v)^T
        Vs = (Vs / 48).reshape(V.shape); rs_ = (rs_ / 48).reshape(rho.shape)
        Vr = V - Vs; er = Q[:, :, :, None, :]
        ks_r = 0.5 * m * (Vs * er).sum(-1) ** 2; ks = 0.5 * m * (Vs ** 2).sum(-1); kn = 0.5 * m * (Vr ** 2).sum(-1); kx = m * (Vs * Vr).sum(-1); kt = 0.5 * m * (V ** 2).sum(-1); ktr = Kr
        v0 = (w[None, :, :, None] * v1).sum((0, 1, 2)) / (6 * w.sum()); Vc = Vs - v0[None, None, None, :, None] * er          # shell-mean radial flow (l = 0) removed from the symmetric field
        k0 = 0.5 * m * v0[None, None, None, :] ** 2; kc = 0.5 * m * (Vc ** 2).sum(-1); kc_r = 0.5 * m * (Vc * er).sum(-1) ** 2
        kn_r = 0.5 * m * (Vr * er).sum(-1) ** 2; sh = np.sqrt((v1 ** 2 + A if not cov else v1 ** 2 + B))
        RES[(a, cyc)] = []
        for lo, hi in list(zip(EDGES[:-1], EDGES[1:])) + [(0, 2000)]:
            q = (r >= lo) & (r < hi)
            if not q.any(): continue
            S = lambda x: float(x[..., q].sum())
            RES[(a, cyc)].append((lo, hi, S(kt), S(ktr), S(kt) - S(ktr), S(ks), S(ks_r), S(kn), S(kn_r), S(kx), float(sh[..., q].max()), S(np.broadcast_to(k0, m.shape)), S(kc), S(kc_r)))
print('\nlateral convention used: %s' % ('B (covariant)' if cov else 'A (contravariant, unit tangents)'))
for cyc in sorted({k[1] for k in RES}):
    print('\n# cycle %d: KE [erg, full sphere]; tot = r + lat; sym = 48-image mean field (rho not symmetrised in the weight), sym_r its radial part; rest = KE of v - v_sym; cross = sum m v_sym.(v - v_sym)' % cyc)
    print('%-12s %-3s %11s %11s %11s %11s %11s %11s %11s %11s %9s' % ('r [Rsun]', 'arm', 'total', 'radial', 'lateral', 'sym', 'sym radial', 'rest', 'rest radial', 'cross', '|v|max'))
    for i in range(len(RES[('Bp', cyc)])):
        for a in ('Bp', 'L'):
            if (a, cyc) in RES:
                x = RES[(a, cyc)][i]; print('%4d-%-6d %-3s %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.3e %9.3e | l0 %.4e cube %.4e cube_r %.4e' % (x[0], x[1], a, *x[2:]))
        if ('L', cyc) in RES:
            b, l = RES[('Bp', cyc)][i], RES[('L', cyc)][i]; print('%4d-%-6d %-3s %+10.1f%% %+10.1f%% %+10.1f%% %+10.1f%% %+10.1f%% %+10.1f%% %+10.1f%%' % (b[0], b[1], 'L-Bp', *[100 * (l[j] / b[j] - 1) for j in (2, 3, 4, 5, 6, 7, 8)]))
            print('%4d-%-6d d   %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e' % (b[0], b[1], *[b[j] - l[j] for j in (2, 3, 4, 5, 6, 7, 8)]))
# ---- hst
nb, ub = rhst(arms['Bp'] + '/run_cs/rsg.user.hst'); nl, ul = rhst(arms['L'] + '/run_cs/rsg.user.hst'); n = min(len(ub), len(ul))
print('\n# hst every 100 cycles (row = cycle, dt = 1000 s): KE_int Bp, L, rel diff | KEr_int rel diff | L_bot, L_mid rel diff | M_tot, Etot rel diff | Lz_tot Bp, L | L_midC / L_top rel diff')
for c in range(100, n, 100):
    g = lambda u, nn, k: u[c, nn.index(k)]
    rel = lambda k: (g(ul, nl, k) / g(ub, nb, k) - 1) if g(ub, nb, k) != 0 else float('nan')
    print('c%-5d t %.3e  %.4e %.4e %+.4f | %+.4f | %+.2e %+.2e | %+.2e %+.2e | %.3e %.3e | %+.2e %+.2e' % (c, ub[c, 0], g(ub, nb, 'KE_int'), g(ul, nl, 'KE_int'), rel('KE_int'), rel('KEr_int'), rel('L_bot'), rel('L_mid'), rel('M_tot'), rel('Etot'), g(ub, nb, 'Lz_tot'), g(ul, nl, 'Lz_tot'), rel('L_midC'), rel('L_top')))
print('KE_int gap vs t (Bp - L): log-log slope of the absolute gap, cycles 100-%d: %.2f' % (n - 1, np.polyfit(np.log(ub[100:n, 0]), np.log(np.maximum(ub[100:n, nb.index('KE_int')] - ul[100:n, nl.index('KE_int')], 1e-99)), 1)[0]))
