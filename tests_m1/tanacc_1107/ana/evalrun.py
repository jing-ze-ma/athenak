"""tanacc-1107 evaluator: python3 evalrun.py <run_dir> [...]  -> one json per run (<run>/eval.json)
Cells: radial shells nearest r = 2.0, 3.6, 4.5, 5.5 (and, for hot, every cell inside the clump), within
0.5 rad of the clump axis.  Errors are normalised by the unshadowed star field J0(r) = 2(1 - sqrt(1-1/r^2))
(E, J_fs, the tangential pressure anisotropy) and by 1/r^2 (F).  Regions by g = J_ex/J0: 'free' |g-1| < 0.02,
'aff' (affected by the clump) otherwise."""
import glob, json, math, os, sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/worktrees/tanacc-1107/vis/python')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc  # noqa: E402
import exact  # noqa: E402

CACHE = '/viper/ptmp2/jinma/openitems_1107/tanacc/ana/cache'


def load(run, var):
    f = sorted(glob.glob(f'{run}/bin/ta.{var}.*.bin') + glob.glob(f'{run}/ta.{var}.*.bin'))[-1]
    d = bc.read_binary(f)
    out = {k: np.concatenate([np.asarray(a).ravel() for a in d['mb_data'][k]])
           for k in d['var_names']}
    if 'mb_x1v' in d:
        x1 = []
        for m, a in enumerate(d['mb_data'][d['var_names'][0]]):
            x1v = np.asarray(d['mb_x1v'][m])
            x1.append(np.broadcast_to(x1v[None, None, :], np.asarray(a).shape).ravel())
        out['r'] = np.concatenate(x1)
    if 'r' not in out:   # uniform r: the volume centroid (3/4)(b^4 - a^4)/(b^3 - a^3)
        x1 = []
        for m, a in enumerate(d['mb_data'][d['var_names'][0]]):
            g = d['mb_geometry'][m]
            n1 = np.asarray(a).shape[-1]
            f = np.linspace(g[0], g[1], n1 + 1)
            rc = 0.75*(f[1:]**4 - f[:-1]**4)/(f[1:]**3 - f[:-1]**3)
            x1.append(np.broadcast_to(rc[None, None, :], np.asarray(a).shape).ravel())
        out['r'] = np.concatenate(x1)
    out['_cycle'] = d['cycle']
    out['_file'] = f
    return out


def cfg(run):
    p = {}
    for l in open(glob.glob(f'{run}/*.athinput')[0]):
        if '=' in l:
            k, v = l.split('#')[0].split('=', 1)
            p[k.strip()] = v.strip()
    return p


def main(run):
    p = cfg(run)
    m = load(run, 'm1')
    t = load(run, 'm1_lat')
    r = t.get('r', m.get('r'))
    rh = np.stack([t['m1_rx'], t['m1_ry'], t['m1_rz']], 1)
    X = r[:, None]*rh
    cs = 'use_cubed_sphere' in p
    clump = None
    if 'clump_a' in p:
        th, ph, rc = float(p['clump_th']), float(p['clump_ph']), float(p['clump_r'])
        xc = rc*np.array([math.sin(th)*math.cos(ph), math.sin(th)*math.sin(ph), math.cos(th)])
        chi = 1e-30*float(p['clump_rho'])
        B = float(p["arad"] if "arad" in p else 1.0)*float(p["clump_temp"])**4
        if B < 1e-60:
            B = 0.0
        clump = (xc, float(p['clump_a']), chi, B)
        ax = xc/rc
    else:
        ax = np.array([1.0, 0.0, 0.0]) if cs else np.array([math.cos(math.pi/4), math.sin(math.pi/4), 0])
    ru = np.unique(np.round(r, 10))
    shells = [ru[np.argmin(np.abs(ru - x))] for x in (2.0, 3.6, 4.5, 5.5)]
    sel = np.isin(np.round(r, 10), shells) & (rh @ ax > math.cos(0.5))
    if clump is not None and clump[3] > 0:
        sel |= np.linalg.norm(X - clump[0][None, :], axis=1) < clump[1]
    idx = np.where(sel)[0]
    key = '%s_%s_%d_%s' % ('cs' if cs else 'sp', p['nx1'], len(idx),
                           'none' if clump is None else '%g_%g' % (clump[2], clump[3]))
    os.makedirs(CACHE, exist_ok=True)
    cf = f'{CACHE}/{key}.npz'
    if os.path.exists(cf):
        z = np.load(cf)
        J, H, K = z['J'], z['H'], z['K']
    else:
        J, H, K = exact.moments(X[idx], clump)
        np.savez(cf, J=J, H=H, K=K)
    rs = r[idx]
    rhs = rh[idx]
    J0 = 2*(1 - np.sqrt(1 - 1/rs**2))
    g = J/J0
    E = m['m1_e'][idx]
    Fr = m['m1_f1'][idx]
    Hr = np.einsum('pi,pi->p', H, rhs)
    Kr = np.einsum('pij,pj->pi', K, rhs)
    Drr_ex = np.einsum('pi,pi->p', Kr, rhs)/np.maximum(J, 1e-300)
    # basis theta, phi
    thc = np.arccos(np.clip(rhs[:, 2], -1, 1))
    phc = np.arctan2(rhs[:, 1], rhs[:, 0])
    et = np.stack([np.cos(thc)*np.cos(phc), np.cos(thc)*np.sin(phc), -np.sin(thc)], 1)
    ep = np.stack([-np.sin(phc), np.cos(phc), 0*phc], 1)
    Ktt = np.einsum('pi,pij,pj->p', et, K, et)
    Kpp = np.einsum('pi,pij,pj->p', ep, K, ep)
    Ktp = np.einsum('pi,pij,pj->p', et, K, ep)
    q1e, q2e = 0.5*(Ktt - Kpp), Ktp                       # E units (x J)
    Qe = np.hypot(q1e, q2e)
    vlat_e = Kr - (np.einsum('pi,pi->p', Kr, rhs))[:, None]*rhs   # D_r,lat x J (Cartesian)
    tan_on = p.get('vet_gd_tangential', 'false') == 'true'
    L = {c: t['m1_lat%d' % c][idx] for c in range(6)}
    Drr = t['m1_tt0'][idx]   # slot 0 = D_rr after the VetLatBuild fold (checked vs exact)
    if cs:
        Qc = np.hypot(L[4], L[5])*E
        vc = np.stack([L[1], L[2], L[3]], 1)*E[:, None]
        dQ = np.abs(Qc - Qe)                                # magnitude only (frame-free)
    else:
        q1c, q2c = 0.5*(L[3] - L[5])*E, L[4]*E
        Qc = np.hypot(q1c, q2c)
        vc = (L[1]*E)[:, None]*et + (L[2]*E)[:, None]*ep
        dQ = np.hypot(q1c - q1e, q2c - q2e)
    if not tan_on:
        Qused = 0*Qc
        dQused = Qe
    else:
        Qused = Qc
        dQused = dQ
    Jfs = t['m1_jg'][idx]
    Hfs = t['m1_hg'][idx]
    res = {'run': run, 'cycle': int(m['_cycle']), 'n': int(len(idx)), 'tan': tan_on}
    regions = {'all': np.ones_like(g, bool), 'free': np.abs(g - 1) < 0.02, 'aff': np.abs(g - 1) >= 0.02}
    if clump is not None:
        regions['aff_far'] = (np.abs(g - 1) >= 0.02) & (rs > 4.0)
    for nm, R in regions.items():
        if R.sum() == 0:
            continue
        w = J0[R]
        d = {}
        d['N'] = int(R.sum())
        d['E_L1'] = float(np.mean(np.abs(E[R] - J[R])/w))
        d['E_bias'] = float(np.mean((E[R] - J[R])/w))
        d['E_max'] = float(np.max(np.abs(E[R] - J[R])/w))
        d['Jfs_L1'] = float(np.mean(np.abs(Jfs[R] - J[R])/w))
        d['E_vs_Jfs_L1'] = float(np.mean(np.abs(E[R] - Jfs[R])/w))
        fw = 1/rs[R]**2
        d['Fr_L1'] = float(np.mean(np.abs(Fr[R] - Hr[R])/fw))
        d['Hfs_L1'] = float(np.mean(np.abs(Hfs[R] - Hr[R])/fw))
        if not cs:
            Ft = m['m1_f2'][idx][R]
            Fp = m['m1_f3'][idx][R]
            Hte = np.einsum('pi,pi->p', H[R], et[R])
            Hpe = np.einsum('pi,pi->p', H[R], ep[R])
            d['Flat_L1'] = float(np.mean(np.hypot(Ft - Hte, Fp - Hpe)/fw))
            d['Flat_ex_mean'] = float(np.mean(np.hypot(Hte, Hpe)/fw))
        d['Drr_L1'] = float(np.mean(np.abs(Drr[R] - Drr_ex[R])*J[R]/w))
        d['Qex_mean'] = float(np.mean(Qe[R]/w))
        d['Qcode_mean'] = float(np.mean(Qc[R]/w))
        d['Q_fs_err'] = float(np.mean(dQ[R]/w))
        d['Q_used_err'] = float(np.mean(dQused[R]/w))
        d['vlat_err'] = float(np.mean(np.linalg.norm(vc[R] - vlat_e[R], axis=1)/w))
        d['vlat_ex'] = float(np.mean(np.linalg.norm(vlat_e[R], axis=1)/w))
        res[nm] = d
    json.dump(res, open(f'{run}/eval.json', 'w'), indent=1)
    a = res.get('aff', res['all'])
    print('%-28s cyc %d tan %d | all E_L1 %.3e | aff N %d E_L1 %.3e bias %+.3e Fr_L1 %.3e Jfs_L1 %.3e '
          'E-Jfs %.3e Qex %.3e Qused_err %.3e' % (os.path.basename(run), res['cycle'], tan_on,
                                                    res['all']['E_L1'], a['N'], a['E_L1'], a['E_bias'],
                                                    a['Fr_L1'], a['Jfs_L1'], a['E_vs_Jfs_L1'],
                                                    a['Qex_mean'], a['Q_used_err']))


if __name__ == '__main__':
    for rdir in sys.argv[1:]:
        try:
            main(rdir.rstrip('/'))
        except Exception as e:  # incomplete run (no final output) etc.
            print('%-28s SKIPPED: %r' % (os.path.basename(rdir.rstrip('/')), e))
