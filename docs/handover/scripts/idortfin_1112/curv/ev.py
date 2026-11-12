"""idortcurv-1112 evaluator: python3 ev.py <run_dir> [...] -> one line per run + <run>/ev.json
T1-T3 (thin star / clump): errors vs the exact ray quadrature (tanacc exact.py) on every cell (star:
analytic) or on the 4 shells near r = 2, 3.6, 4.5, 5.5 within 0.5 rad of the clump axis (clump);
cs cells split by panel position (interior / seam = within 2 cells of a panel edge / corner = both).
T4 (core): homogeneous sphere exact moments (1-D quadrature), every cell.
Counters from run.log."""
import glob, json, math, os, re, sys
import numpy as np
sys.path.insert(0, os.environ.get('ATHENAK_VIS', '/orion/ptmp/jinma/Athenak/wt_fin/vis/python'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc  # noqa: E402
import exact  # noqa: E402

CACHE = os.environ.get('EV_CACHE', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))


def load(run, var):
    fs = sorted(glob.glob(f'{run}/bin/ta.{var}.*.bin'))
    d = bc.read_binary(fs[-1])
    out = {k: np.concatenate([np.asarray(a).ravel() for a in d['mb_data'][k]]) for k in d['var_names']}
    x1 = []
    for m, a in enumerate(d['mb_data'][d['var_names'][0]]):
        g = d['mb_geometry'][m]
        n1 = np.asarray(a).shape[-1]
        f = np.linspace(g[0], g[1], n1 + 1)
        rc = 0.75*(f[1:]**4 - f[:-1]**4)/(f[1:]**3 - f[:-1]**3)
        x1.append(np.broadcast_to(rc[None, None, :], np.asarray(a).shape).ravel())
    out['r'] = np.concatenate(x1)
    out['_cycle'] = d['cycle']
    out['_time'] = d['time']
    return out


def cfg(run):
    p = {}
    for ln in open(f'{run}/in.athinput'):
        if '=' in ln:
            k, v = ln.split('#')[0].split('=', 1)
            p[k.strip()] = v.strip()
    return p


PATS = {'notadm': r'NOT ADMISSIBLE', 'notpos': r'NOT POSITIVE', 'fatal': r'FATAL'}


def counters(run):
    c = {}
    txt = open(f'{run}/run.log', errors='replace').read()
    for k, p in PATS.items():
        c[k] = len(re.findall(p, txt))
    m = re.findall(r'stage steps=([0-9.e+]+).*?backward-Euler steps=([0-9.e+]+)', txt)
    if m:
        c['stage'], c['be'] = m[-1]
    for k, p in (('fclip', r'FCLIPDIAG final: n\(metric\)=([0-9.e+-]+)'),
                 ('fclipm', r'FCLIPDIAG final: n\(metric\)=[0-9.e+-]+ mean\(\|F\|/cE-1\)=([0-9.e+-]+)'),
                 ('realis', r'rescale: lag \(beta.F of E0\) cell-passes=([0-9.e+-]+)'),
                 ('drag', r'write-back drag cell-solves=([0-9.e+-]+)'),
                 ('ncv', r'NON-CONVERGED=([0-9.e+-]+)'),
                 ('picard', r'Picard iterations mean=([0-9.e+-]+)'),
                 ('floorE', r'solved E <= e_floor cell-solves=([0-9.e+-]+)'),
                 ('flooreint', r'gas eint <= 0 cell-solves=([0-9.e+-]+)'),
                 ('stfb', r'stage fallbacks=([0-9.e+-]+)'),
                 ('posfb', r'positivity fallbacks=([0-9.e+-]+)')):
        m = re.findall(p, txt)
        if m:
            c[k] = '%g' % float(m[-1])
    m = re.findall(r'cycle=(\d+)', txt)
    c['cyc'] = int(m[-1]) if m else -1
    m = re.findall(r'zone-cycles/cpu_second = ([0-9.e+-]+)', txt)
    if m:
        c['zcps'] = float(m[-1])
    m = re.findall(r'cpu time used\s*=\s*([0-9.e+-]+)', txt)
    if m:
        c['wall'] = float(m[-1])
    return c


def panel_region(rh, n2):
    a = np.sort(np.abs(rh), axis=1)[:, ::-1]
    xi = np.arctan(a[:, 1]/a[:, 0])
    eta = np.arctan(a[:, 2]/a[:, 0])
    d = 2*(math.pi/2)/n2
    s1 = xi > math.pi/4 - d
    s2 = eta > math.pi/4 - d
    return np.where(s1 & s2, 2, np.where(s1 | s2, 1, 0))


def sphere_moments(r, R, chi, B, n=400):
    """homogeneous sphere (radius R, extinction chi, source B): J, H_r, K_rr at radii r."""
    J, H, K = np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)
    x, w = np.polynomial.legendre.leggauss(n)
    for i, ri in enumerate(r):
        if ri < R:
            segs = [(-1.0, 1.0)]
        else:
            segs = [(math.sqrt(1 - (R/ri)**2), 1.0)]
        for a, b in segs:
            mu = a + (b - a)*0.5*(x + 1)
            ww = (b - a)*0.5*w
            disc = np.clip(ri**2*mu**2 - ri**2 + R**2, 0, None)
            L = ri*mu + np.sqrt(disc) if ri < R else 2*np.sqrt(disc)
            I = B*(1 - np.exp(-chi*L))
            J[i] += 0.5*np.sum(ww*I)
            H[i] += 0.5*np.sum(ww*I*mu)
            K[i] += 0.5*np.sum(ww*I*mu**2)
    return J, H, K


def stats(err, w, reg, names):
    o = {}
    for k, nm in enumerate(names):
        s = reg == k
        if s.sum():
            o[nm] = (int(s.sum()), float(np.mean(np.abs(err[s])/w[s])), float(np.mean(err[s]/w[s])),
                     float(np.max(np.abs(err[s])/w[s])))
    return o


def thin(run, p, res):
    m = load(run, 'm1')
    t = load(run, 'm1_lat')
    r = m['r']
    rh = np.stack([t['m1_rx'], t['m1_ry'], t['m1_rz']], 1)
    cs = p.get('use_cubed_sphere') == 'true'
    n2 = int(p['nx2'])
    E, Fr = m['m1_e'], m['m1_f1']
    J0 = 2*(1 - np.sqrt(1 - 1/r**2))
    reg = panel_region(rh, n2) if cs else np.zeros(len(r), int)
    names = ['int', 'seam', 'corner'] if cs else ['all']
    res['fr_over_cE_max'] = float(np.max(Fr/np.maximum(E, 1e-300)))
    res['n_fr_gt_cE'] = int(np.sum(Fr > E*(1 + 1e-12)))
    if 'clump_a' not in p:
        H = 1/r**2
        mu0 = np.sqrt(1 - 1/r**2)
        Drr = (1 + mu0 + mu0**2)/3
        res['E'] = stats(E - J0, J0, reg, names)
        res['Fr'] = stats(Fr - H, H, reg, names)
        res['Drr'] = stats(t['m1_tt0'] - Drr, Drr, reg, names)
        far = r > 5.0
        res['f_far'] = float(np.mean(Fr[far]/E[far]))
        res['f_far_ex'] = float(np.mean(H[far]/J0[far]))
        res['L_far'] = float(np.mean(Fr[far]*r[far]**2))
        line = 'E_L1 ' + ' '.join('%s %.4f/%+.4f/%.3f' % (k, v[1], v[2], v[3]) for k, v in res['E'].items())
        line += ' | Fr_L1 ' + ' '.join('%s %.4f' % (k, v[1]) for k, v in res['Fr'].items())
        line += ' | Drr_L1 %.4f | f_far %.4f (ex %.4f) L_far %.4f' % (
            np.mean([v[1] for v in res['Drr'].values()]), res['f_far'], res['f_far_ex'], res['L_far'])
        return line
    th, ph, rc = float(p['clump_th']), float(p['clump_ph']), float(p['clump_r'])
    xc = rc*np.array([math.sin(th)*math.cos(ph), math.sin(th)*math.sin(ph), math.cos(th)])
    chi = 1e-30*float(p['clump_rho'])
    B = float(p['arad'])*float(p['clump_temp'])**4
    B = 0.0 if B < 1e-60 else B
    clump = (xc, float(p['clump_a']), chi, B)
    ax = xc/rc
    ru = np.unique(np.round(r, 10))
    shells = [ru[np.argmin(np.abs(ru - x))] for x in (2.0, 3.6, 4.5, 5.5)]
    sel = np.isin(np.round(r, 10), shells) & (rh @ ax > math.cos(0.5))
    if B > 0:
        sel |= np.linalg.norm(r[:, None]*rh - xc[None, :], axis=1) < clump[1]
    idx = np.where(sel)[0]
    key = '%s_%s_%d_%.6f_%.6f_%g_%g' % ('cs' if cs else 'sp', p['nx1'], len(idx), th, ph, chi, B)
    os.makedirs(CACHE, exist_ok=True)
    cf = f'{CACHE}/{key}.npz'
    if os.path.exists(cf):
        z = np.load(cf)
        J, Hv = z['J'], z['H']
    else:
        J, Hv, K = exact.moments(r[idx, None]*rh[idx], clump)
        np.savez(cf, J=J, H=Hv)
    rs = r[idx]
    w = J0[idx]
    g = J/w
    Hr = np.einsum('pi,pi->p', Hv, rh[idx])
    Ei, Fi = E[idx], Fr[idx]
    R = {'all': np.ones(len(idx), bool), 'free': np.abs(g - 1) < 0.02, 'aff': np.abs(g - 1) >= 0.02,
         'aff_far': (np.abs(g - 1) >= 0.02) & (rs > 4.0), 'umbra': g < 0.05}
    if B > 0:
        R['inside'] = np.linalg.norm(rs[:, None]*rh[idx] - xc[None, :], axis=1) < clump[1]
    o = {}
    for nm, s in R.items():
        if s.sum() == 0:
            continue
        o[nm] = {'N': int(s.sum()), 'E_L1': float(np.mean(np.abs(Ei[s] - J[s])/w[s])),
                 'E_bias': float(np.mean((Ei[s] - J[s])/w[s])),
                 'Fr_L1': float(np.mean(np.abs(Fi[s] - Hr[s])*rs[s]**2)),
                 'E_over_J0': float(np.mean(Ei[s]/w[s])), 'Jex_over_J0': float(np.mean(J[s]/w[s]))}
    res['reg'] = o
    a = o['aff']
    u = o.get('umbra', {'N': 0, 'E_over_J0': float('nan'), 'Jex_over_J0': float('nan')})
    return ('all E_L1 %.4f | aff N %d E_L1 %.4f bias %+.4f Fr_L1 %.4f | aff_far E_L1 %.4f | umbra N %d E/J0 %.4f '
            '(ex %.4f)%s' % (o['all']['E_L1'], a['N'], a['E_L1'], a['E_bias'], a['Fr_L1'],
                             o.get('aff_far', {'E_L1': float('nan')})['E_L1'], u['N'], u['E_over_J0'],
                             u['Jex_over_J0'],
                             (' | inside E_L1 %.4f' % o['inside']['E_L1']) if 'inside' in o else ''))


def core(run, p, res):
    m = load(run, 'm1')
    t = load(run, 'm1_lat')
    r = m['r']
    rh = np.stack([t['m1_rx'], t['m1_ry'], t['m1_rz']], 1)
    chi = 1e-30*float(p['core_rho'])
    B = float(p['arad'])*float(p['core_temp'])**4
    ru = np.unique(np.round(r, 10))
    J, H, K = sphere_moments(ru, 1.0, chi, B)
    ii = np.searchsorted(ru, np.round(r, 10))
    Je, He, De = J[ii], H[ii], (K/J)[ii]
    E, Fr = m['m1_e'], m['m1_f1']
    reg = panel_region(rh, int(p['nx2']))
    out = r > 1.0
    res['E_out'] = stats((E - Je)[out], Je[out], reg[out], ['int', 'seam', 'corner'])
    res['E_in'] = stats((E - Je)[~out], Je[~out], reg[~out], ['int', 'seam', 'corner'])
    res['Fr_out'] = stats((Fr - He)[out], He[out], reg[out], ['int', 'seam', 'corner'])
    res['Drr_out'] = stats((t['m1_tt0'] - De)[out], De[out], reg[out], ['int', 'seam', 'corner'])
    Lr = []
    for x in (1.5, 3.0, 6.5):
        s = np.round(r, 10) == ru[np.argmin(np.abs(ru - x))]
        Lr.append((x, float(np.mean(Fr[s]*r[s]**2)), float(np.mean(He[s]*r[s]**2)),
                   float((Fr[s].max() - Fr[s].min())/np.mean(Fr[s]))))
    res['L'] = Lr
    return ('E_out L1 ' + ' '.join('%s %.4f' % (k, v[1]) for k, v in res['E_out'].items())
            + ' | E_in L1 %.4f' % np.mean([v[1] for v in res['E_in'].values()])
            + ' | Fr_out L1 ' + ' '.join('%s %.4f' % (k, v[1]) for k, v in res['Fr_out'].items())
            + ' | Drr_out L1 %.4f' % np.mean([v[1] for v in res['Drr_out'].values()])
            + ' | L/Lex ' + ' '.join('r%.1f %.4f (spread %.3f)' % (x, a/b, s) for x, a, b, s in Lr))


def main(run):
    p = cfg(run)
    res = {'run': run, 'cnt': counters(run)}
    nm = os.path.basename(run)
    try:
        if nm.startswith(('T1', 'T2', 'T3')):
            line = thin(run, p, res)
        elif nm.startswith('T4'):
            line = core(run, p, res)
        else:
            line = ''
    except Exception as e:  # noqa
        line = 'EVAL FAILED %r' % e
    json.dump(res, open(f'{run}/ev.json', 'w'), indent=1)
    c = res['cnt']
    print('%-20s cyc %s NA %s NP %s NC %s F %s st/BE %s/%s stfb %s fclip %s (%s) realis %s/%s floor %s/%s Pic %s | %s' % (
        nm, c['cyc'], c['notadm'], c['notpos'], c.get('ncv', '-'), c['fatal'], c.get('stage', '-'),
        c.get('be', '-'), c.get('stfb', '-'), c.get('fclip', '-'), c.get('fclipm', '-'), c.get('realis', '-'),
        c.get('drag', '-'), c.get('floorE', '-'), c.get('flooreint', '-'), c.get('picard', '-'), line))


if __name__ == '__main__':
    for d in sys.argv[1:]:
        main(d.rstrip('/'))
