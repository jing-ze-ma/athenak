#!/usr/bin/env python3
"""lowbeta-tests-1102: per-run metrics -> metrics.json (numpy only).  usage: analyze.py <dir>"""
import glob, json, math, os, re, sys
import numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/lowbeta_1102/bip/src/vis/python')
import bin_convert as bc

d = sys.argv[1]
meta = json.load(open(f'{d}/meta.json'))
st = json.load(open(f'{d}/status.json')) if os.path.exists(f'{d}/status.json') else {}
test = meta['test']
deck = open(f'{d}/run.athinput').read()
gam = float(re.search(r'^gamma = (\S+)', deck, re.M).group(1))
gm1 = gam - 1.0
M = dict(meta)
M.update(st)


def rd(f):
    x = bc.read_binary(f)
    out = {k: v[0] for k, v in x['mb_data'].items()}   # single MeshBlock
    out['time'] = x['time']
    out['x1min'], out['x1max'] = x['x1min'], x['x1max']
    out['x2min'], out['x2max'] = x['x2min'], x['x2max']
    return out


def squeeze(a):
    return a[0] if a.shape[0] == 1 else a


# ---- event log: total floor/limiter counts and first cycle with a floor ----
cnt = dict(dfloor=0, efloor=0, tfloor=0, vceil=0, fail=0, fofc=0)
first = None
if os.path.exists(f'{d}/run.log'):
    for ln in open(f'{d}/run.log'):
        if ln.startswith('#') or not ln.strip():
            continue
        v = ln.split()
        c = int(v[0])
        n = [int(v[1]), int(v[2]), int(v[3]), int(v[4]), int(v[5]), int(v[7])]
        for k, x in zip(cnt, n):
            cnt[k] += x
        if first is None and (n[0] or n[1] or n[2]):
            first = c
M['counts'] = cnt
M['first_floor_cycle'] = first

# ---- out.log: final time / cycles / crash ----
out = open(f'{d}/out.log', errors='replace').read() if os.path.exists(f'{d}/out.log') else ''
m = re.findall(r'cycle=(\d+) time=(\S+) dt=(\S+)', out)
if m:
    M['last_cycle'], M['last_time'], M['last_dt'] = int(m[-1][0]), float(m[-1][1]), float(m[-1][2])
M['fatal'] = ('FATAL' in out) or ('non-finite' in out) or ('Segmentation' in out)
mt = re.search(r'time=(\S+) cycle=(\d+)\s*\n.*tlim', out)
fin = re.findall(r'^time=(\S+) cycle=(\d+)', out, re.M)
if fin:
    M['end_time'] = float(fin[-1][0]); M['end_cycle'] = int(fin[-1][1])
M['completed'] = ('Terminating on time limit' in out) or ('time limit' in out.lower() and not M['fatal'])

# ---- history ----
hst = sorted(glob.glob(f'{d}/run.*hst'))
H = None
if hst:
    hd = None
    for f in hst:
        lines = open(f).read().splitlines()
        hdr = [l for l in lines if l.startswith('#  [')]
        labs = re.findall(r'\[\d+\]=(\S+)', hdr[0]) if hdr else []
        try:
            a = np.loadtxt(f, ndmin=2)
        except Exception:
            a = np.zeros((0, len(labs)))
        if a.size:
            H = H or {}
            for i, l in enumerate(labs):
                H[l] = a[:, i]
if H is not None and 'pmin' in H:
    M['hist'] = {k: H[k].tolist() for k in ('time', 'pmin', 'emin', 'bmin', 'dmin') if k in H}
    M['pmin_all'] = float(np.min(H['pmin']))
    M['emin_all'] = float(np.min(H['emin']))
    M['hist_tmax'] = float(H['time'][-1])

files = sorted(glob.glob(f'{d}/bin/run.mhd_w_bcc.*.bin'))
pfiles = sorted(glob.glob(f'{d}/bin/run.mhd_w_e.*.bin'))


def pseries(p0=None):
    ts, mins, l1s = [], [], []
    for f in pfiles:
        try:
            x = rd(f)
        except Exception:
            continue
        p = gm1*x['eint']
        ts.append(x['time']); mins.append(float(p.min()))
        if p0 is not None:
            l1s.append(float(np.mean(np.abs(p - p0))))
    return ts, mins, l1s


if test.startswith('cpaw'):
    p0 = meta['pres']
    ts, mins, l1s = pseries(p0)
    M['p_t'] = ts; M['pmin_t'] = mins; M['pL1_t'] = l1s
    if len(files) >= 2:
        a, b = rd(files[0]), rd(files[-1])
        M['t_final'] = b['time']
        if meta['dim'] == 2:
            k = np.array([1/math.sqrt(5), 2/math.sqrt(5), 0.0])
        else:
            k = np.array([1.0, 0.0, 0.0])
        def perp(x, n):
            v = np.stack([x[n[0]], x[n[1]], x[n[2]]])
            par = np.tensordot(k, v, axes=1)
            return v - k[:, None, None, None]*par, par
        bp0, bpar0 = perp(a, ('bcc1', 'bcc2', 'bcc3'))
        bp1, bpar1 = perp(b, ('bcc1', 'bcc2', 'bcc3'))
        vp0, _ = perp(a, ('velx', 'vely', 'velz'))
        vp1, _ = perp(b, ('velx', 'vely', 'velz'))
        M['L1_Bperp'] = float(np.mean(np.sqrt(((bp1 - bp0)**2).sum(0))))
        M['L1_vperp'] = float(np.mean(np.sqrt(((vp1 - vp0)**2).sum(0))))
        p1 = gm1*b['eint']
        M['L1_p'] = float(np.mean(np.abs(p1 - p0)))
        M['Linf_p'] = float(np.max(np.abs(p1 - p0)))
        M['p0_ic_check'] = float(np.mean(gm1*a['eint']))
        M['L1_rho'] = float(np.mean(np.abs(b['dens'] - 1.0)))
        bm1 = np.sqrt(b['bcc1']**2 + b['bcc2']**2 + b['bcc3']**2)
        bm0 = math.sqrt(1.0 + 0.01)
        M['L1_absB'] = float(np.mean(np.abs(bm1 - bm0)))
        M['pmin_final'] = float(p1.min())

elif test == 'vortex':
    ts, mins, _ = pseries()
    M['p_t'] = ts; M['pmin_t'] = mins
    if len(files) >= 2:
        b = rd(files[-1])
        t = b['time']
        vt, bt = meta['vtil'], meta['btil']
        va = vt/math.sqrt(2.0)
        nx = b['dens'].shape[-1]; ny = b['dens'].shape[-2]
        dx = (b['x1max'] - b['x1min'])/nx; dy = (b['x2max'] - b['x2min'])/ny
        xc = b['x1min'] + (np.arange(nx) + 0.5)*dx
        yc = b['x2min'] + (np.arange(ny) + 0.5)*dy
        gp, gw = np.polynomial.legendre.leggauss(4)
        L = b['x1max'] - b['x1min']
        acc = {k: np.zeros((ny, nx)) for k in ('vx', 'vy', 'bx', 'by', 'p')}
        for ia in range(4):
            for ib in range(4):
                X = xc[None, :] + 0.5*dx*gp[ia]
                Y = yc[:, None] + 0.5*dy*gp[ib]
                Xs = (X - va*t - b['x1min']) % L + b['x1min']
                Ys = (Y - va*t - b['x2min']) % L + b['x2min']
                r2 = Xs**2 + Ys**2
                e = np.exp(0.5*(1 - r2))
                w = gw[ia]*gw[ib]/4
                acc['vx'] += w*(va - vt*e*Ys); acc['vy'] += w*(va + vt*e*Xs)
                acc['bx'] += w*(-bt*e*Ys); acc['by'] += w*(bt*e*Xs)
                acc['p'] += w*(1.0 + (0.5*bt*bt*(1 - r2) - 0.5*vt*vt)*e*e)
        sq = lambda a: a.reshape(ny, nx)
        p1 = gm1*sq(b['eint'])
        M['t_final'] = t
        M['L1_p'] = float(np.mean(np.abs(p1 - acc['p'])))
        M['L1_v'] = float(np.mean(np.abs(sq(b['velx']) - acc['vx']) + np.abs(sq(b['vely']) - acc['vy'])))
        M['L1_B'] = float(np.mean(np.abs(sq(b['bcc1']) - acc['bx']) + np.abs(sq(b['bcc2']) - acc['by'])))
        M['pmin_final'] = float(p1.min())
        M['pmin_exact_cellavg'] = float(acc['p'].min())

elif test == 'lw':
    f = f'{d}/run-errs.dat'
    if os.path.exists(f):
        v = [l for l in open(f) if not l.startswith('#')][-1].split()
        names = ['nx1', 'nx2', 'nx3', 'ncycle', 'rms', 'linf', 'd', 'M1', 'M2', 'M3', 'E', 'B1', 'B2', 'B3']
        for n, x in zip(names, v):
            M['err_' + n] = float(x)

elif test == 'rd':
    if files:
        b = rd(files[-1])
        x = b['x1min'] + (np.arange(b['dens'].size) + 0.5)*(b['x1max'] - b['x1min'])/b['dens'].size
        rho = b['dens'].ravel(); p = gm1*b['eint'].ravel()
        M['t_final'] = b['time']
        M['rho_min'] = float(rho.min()); M['rho_max'] = float(rho.max())
        M['eint_max'] = float(b['eint'].max()); M['eint_min'] = float(b['eint'].min())
        np.savez_compressed(f'{d}/profile.npz', x=x, rho=rho, eint=b['eint'].ravel(),
                            by=b['bcc2'].ravel(), vx=b['velx'].ravel(), vy=b['vely'].ravel())

elif test in ('rotor', 'blast'):
    if len(files) >= 2:
        b = rd(files[-1])
        M['t_final'] = b['time']
        p1 = gm1*b['eint']
        M['pmin_final'] = float(p1.min())
        nfl = int(np.sum(p1 < 1e-30))
        M['ncell_p_floor_final'] = nfl
        if b['dens'].shape[0] == 1:
            sl = p1[0]; dn = b['dens'][0]
        else:
            k = b['dens'].shape[0]//2
            sl = p1[k]; dn = b['dens'][k]
        np.savez_compressed(f'{d}/final.npz', p=sl.astype(np.float32), d=dn.astype(np.float32))

json.dump(M, open(f'{d}/metrics.json', 'w'))
