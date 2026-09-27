# usage: order.py <g6 dir>: space and time convergence orders at t = 0.2 (final dump), L1
# volume-weighted, per variable (T = p/rho, rho, E, v = (v1, v2, v3) summed); space: raw fields
# and spot-minus-twin perturbations, fine restricted to coarse (2x2x2 volume average).
import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/sprhd_0926/wt/vis/python')
import bin_convert as bc
R = sys.argv[1]; gm1 = 2/3


def load(arm):
    fh = sorted(glob.glob(f'{R}/{arm}/bin/*.hydro_w.*.bin'))[-1]
    h = bc.read_binary(fh); m = bc.read_binary(fh.replace('hydro_w', 'm1'))
    geo = np.array(h['mb_geometry']); n1, n2, n3 = h['nx1_out_mb'], h['nx2_out_mb'], h['nx3_out_mb']
    x1 = sorted(set(geo[:, 0])); x2 = sorted(set(geo[:, 2])); x3 = sorted(set(geo[:, 4]))
    sh = (len(x3)*n3, len(x2)*n2, len(x1)*n1)
    out = {k: np.zeros(sh) for k in ['T', 'rho', 'E', 'v1', 'v2', 'v3']}
    for b in range(h['n_mbs']):
        g = geo[b]; i = x1.index(g[0]); j = x2.index(g[2]); k = x3.index(g[4])
        s = (slice(k*n3, (k+1)*n3), slice(j*n2, (j+1)*n2), slice(i*n1, (i+1)*n1))
        d = h['mb_data']['dens'][b]
        out['rho'][s] = d; out['T'][s] = gm1*h['mb_data']['eint'][b]/d
        out['E'][s] = m['mb_data']['m1_e'][b]
        out['v1'][s] = h['mb_data']['velx'][b]; out['v2'][s] = h['mb_data']['vely'][b]
        out['v3'][s] = h['mb_data']['velz'][b]
    r0, r1 = geo[:, 0].min(), geo[:, 1].max(); t0, t1 = geo[:, 2].min(), geo[:, 3].max()
    r = r0 + (np.arange(sh[2])+0.5)*(r1-r0)/sh[2]; th = t0 + (np.arange(sh[1])+0.5)*(t1-t0)/sh[1]
    out['dv'] = (r**2)[None, None, :]*np.sin(th)[None, :, None]*np.ones((sh[0], 1, 1))
    out['t'] = h['time']
    return out


def restrict(a, dv):
    w = a*dv
    s = lambda x: x.reshape(x.shape[0]//2, 2, x.shape[1]//2, 2, x.shape[2]//2, 2).sum(axis=(1, 3, 5))
    return s(w)/s(dv), s(dv)


VARS = ['T', 'rho', 'E', 'v']


def err(a, b, dv, norm):
    out = {}
    for v in ['T', 'rho', 'E']:
        out[v] = np.sum(abs(a[v]-b[v])*dv)/np.sum(abs(norm[v])*dv)
    out['v'] = sum(np.sum(abs(a[c]-b[c])*dv) for c in ['v1', 'v2', 'v3']) / \
        sum(np.sum(abs(norm[c])*dv) for c in ['v1', 'v2', 'v3'])
    return out


def rs(o):
    q = {}
    for v in ['T', 'rho', 'E', 'v1', 'v2', 'v3']:
        q[v], dvc = restrict(o[v], o['dv'])
    q['dv'] = dvc
    return q


def delta(a, b):
    return {k: (a[k]-b[k] if k not in ('dv', 't') else a[k]) for k in a}


print('=== SPACE (cfl 0.1, t = 0.2): relative L1 of level - restricted next finer level')
try:
    L = [load('s1'), load('t_c0.1'), load('s3')]; N = [load('s1ns'), load('s2ns'), load('s3ns')]
    print('t =', [x['t'] for x in L])
    for tag, F in (('raw fields (norm: field)', L), ('perturbation spot - twin (norm: finest pert.)',
                   [delta(L[i], N[i]) for i in range(3)])):
        e = []
        for i in range(2):
            fr = rs(F[i+1])
            norm = F[i] if tag.startswith('raw') else rs(rs(F[2])) if i == 0 else rs(F[2])
            e.append(err(F[i], fr, F[i]['dv'], norm))
        print(' ', tag)
        for v in VARS:
            print('    %-4s e(48-96) %.3e  e(96-192) %.3e  order %.2f' % (v, e[0][v], e[1][v],
                  np.log2(e[0][v]/e[1][v])))
except Exception as ex:
    print('space failed:', ex)
print('=== TIME (96x32x32, t = 0.2): relative L1 vs cfl 0.025 (every 1), norm = the perturbation')
ref = load('t_c0.025'); tw = load('s2ns'); pn = delta(ref, tw)
for pre, lab in (('t_c', 'vet_col_every 1'), ('t4_c', 'vet_col_every 4 + lag')):
    cs = ['0.8', '0.4', '0.2', '0.1']; E = []
    for c in cs:
        try:
            E.append(err(load(pre+c), ref, ref['dv'], pn))
        except Exception:
            E.append(None)
    print(' ', lab)
    for v in VARS:
        s = '    %-4s' % v
        for i, c in enumerate(cs):
            s += '  cfl %s %.3e' % (c, E[i][v]) if E[i] else '  cfl %s  n/a' % c
            if i > 0 and E[i] and E[i-1]:
                s += ' (p %.2f)' % np.log2(E[i-1][v]/E[i][v])
        print(s)
