# test (C) analysis: usage anac.py <g3 run dir>; arms vs c01 (reference) and vs c03ns (no spot)
import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/sprhd_0926/wt/vis/python')
import bin_convert as bc
R = sys.argv[1]; gam = 5/3


def load(arm, n):
    h = bc.read_binary(glob.glob(f'{R}/{arm}/bin/*.hydro_w.{n:05d}.bin')[0])
    m = bc.read_binary(glob.glob(f'{R}/{arm}/bin/*.m1.{n:05d}.bin')[0])
    out = {}
    blocks = []
    for b in range(h['n_mbs']):
        g = h['mb_geometry'][b]
        blocks.append((tuple(np.round(g, 9)), b))
    blocks.sort()
    keys = ['dens', 'eint', 'velx', 'vely', 'velz']
    for k in keys:
        out[k] = np.concatenate([h['mb_data'][k][b].ravel() for _, b in blocks])
    out['E'] = np.concatenate([m['mb_data']['m1_e'][b].ravel() for _, b in blocks])
    dv = []
    for g, b in blocks:
        nx1, nx2, nx3 = h['nx1_out_mb'], h['nx2_out_mb'], h['nx3_out_mb']
        r = g[0] + (np.arange(nx1)+0.5)*(g[1]-g[0])/nx1
        th = g[2] + (np.arange(nx2)+0.5)*(g[3]-g[2])/nx2
        v = (r[None, None, :]**2*np.sin(th)[None, :, None]*np.ones((nx3, 1, 1))
             * (g[1]-g[0])/nx1*(g[3]-g[2])/nx2*(g[5]-g[4])/nx3)
        dv.append(v.ravel())
    out['dv'] = np.concatenate(dv)
    out['t'] = h['time']
    return out


arms = ['c01', 'c03', 'c09', 'c03v4', 'c09v4', 'c03g1']
for n in [2, 4, 6]:
    ns = load('c03ns', n)
    Tns = ns['eint']/ns['dens']
    res = {}
    for a in arms:
        try:
            s = load(a, n)
        except Exception:
            continue
        dT = (s['eint']/s['dens'])/Tns - 1.0
        de = (s['eint'] + s['E'] + 0.5*s['dens']*(s['velx']**2+s['vely']**2+s['velz']**2)
              - ns['eint'] - ns['E'])
        res[a] = (s, dT, de)
    if 'c01' not in res:
        print('no reference at', n); continue
    s0, dT0, _ = res['c01']
    print('t = %.2f' % s0['t'])
    for a, (s, dT, de) in res.items():
        dv = s['dv']
        eT = np.sum(abs(dT - dT0)*dv)/np.sum(abs(dT0)*dv)
        vv = np.sqrt((s['velx']-s0['velx'])**2+(s['vely']-s0['vely'])**2+(s['velz']-s0['velz'])**2)
        v0 = np.sqrt(s0['velx']**2+s0['vely']**2+s0['velz']**2)
        ev = np.sum(vv*s['dens']*dv)/np.sum(v0*s0['dens']*dv)
        print('  %-6s peak dT %.4e  excess energy %.5e  L1(dT - ref)/L1(dT_ref) %.3e  '
              'L1(rho|v - v_ref|)/L1(rho|v_ref|) %.3e  max|v| %.3e'
              % (a, dT.max(), np.sum(de*dv), eT, ev, np.sqrt(s['velx']**2+s['vely']**2+s['velz']**2).max()))
