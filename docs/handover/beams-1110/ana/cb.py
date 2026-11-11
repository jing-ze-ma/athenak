"""littests-1011 crossing beams (Ma+2025 4.1.1 / J21 4.2), rotated: x1 = paper y, x2 = paper x.
usage: cb.py <rundir> ...   prints metrics; exact = cell-averaged analytic straight beams."""
import sys
import numpy as np
from lcommon import load, runinfo, W

L2 = 1.0


def readin(fn):
    p, blk = {}, ''
    for ln in open(fn):
        ln = ln.split('#')[0].strip()
        if ln.startswith('<'):
            blk = ln.strip('<>')
        elif '=' in ln:
            k, v = [s.strip() for s in ln.split('=', 1)]
            p[blk + '/' + k] = v
    return p


def beams(p):
    out = []
    for b in range(int(p.get('rad_m1/implicit_beam_n', 0))):
        q = 'rad_m1/implicit_beam%d_' % b
        out.append(dict(e=float(p[q + 'e']), ang=np.radians(float(p[q + 'ang'])),
                        y0=float(p[q + 'y0']), y1=float(p[q + 'y1'])))
    return out


def exact(p, x1, x2, each=False, ns=32, nsy=32):
    """cell average of straight periodic strips (x1 sub-samples x x2 sub-samples)"""
    x1a, x2a = float(p['mesh/x1min']), float(p['mesh/x2min'])
    L = float(p['mesh/x2max']) - x2a
    dx1, dx2 = x1[1] - x1[0], x2[1] - x2[0]
    xs = (x1[:, None] - 0.5*dx1 + (np.arange(ns)[None, :] + 0.5)*dx1/ns)      # (nx1, ns)
    ys = (x2[:, None] - 0.5*dx2 + (np.arange(nsy)[None, :] + 0.5)*dx2/nsy)    # (nx2, nsy)
    res = []
    for b in beams(p):
        s = (xs - x1a)*np.tan(b['ang'])                                       # (nx1, ns)
        yy = ys[:, :, None, None] - s[None, None, :, :]                         # (nx2,nsy,nx1,ns)
        yy = (yy - b['y0']) % L
        inside = yy < (b['y1'] - b['y0'])
        res.append(b['e']*inside.mean(axis=(1, 3)))
    return res if each else sum(res)


def centre(p, b, x1):
    x1a, x2a = float(p['mesh/x1min']), float(p['mesh/x2min'])
    L = float(p['mesh/x2max']) - x2a
    c = 0.5*(b['y0'] + b['y1']) + (x1 - x1a)*np.tan(b['ang'])
    return (c - x2a) % L + x2a


def fwhm(x, v):
    m = v.max()
    k = np.argmax(v)
    h = 0.5*m
    lo = k
    while lo > 0 and v[lo] > h:
        lo -= 1
    hi = k
    while hi < len(v) - 1 and v[hi] > h:
        hi += 1
    xl = x[lo] + (h - v[lo])*(x[lo+1] - x[lo])/(v[lo+1] - v[lo])
    xh = x[hi-1] + (h - v[hi-1])*(x[hi] - x[hi-1])/(v[hi] - v[hi-1])
    return xh - xl


def metrics(d, rows=(-1.35, -0.35, 0.65, 1.65)):
    p = readin(d + '/in.athinput')
    v = load(d)
    x1, x2 = v['x1'], v['x2']
    dx2 = x2[1] - x2[0]
    L = x2[-1] - x2[0] + dx2
    Eex = exact(p, x1, x2)
    bm = beams(p)
    out = {'run': d.split('runs/')[-1]}
    out.update(runinfo(d))
    for nm in ('E', 'J'):
        a = v.get(nm)
        if a is None:
            continue
        out[nm + '_L1'] = np.abs(a - Eex).sum()/Eex.sum()
        # row luminosity (energy per row) at the top relative to exact
        out[nm + '_lum_top'] = a[:, -3].sum()/Eex[:, -3].sum()
        ws, pk, off = [], [], []
        for r in rows:
            i = int(np.argmin(np.abs(x1 - r)))
            b = bm[1]
            c = centre(p, b, x1[i])
            dd = (x2 - c + 0.5*L) % L - 0.5*L
            s = np.abs(dd) < 0.24
            o = np.argsort(dd[s])
            xx, vv = dd[s][o], a[s, i][o]
            ws.append(fwhm(xx, vv)/dx2)
            pk.append(vv.max()/b['e'])
            off.append((xx*vv).sum()/vv.sum()/dx2)
        out[nm + '_fwhm'] = ws
        out[nm + '_peak'] = pk
        out[nm + '_off'] = off
    # exact FWHM / peak at the same rows (cell averaging)
    ws, pk = [], []
    for r in rows:
        i = int(np.argmin(np.abs(x1 - r)))
        c = centre(p, bm[1], x1[i])
        dd = (x2 - c + 0.5*L) % L - 0.5*L
        s = np.abs(dd) < 0.24
        o = np.argsort(dd[s])
        ws.append(fwhm(dd[s][o], Eex[s, i][o])/dx2)
        pk.append(Eex[s, i].max()/bm[1]['e'])
    out['ex_fwhm'] = ws
    out['ex_peak'] = pk
    out['minE'] = float(v['E'].min())
    return out


if __name__ == '__main__':
    for d in sys.argv[1:]:
        m = metrics(d)
        f = lambda a: ' '.join('%.2f' % q for q in a)  # noqa: E731
        print('%-22s L1 E %.3f J %.3f | lum_top E %.3f J %.3f | FWHM(cells) E %s J %s ex %s | '
              'peak E %s J %s ex %s | off E %s | picard %.1f nc %d %.0fs' % (
                  m['run'], m['E_L1'], m['J_L1'], m['E_lum_top'], m['J_lum_top'],
                  f(m['E_fwhm']), f(m['J_fwhm']), f(m['ex_fwhm']), f(m['E_peak']),
                  f(m['J_peak']), f(m['ex_peak']), f(m['E_off']), m.get('picard', -1),
                  m.get('nonconv', -1), m.get('sec', -1)))
