"""blendall-1009 collimated-beam tests (xb20, ba0, ba20): exact moments and metrics.
Every object is OPAQUE (emitter T > 10 T_bg: I = a T^4; absorber: I = 0); infinite along z,
so J = (1/2pi) int I dphi and H = F/c = (1/8) int I Omega dphi (E units), as cylref.py
with tau -> inf.  In-plane rays from the cell centre back to the first object hit, periodic
images in y, vacuum outside x1min..x1max.  Region: x in [xlo, 0.95] (beyond the channel
mouths).  usage: collref.py <input> <rundir> [rundir ...]"""
import glob
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.environ['REPO'], 'vis/python'))
import bin_convert as bc  # noqa: E402

NPHI = 4096
XLO = 0.45


def readin(fn):
    p = {}
    blk = ''
    for ln in open(fn):
        ln = ln.split('#')[0].strip()
        if ln.startswith('<'):
            blk = ln.strip('<>')
        elif '=' in ln:
            k, v = [s.strip() for s in ln.split('=', 1)]
            p[blk + '/' + k] = v
    return p


def objects(p):
    out = []
    for c in range(int(p['problem/cyl_n'])):
        b = 'problem/cyl%d_' % c
        o = dict(x=float(p[b + 'x']), y=float(p[b + 'y']), r=float(p[b + 'r']),
                 t=float(p[b + 't']), rect=p.get(b + 'shape', 'cyl') == 'rect')
        if o['rect']:
            o.update(hl=float(p[b + 'hl']), hw=float(p[b + 'hw']),
                     ang=math.radians(float(p[b + 'ang'])))
        out.append(o)
    return out


def hit(o, x, y, ux, uy):
    """distance along (ux, uy) from (x, y) to object o (inf if missed)"""
    if not o['rect']:
        dx, dy = o['x'] - x, o['y'] - y
        b = dx*ux + dy*uy
        disc = b*b - (dx*dx + dy*dy - o['r']**2)
        ok = (disc >= 0) & (b > 0)
        return np.where(ok, b - np.sqrt(np.maximum(disc, 0)), np.inf)
    ca, sa = math.cos(o['ang']), math.sin(o['ang'])
    px, py = (x - o['x'])*ca + (y - o['y'])*sa, -(x - o['x'])*sa + (y - o['y'])*ca
    vx, vy = ux*ca + uy*sa, -ux*sa + uy*ca
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = (-o['hl'] - px)/vx
        t2 = (o['hl'] - px)/vx
        t3 = (-o['hw'] - py)/vy
        t4 = (o['hw'] - py)/vy
    tmin = np.maximum(np.minimum(t1, t2), np.minimum(t3, t4))
    tmax = np.minimum(np.maximum(t1, t2), np.maximum(t3, t4))
    tmin = np.where(np.isnan(tmin), -np.inf, tmin)
    ok = (tmax >= np.maximum(tmin, 0)) & (tmax > 0)
    return np.where(ok, np.maximum(tmin, 0), np.inf)


def exact(p, xs, ys):
    ob = objects(p)
    ar = float(p['rad_m1/arad'])
    tb = float(p['problem/cyl_t_bg'])
    x0, x1 = float(p['mesh/x1min']), float(p['mesh/x1max'])
    ly = float(p['mesh/x2max']) - float(p['mesh/x2min'])
    phi = (np.arange(NPHI) + 0.5)*2*np.pi/NPHI
    ux, uy = np.cos(phi), np.sin(phi)
    J = np.zeros((len(ys), len(xs)))
    Hx = np.zeros_like(J)
    Hy = np.zeros_like(J)
    for jj, y in enumerate(ys):
        for ii, x in enumerate(xs):
            with np.errstate(divide='ignore'):
                tbest = np.where(ux > 0, (x1 - x)/ux, np.where(ux < 0, (x0 - x)/ux, np.inf))
            ival = np.zeros(NPHI)
            for o in ob:
                bri = ar*o['t']**4 if o['t'] > 10*tb else 0.0
                for img in (-1, 0, 1):
                    oo = dict(o)
                    oo['y'] = o['y'] + img*ly
                    t = hit(oo, x, y, ux, uy)
                    closer = t < tbest
                    tbest = np.where(closer, t, tbest)
                    ival = np.where(closer, bri, ival)
            dphi = 2*np.pi/NPHI
            J[jj, ii] = ival.sum()*dphi/(2*np.pi)
            Hx[jj, ii] = -(ival*ux).sum()*dphi/8.0
            Hy[jj, ii] = -(ival*uy).sum()*dphi/8.0
    return J, Hx, Hy


def main():
    inp = sys.argv[1]
    p = readin(inp)
    name = os.path.basename(inp).split('.')[0]
    for d in sys.argv[2:]:
        m = bc.read_binary_as_athdf(sorted(glob.glob(d + '/bin/*.m1.*.bin'))[-1])
        xs, ys = np.asarray(m['x1v']), np.asarray(m['x2v'])
        sel = (xs >= XLO) & (xs <= 0.95)
        exf = os.path.join(os.path.dirname(d.rstrip('/')), 'exact_%s.npz' % name)
        if os.path.exists(exf):
            z = np.load(exf)
            J, Hx, Hy = z['J'], z['Hx'], z['Hy']
        else:
            J, Hx, Hy = exact(p, xs[sel], ys)
            np.savez(exf, J=J, Hx=Hx, Hy=Hy, x=xs[sel], y=ys)
        cl = float(p['rad_m1/c_light'])
        for w in open(d + '/cmd.txt').read().split():
            if w.startswith('rad_m1/c_light='):
                cl = float(w.split('=')[1])
        fc = bc.read_binary_as_athdf(sorted(glob.glob(d + '/bin/*.m1_face.*.bin'))[-1])
        E = m['m1_e'][0][:, sel]
        F1 = (0.5*(fc['m1_fx1'][0] + fc['m1_fx1o'][0])/cl)[:, sel]
        F2 = (0.5*(fc['m1_fx2'][0] + np.roll(fc['m1_fx2'][0], -1, axis=0))/cl)[:, sel]
        # face-level |F|/cE: x1 and x2 faces separately (the superluminal test)
        Ex = m['m1_e'][0]
        # against the larger of the two cells' E: an upwind flux c f E_up, |f| <= 1, never
        # exceeds it; the x1 face i is between cells i-1 and i (the i = 0 face: cell 0)
        exm = np.concatenate([Ex[:, :1], Ex[:, :-1]], axis=1)
        fx1 = np.abs(fc['m1_fx1'][0])/(cl*np.maximum(np.maximum(Ex, exm), 1e-300))
        fx2 = np.abs(fc['m1_fx2'][0])/(cl*np.maximum(np.maximum(Ex, np.roll(Ex, 1, axis=0)),
                                                   1e-300))
        # faces whose larger cell E is at the floor (inside an absorber: E solved
        # negative and raised to e_floor) are excluded (their ratio is 1e28)
        emx = Ex.max()
        fx1 = np.where(np.maximum(Ex, exm) > 1e-8*emx, fx1, 0.0)
        fx2 = np.where(np.maximum(Ex, np.roll(Ex, 1, axis=0)) > 1e-8*emx, fx2, 0.0)
        fface = np.maximum(fx1, fx2)[:, sel]
        hm = np.hypot(Hx, Hy)
        jmax = J.max()
        beam = J > 0.05*jmax
        dark = J < 1.0e-3*jmax
        ang_ex = np.degrees(np.arctan2(Hy, Hx))
        dang = np.abs((np.degrees(np.arctan2(F2, F1)) - ang_ex + 180) % 360 - 180)
        fm = np.hypot(F1, F2)
        wgt = J[beam]/J[beam].sum()
        print('%s t=%.3e  [%s, x in %.2f..0.95, beam = J > 0.05 Jmax: %d cells]' % (
            d, m['Time'], name, XLO, beam.sum()))
        print('  beam: E/J median %.3f  sum E/sum J %.3f  |F|/|cH| median %.3f  '
              'flux-weighted <|F|>/<|cH|> %.3f' % (
                  np.median(E[beam]/J[beam]), E[beam].sum()/J[beam].sum(),
                  np.median(fm[beam]/hm[beam]), fm[beam].sum()/hm[beam].sum()))
        print('  beam: direction error J-weighted mean %.1f deg  median %.1f  p90 %.1f   '
              'face |F|/cE > 1: frac %.4f (beam) %.4f (all)  max %.2f' % (
                  (wgt*dang[beam]).sum(), np.median(dang[beam]),
                  np.percentile(dang[beam], 90), np.mean(fface[beam] > 1.0),
                  np.mean(fface > 1.0), fface.max()))
        print('  dark (J < 1e-3 Jmax, %d cells): <E>/Jmax %.4f  max E/Jmax %.4f' % (
            dark.sum(), E[dark].mean()/jmax if dark.any() else np.nan,
            E[dark].max()/jmax if dark.any() else np.nan))
        xsel = xs[sel]
        for xc in (0.70, 0.90):
            i = int(np.argmin(np.abs(xsel - xc)))
            jl = J[:, i]
            e = E[:, i]
            pk = jl > 0.5*jl.max()
            en, jn = e/max(e.sum(), 1e-300), jl/jl.sum()
            print('  x=%.2f shape: L1 |E/sumE - J/sumJ| %.3f' % (xsel[i], np.abs(en - jn).sum()))
            print('  x=%.2f profile: overlap sum min(E,J)/sum J %.3f  L1 |E-J|/sum J %.3f  '
                  'E-centroid %.4f (exact %.4f)  E outside the J>5%% band / total %.3f' % (
                      xsel[i], np.minimum(e, jl).sum()/jl.sum(),
                      np.abs(e - jl).sum()/jl.sum(), (e*ys).sum()/e.sum(),
                      (jl*ys).sum()/jl.sum(), e[jl < 0.05*jl.max()].sum()/e.sum()))
            # peaks / edges: list the y of local maxima of E and of J above 30 % max
            def peaks(v):
                k = [q for q in range(1, len(v)-1) if v[q] >= v[q-1] and v[q] > v[q+1]
                     and v[q] > 0.3*v.max()]
                return ' '.join('%.3f' % ys[q] for q in k)
            dy = ys[1] - ys[0]
            print('    peaks y: E [%s]  J [%s]   E max/J max %.3f  width(>50%% max) E %.3f '
                  'J %.3f' % (peaks(e), peaks(jl), e.max()/jl.max(),
                              (e > 0.5*e.max()).sum()*dy, pk.sum()*dy))
        np.savez(d + '/fields_coll.npz', E=E, F1=F1, F2=F2, x=xsel, y=ys)


main()
