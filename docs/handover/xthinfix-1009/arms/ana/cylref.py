"""fsanchor-1009 cyl test: exact moments and metrics.
Infinite cylinders along z: I(phi) is independent of the polar angle, so
  J = (1/2pi) int I dphi,   H = (1/8) int I (cos phi, sin phi) dphi   (E units, F = c H)
I(-phi): photons arriving FROM direction phi carry B = a T^4 if the ray toward phi first meets an emitter (absorbers
block, images in periodic y, x outside [0, 1] is vacuum), else 0.
usage: cylref.py <input> <rundir> [rundir ...]"""
import glob
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.environ['REPO'], 'vis/python'))
import bin_convert as bc  # noqa: E402

NPHI = 4096
LY = 1.0


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


def cyls(p):
    n = int(p['problem/cyl_n'])
    out = []
    for c in range(n):
        g = [float(p['problem/cyl%d_%s' % (c, q)]) for q in ('x', 'y', 'r', 'rho', 't')]
        out.append(g)
    return out


def exact(p, xs, ys):
    cy = cyls(p)
    global LY
    LY = float(p['mesh/x2max']) - float(p['mesh/x2min'])
    ar = float(p['rad_m1/arad'])
    tb = float(p['problem/cyl_t_bg'])
    phi = (np.arange(NPHI) + 0.5)*2*np.pi/NPHI
    ux, uy = np.cos(phi), np.sin(phi)
    kap = float(p['rad_m1/kappa_f'])
    xg, wt = np.polynomial.legendre.leggauss(48)
    th = 0.5*np.pi*(xg + 1.0)
    wt = 0.5*np.pi*wt
    st = np.sin(th)
    J = np.zeros((len(ys), len(xs)))
    Hx = np.zeros_like(J)
    Hy = np.zeros_like(J)
    inside = np.zeros(J.shape, bool)
    for jj, y in enumerate(ys):
        for ii, x in enumerate(xs):
            if any((x - c[0])**2 + (y - c[1])**2 <= c[2]**2 for c in cy):
                inside[jj, ii] = True
                continue
            # distance to leave the box in x
            with np.errstate(divide='ignore'):
                tx = np.where(ux > 0, (1.0 - x)/ux, np.where(ux < 0, -x/ux, np.inf))
            tbest = tx.copy()
            ival = np.zeros(NPHI)
            tauc = np.full(NPHI, np.inf)
            for c in cy:
                bri = ar*c[4]**4 if c[4] > 10*tb else 0.0
                chi = c[3]*kap
                for img in (-2, -1, 0, 1, 2):
                    dx, dy = c[0] - x, c[1] + img*LY - y
                    b = dx*ux + dy*uy
                    disc = b*b - (dx*dx + dy*dy - c[2]**2)
                    hit = (disc >= 0) & (b > 0)
                    t = np.where(hit, b - np.sqrt(np.maximum(disc, 0)), np.inf)
                    closer = t < tbest
                    tbest = np.where(closer, t, tbest)
                    ival = np.where(closer, bri, ival)
                    tauc = np.where(closer, chi*2*np.sqrt(np.maximum(disc, 0)), tauc)
            # polar integrals of a finite in-plane chord optical depth tau_xy
            ta = tauc[:, None]/st[None, :]
            ea = 1.0 - np.exp(-np.minimum(ta, 700.0))
            aw = 0.5*(ea*st[None, :]*wt[None, :]).sum(axis=1)       # -> 1
            bw = (ea*st[None, :]**2*wt[None, :]).sum(axis=1)        # -> pi/2
            dphi = 2*np.pi/NPHI
            J[jj, ii] = (ival*aw).sum()*dphi/(2*np.pi)
            Hx[jj, ii] = -(ival*bw*ux).sum()*dphi/(4*np.pi)
            Hy[jj, ii] = -(ival*bw*uy).sum()*dphi/(4*np.pi)
    return J, Hx, Hy, inside


def main():
    p = readin(sys.argv[1])
    ref = None
    for d in sys.argv[2:]:
        f = sorted(glob.glob(d + '/bin/*.m1.*.bin'))[-1]
        m = bc.read_binary_as_athdf(f)
        xs, ys = np.asarray(m['x1v']), np.asarray(m['x2v'])
        exf = d + '/../exact_%s.npz' % sys.argv[1].split('/')[-1].split('.')[0]
        if ref is None and glob.glob(exf):
            z = np.load(exf)
            X0, Y0 = np.meshgrid(xs, ys)
            ins0 = np.zeros(z['J'].shape, bool)
            for c in cyls(p):
                ins0 |= (X0 - c[0])**2 + (Y0 - c[1])**2 <= c[2]**2
            ref = (z['J'], z['Hx'], z['Hy'], ins0)
        if ref is None:
            ref = exact(p, xs, ys)
            np.savez(d + '/../exact_%s.npz' % sys.argv[1].split('/')[-1].split('.')[0],
                     J=ref[0], Hx=ref[1], Hy=ref[2], x=xs, y=ys)
        J, Hx, Hy, ins = ref
        E = m['m1_e'][0]
        cl = float(p['rad_m1/c_light'])
        for w in open(d + '/cmd.txt').read().split():
            if w.startswith('rad_m1/c_light='):
                cl = float(w.split('=')[1])
        # the TRANSPORTED face fluxes (m1_face), centred: F_x = mean of the two x1 faces,
        # F_y = mean of faces j and j+1 (periodic y, one block); the cell flux m1_f1 is
        # derived and clipped to |F| <= cE
        fc = bc.read_binary_as_athdf(sorted(glob.glob(d + '/bin/*.m1_face.*.bin'))[-1])
        F1 = 0.5*(fc['m1_fx1'][0] + fc['m1_fx1o'][0])/cl
        F2 = 0.5*(fc['m1_fx2'][0] + np.roll(fc['m1_fx2'][0], -1, axis=0))/cl
        C1, C2 = m['m1_f1'][0]/cl, m['m1_f2'][0]/cl
        X, Y = np.meshgrid(xs, ys)
        out = ~ins & (X > 0.3) & (X < 0.95)
        ej = E[out]/np.maximum(J[out], 1e-30)
        hm = np.hypot(Hx, Hy)
        ang_ex = np.degrees(np.arctan2(Hy, Hx))
        ang = np.degrees(np.arctan2(F2, F1))
        dang = np.abs((ang - ang_ex + 180) % 360 - 180)
        fm, fe = np.hypot(F1, F2)/np.maximum(E, 1e-300), hm/np.maximum(J, 1e-300)
        print('%s t=%.3e' % (d, m['Time']))
        print('  x in 0.3-0.95: <|E/J-1|> %.3f  max|E/J-1| %.3f  <|F|/|cH|> %.3f  '
              'median dangle %.1f deg  p90 %.1f  <f_M1> %.3f <f_ex> %.3f' % (
                  np.mean(np.abs(ej - 1)), np.max(np.abs(ej - 1)),
                  np.mean(np.hypot(F1, F2)[out]/np.maximum(hm[out], 1e-30)),
                  np.median(dang[out]), np.percentile(dang[out], 90),
                  np.mean(fm[out]), np.mean(fe[out])))
        cm = np.hypot(C1, C2)
        dangc = np.abs((np.degrees(np.arctan2(C2, C1)) - ang_ex + 180) % 360 - 180)
        print('  cell flux (derived, clipped): <|F|/|cH|> %.3f  median dangle %.1f  '
              'face max |F|/cE %.2f  frac(|F_face| > cE) %.3f' % (
                  np.mean(cm[out]/np.maximum(hm[out], 1e-30)), np.median(dangc[out]),
                  np.max(fm[out]), np.mean(fm[out] > 1.0)))
        # x-cuts: profile in y of E/J and F_x/(cH_x) at x = 0.6, 0.85
        for xc in (0.6, 0.85):
            i = int(np.argmin(np.abs(xs - xc)))
            jj = np.arange(0, len(ys), max(1, len(ys)//16))
            print('  x=%.2f y:   ' % xs[i] + ' '.join('%6.3f' % ys[j] for j in jj))
            print('    J_ex     ' + ' '.join('%6.3f' % J[j, i] for j in jj))
            print('    E        ' + ' '.join('%6.3f' % E[j, i] for j in jj))
            print('    ang_ex   ' + ' '.join('%6.1f' % ang_ex[j, i] for j in jj))
            print('    ang      ' + ' '.join('%6.1f' % ang[j, i] for j in jj))
        np.savez(d + '/fields.npz', E=E, F1=F1, F2=F2)


main()
