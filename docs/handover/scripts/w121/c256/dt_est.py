"""Offline CFL estimate of a remap.dat (dhj-remap-h): min over cells of dx_d/(|v_d|+cs)
per direction, cs^2 = gam (gam-1) e/rho with gam = 1 + p/e from the source restart
(median per old radial cell, interpolated in r).  usage: dt_est.py REMAPDIR RST CFL"""
import sys, struct, numpy as np
sys.path.insert(0, 'wt/docs/handover/scripts')
import dhj_remap as d
rdir, rst, cfl = sys.argv[1], sys.argv[2], float(sys.argv[3])
par = d.parse_input(open(rdir + '/remap.athinput').read())
h = d.read_rst(rst)
ng = h['ng']; nx1o, r0o, r1o, co = d.grid_of(h['par'])
ro = d.centroids(d.edges(nx1o, r0o, r1o, co, 0))
s3 = slice(ng, -ng); act = slice(ng, ng + nx1o)
gam_o = 1 + np.median((h['p'][:, s3, s3, act]/h['eint'][:, s3, s3, act]).reshape(-1, nx1o), 0)
nx1, r0, r1, c = d.grid_of(par)
en = d.edges(nx1, r0, r1, c, 0); rc = d.centroids(en); dr = np.diff(en)
gam = np.interp(rc, ro, gam_o)
n2 = int(par['mesh']['nx2']); n3 = int(par['mesh']['nx3'])
f = open(rdir + '/remap.dat', 'rb'); f.read(8)
nmb, n1, b2, b3 = struct.unpack('<4i', f.read(16))
best = {1: (1e99,), 2: (1e99,), 3: (1e99,)}
for m in range(nmb):
    ll = np.frombuffer(f.read(20), '<i4')
    a = np.frombuffer(f.read(8*5*b3*b2*n1), '<f8').reshape(5, b3, b2, n1)
    rho, u1, u2, u3, e = a
    j0, k0, pn = ll[1]*b2, ll[2]*b3, ll[4]
    xe = np.tan(0.25*np.pi*(-1 + 2.0*np.arange(j0, j0 + b2 + 1)/n2))
    ye = np.tan(0.25*np.pi*(-1 + 2.0*np.arange(k0, k0 + b3 + 1)/n3))
    X = np.tan(0.25*np.pi*(-1 + 2.0*(np.arange(j0, j0 + b2) + 0.5)/n2))[None, :]
    Y = np.tan(0.25*np.pi*(-1 + 2.0*(np.arange(k0, k0 + b3) + 0.5)/n3))[:, None]
    xl, xr, yl, yr = xe[None, :-1], xe[None, 1:], ye[:-1, None], ye[1:, None]
    dxi = np.arccos((1 + xl*xr + Y*Y)/np.sqrt(1 + xl*xl + Y*Y)/np.sqrt(1 + xr*xr + Y*Y))
    det = np.arccos((1 + X*X + yl*yr)/np.sqrt(1 + X*X + yl*yl)/np.sqrt(1 + X*X + yr*yr))
    C, D = np.sqrt(1 + X*X), np.sqrt(1 + Y*Y); cc = (-X*Y/(C*D))[..., None]
    v2 = (u2 - cc*u3)/(1 - cc*cc); v3 = (u3 - cc*u2)/(1 - cc*cc)
    cs = np.sqrt(gam*(gam - 1)*e/rho)
    for dd, t in ((1, dr/(np.abs(u1) + cs)), (2, rc*dxi[..., None]/(np.abs(v2) + cs)),
                  (3, rc*det[..., None]/(np.abs(v3) + cs))):
        i = np.unravel_index(np.argmin(t), t.shape)
        if t[i] < best[dd][0]:
            best[dd] = (t[i], pn, k0 + i[0], j0 + i[1], i[2], rc[i[2]], cs[i], rho[i],
                        [u1[i], v2[i], v3[i]][dd - 1])
for dd in (1, 2, 3):
    b = best[dd]
    print('dir %d: cfl*dt = %.3f s  panel %d k %d j %d i %d r %.4e cs %.3e rho %.2e v %.3e'
          % (dd, cfl*b[0], *b[1:]))
