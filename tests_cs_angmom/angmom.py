"""Axial angular-momentum budget of the deep-hot-Jupiter production runs, cs vs sp.

Reads the mhd_w_bcc binary dumps (full 3-D, single precision) and forms

    L_z(inertial) = int rho * [ u_zonal * R + Omega * R^2 ] dV ,   R = r cos(lat)

    L_rel   = int rho u_zonal R dV          (relative-flow term)
    L_frame = Omega int rho R^2 dV          (mass-redistribution term)

Velocity bases (checked in src/):
  cs  w0 (IVX,IVY,IVZ) = CONTRAVARIANT components on the gnomonic tangent basis of
      UNIT but NON-ORTHOGONAL vectors (e_r, e_xi, e_eta); see
      Coordinates::GnomonicEquiangleRaiseVel.  Physical v = v1 e_r + v2 e_xi + v3 e_eta,
      so projecting v2,v3 on the unit zonal direction (dhjcs.winds) is exact.
  sp  w0 = orthonormal (v_r, v_theta, v_phi); zonal = velz.

Writes angmom.npz with one row per dump.
"""
import sys
import glob
import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhjcs  # noqa: E402

B = '/viper/u2/jinma/ATHENAK/bench/'
OMEGA = 2.06e-5
PROT = 2.0*np.pi/OMEGA          # 3.0501e5 s
CPOLY = [-0.068392, -2.191487, 2.464818, -1.366698]
FSTR_TH = 3.0


def stretch_r(r, r0, r1):
    xi = (r-r0)/(r1-r0)
    u = xi.copy()
    xik = xi.copy()
    for k in range(1, 5):
        u = u + CPOLY[k-1]*xik*(1-xi)
        xik = xik*xi
    return r0 + (r1-r0)*u


def radial_moments(rf):
    """volume element and volume-weighted <r>, <r^2> for each radial cell."""
    rl, rr = rf[:-1], rf[1:]
    dV = (rr**3 - rl**3)/3.0
    r1 = 0.25*(rr**4 - rl**4)/dV/3.0*3.0/1.0
    r1 = (rr**4 - rl**4)/4.0/dV
    r2 = (rr**5 - rl**5)/5.0/dV
    return dV, r1, r2


def geom(raw, kind):
    """returns dV (nmb,nk,nj,ni), coslat (nmb,nk,nj), rmom1, rmom2 (ni,), rc (ni,)"""
    g = np.asarray(raw['mb_geometry'])
    nmb = g.shape[0]
    ni = raw['nx1_mb']
    nj = raw['nx2_mb']
    nk = raw['nx3_mb']
    r0, r1 = raw['x1min'], raw['x1max']
    A = np.zeros((nmb, nk, nj))
    if kind == 'cs':
        gg = dhjcs.cs_geometry(raw)
        lat = gg['lat']
        # EXACT solid angle of an equiangular gnomonic cell:
        #   S = F(xi2,eta2) - F(xi1,eta2) - F(xi2,eta1) + F(xi1,eta1),
        #   F(xi,eta) = atan( tan(xi) tan(eta) / sqrt(1 + tan^2 xi + tan^2 eta) )
        def F(a, b):
            X = np.tan(a)
            Y = np.tan(b)
            return np.arctan(X*Y/np.sqrt(1.0+X*X+Y*Y))
        for m in range(nmb):
            xif = np.pi/4*(g[m, 2]+(g[m, 3]-g[m, 2])*np.arange(nj+1)/nj)
            etaf = np.pi/4*(g[m, 4]+(g[m, 5]-g[m, 4])*np.arange(nk+1)/nk)
            XI, ETA = np.meshgrid(xif, etaf)
            Ff = F(XI, ETA)
            A[m] = Ff[1:, 1:] - Ff[1:, :-1] - Ff[:-1, 1:] + Ff[:-1, :-1]
    else:
        gg = dhjcs.sp_geometry(raw, FSTR_TH)
        lat = gg['lat']
        for m in range(nmb):
            thf = g[m, 2]+(g[m, 3]-g[m, 2])*np.arange(nj+1)/nj
            xi = thf/np.pi
            thf = np.pi/2*(1+np.sinh(FSTR_TH*(2*xi-1))/np.sinh(FSTR_TH))
            dph = (g[m, 5]-g[m, 4])/nk
            A[m] = np.tile((np.cos(thf[:-1])-np.cos(thf[1:]))*dph, (nk, 1))
    rf = stretch_r(g[0, 0]+(g[0, 1]-g[0, 0])*np.arange(ni+1)/ni, r0, r1)
    dVr, rm1, rm2 = radial_moments(rf)
    dV = A[..., None]*dVr[None, None, None, :]
    return dV, np.cos(lat), rm1, rm2, 0.5*(rf[1:]+rf[:-1]), rf


def one(path, kind, cache={}):
    raw = dhjcs.bin_convert.read_binary(path)
    key = (kind, raw['mb_geometry'].shape if hasattr(raw['mb_geometry'], 'shape')
           else len(raw['mb_geometry']))
    if key not in cache:
        cache[key] = geom(raw, kind)
    dV, coslat, rm1, rm2, rc, rf = cache[key]
    mb = raw['mb_data']
    rho = np.asarray(mb['dens'], dtype=np.float64)
    if kind == 'cs':
        gg = dhjcs.cs_geometry(raw)
        u, _ = dhjcs.winds(raw, gg)
    else:
        u = np.asarray(mb['velz'], dtype=np.float64)
    cl = coslat[..., None]
    dm = rho*dV
    M = dm.sum()
    # relative-flow term: rho u R dV  with R = <r>_V cos(lat)
    Lrel = (dm*u*(rm1[None, None, None, :]*cl)).sum()
    # frame term: Omega rho R^2 dV
    Lfrm = OMEGA*(dm*(rm2[None, None, None, :]*cl*cl)).sum()
    # mass-weighted mean zonal wind on 3 fixed radial shells
    ish = [int(0.05*len(rc)), int(0.5*len(rc)), int(0.9*len(rc))]
    ub = [float((dm[..., i]*u[..., i]).sum()/dm[..., i].sum()) for i in ish]
    return dict(t=raw['time'], rot=raw['time']/PROT, M=M, Lrel=Lrel, Lfrm=Lfrm,
                ubar=ub, V=dV.sum(),
                Vex=4/3*np.pi*(raw['x1max']**3-raw['x1min']**3), r=rc[ish])


RUNS = [('cs_mhd_prod3', 'cs'), ('sp_mhd_prod3', 'sp'), ('sp_mhd_nopole', 'sp')]

if __name__ == '__main__':
    out = {}
    for name, kind in RUNS:
        files = sorted(glob.glob(B+name+'/bin/*.bin'))
        if len(sys.argv) > 1:
            files = files[::int(sys.argv[1])]
        rows = []
        for f in files:
            try:
                r = one(f, kind)
            except Exception as e:
                print('  skip', f, repr(e))
                continue
            rows.append([r['rot'], r['M'], r['Lrel'], r['Lfrm']] + r['ubar'])
            if len(rows) == 1:
                print('%-14s %s  nfiles=%d  volume sum/exact = %.6f  shells r=%s'
                      % (name, kind, len(files), r['V']/r['Vex'],
                         np.array2string(r['r'], precision=4)))
        a = np.array(rows)
        out[name] = a
        print('%-14s rows=%d  rot %.1f -> %.1f' % (name, len(a), a[0, 0], a[-1, 0]))
    np.savez('/viper/u2/jinma/ATHENAK/athenak/tests_cs_angmom/angmom.npz', **out)
