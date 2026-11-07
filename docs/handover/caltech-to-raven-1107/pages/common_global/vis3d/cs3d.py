#!/usr/bin/env python3
"""3-D cut-away movie for the GLOBAL cubed-sphere runs (gnomonic equiangular, 6 panels N x N, full sphere).

Frame (USER 10-21): the whole star as its CORRUGATED tau_R = 2/3 photosphere (real r_ph per column, exact crossing as
vis3d/render3d.py v7: tau_R integrated inward from the outer boundary on radial faces, kappa_R rho constant per cell, the
run's own Rosseland table at the general-EOS T(rho, e_int) of the run's composition), coloured by
T_eff = (F_r(r_ph)/sigma)^1/4 (F_r = lab m1_f1 linear in r between cell centroids, no smoothing of T_eff; colour at mesh
points = mean of the adjacent cells as render3d 'cap_points'; geometry: 2 light cell->point->cell passes).
view 'single' (USER 10-21): a 90 x 90 deg wedge (lon 0..90 x lat 0..90, apex at the centre) is CUT from the sphere;
its 3 cut faces (quarter discs z=0, y=0, x=0; all on cubed-sphere cell faces) show v_r (coolwarm) out to the local photosphere.
view 'layer' (USER option): nested spheres -- photosphere cut by lat 0..90 x lon +-90, a v_r sphere at shells[0] R_ph,IC
cut by lon +-67.5, a v_r sphere at shells[1] cut by lon +-45; cut faces v_r.  Black background, fixed camera looking into the wedge, time in days.
Colour ranges = 1-99 % of the NEWEST dump (T_eff map; v_r on the cut faces inside r_ph, symmetric).

Incremental: per-dump reductions in <dir>/red/r_NNNNN.npz (key: dump path+size+mtime + REDUCE_VERSION), frames in
<dir>/frames/f_NNNNN.png + .json (key: reduction key + render config incl. colour ranges); only missing/changed frames
are rendered (a new newest dump changes the colour ranges -> all frames re-render, ~s each).

usage: pvpy.sh cs3d.py RUN [--workers W] [--only N,..] [--tag _x] [--no-encode] [--check]
RUN = rsgE_glob | hggl_glob.  Python: /resnick/groups/carnegie_poc/jingze/venv_vtk/pvpy.sh (vtk-osmesa 9.3.1)."""
import argparse
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

VIS = '/resnick/groups/carnegie_poc/jingze/vis3d'
G = '/resnick/groups/carnegie_poc/jingze'
FFMPEG = ('/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/'
          'ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin/ffmpeg')
sys.dont_write_bytecode = True                       # do not drop __pycache__ into the page source dirs
sys.path.insert(0, G + '/rsgE_page/vispy_4de8168d')   # bin reader, format 1.2
sys.path.insert(0, G + '/rsgE_page/src')              # eos_gen (numpy port of the general EOS)
import bin_convert as bc                             # noqa: E402
import eos_gen                                       # noqa: E402

RSUN, DAY, SIGMA_SB = 6.957e10, 86400.0, 5.670374419e-5
SIZE = tuple(int(v) for v in os.environ.get('SIZE3D', '2000x1520').split('x'))   # 10-23 USER: was 1000 x 760 (coarser than the cells)
FS = SIZE[0]/1000.0
MAXFRAMES = 300
MAXBYTES = 15.0e6   # artifact per-file cap 16 MB
REDUCE_VERSION = 'c3'
RENDER_VERSION = 'c3'   # 10-23: 2000 x 1520; T_eff smooth again (USER: per-cell c2 looked blocky)

RUNS = {
    'rsgE_glob': dict(
        title='RSG E (global cubed sphere), Caltech', movie='rsgE_3d',
        hydro=G + '/rsg_E_1020/pageH/bin/rsg.hydro_w.*.bin',   # 10-29 pageH (prodG < 1676 d + prodH)
        table=G + '/rsg_E_1020/files/rsg-g25-prod/rosseland_mesakap_4c0e8b64d26a.txt', ldmin=-21.0,
        comp=dict(xh=0.653848, yhe=0.326395, a_metal=16.0),
        cmap_teff=('hot', 0.22, 1.0),                 # USER 10-21: heat-like for the RSG; dark end cut (black bg)
        rph0=5.80019e13, shells=(0.7, 0.4), relief=1,   # relief: radial exaggeration of r_ph - median (USER 10-21)          # layered view: inner spheres at 0.7 / 0.4 R_ph,IC (USER example)
        publish=G + '/rsgE_page/out/publish/movies'),
    'hggl_glob': dict(
        title='He giant global (cubed sphere), Caltech', movie='hggl_3d',
        hydro=G + '/hgglob_1021/prodG/bin/hegiant.hydro_w.*.bin',   # 10-22 fresh N^2-seed run
        table=G + '/hgglob_1021/files/hg145-prod-1018/rosseland_mesakap_afe2dd8da2bd.txt', ldmin=-21.0,
        comp=dict(xh=0.0009496, yhe=0.97752, a_metal=1.8047),
        cmap_teff=('Blues', 0.0, 1.0),
        rph0=7.94653e12, shells=(0.75, 0.35), relief=1,        # IC: F_MLT/F peaks at ~0.33 and ~0.95 R_ph, CZ 0.07-0.99 R_ph
        publish=G + '/hgglob_page/out/publish/movies'),
}
VIEWS = {   # single: one 90 (lon) x 90 (lat) wedge (octant) cut; layer: nested spheres cut by 180/135/90 (lon) x 90 (lat) wedges
    'single': dict(az=45.0, el=30.0, suffix=''),
    'layer': dict(az=25.0, el=33.0, suffix='_layers', zoom=0.86, wcenter=(0.0, 0.03)),
}
LAYER_L = (90.0, 67.5, 45.0)          # half-widths in longitude of the 3 nested wedges (lat 0..90, centred on lon 0)
for _c in RUNS.values():
    _c.setdefault('view', 'single')
    _c.setdefault('zoom', 0.93)
    _c.setdefault('wcenter', (0.16, 0.03))

# cubed-sphere panel frames (AthenaK cubed_sphere.hpp; = rsgE_page/src/csplot.py): direction = (a tan xi + b tan eta + n)/delta
FR = {0: ((0, 1, 0), (0, 0, 1), (1, 0, 0)), 1: ((-1, 0, 0), (0, 0, 1), (0, 1, 0)), 2: ((0, -1, 0), (0, 0, 1), (-1, 0, 0)),
      3: ((0, 1, 0), (-1, 0, 0), (0, 0, 1)), 4: ((1, 0, 0), (0, 0, 1), (0, -1, 0)), 5: ((0, 1, 0), (1, 0, 0), (0, 0, -1))}
FA = np.array([FR[p][0] for p in range(6)], float)
FB = np.array([FR[p][1] for p in range(6)], float)
FN = np.array([FR[p][2] for p in range(6)], float)
S45 = np.sqrt(0.5)


# ---------------------------------------------------------------- data
def read_ross(fn):
    lT0 = None
    with open(fn) as fh:
        for ln in fh:
            if not ln.startswith('#'):
                break
            w = ln[1:].split()
            if len(w) == 6:
                try:
                    nT, nD, lT0, dlT, lD0, dlD = [float(v) for v in w]
                except ValueError:
                    pass
    nT, nD = int(nT), int(nD)
    v = np.loadtxt(fn, comments='#').reshape(nT, nD)
    return lT0, dlT, nT, lD0, dlD, nD, v


def extend_ldmin(tab, ldmin):
    lT0, dlT, nT, lD0, dlD, nD, v = tab
    if ldmin is None or ldmin >= lD0 - 1e-9*dlD:
        return tab
    nadd = int(np.ceil((lD0 - ldmin)/dlD - 1e-9))
    sl = np.clip((v[:, 1] - v[:, 0])/dlD, 0.0, 1.0)
    ext = v[:, :1] - sl[:, None]*(nadd - np.arange(nadd))[None, :]*dlD
    return lT0, dlT, nT, lD0 - nadd*dlD, dlD, nD + nadd, np.concatenate([ext, v], axis=1)


def kappa(tab, rho, T):
    """bilinear in (log T, log rho), clamped to the table box (as the code / render3d.py)"""
    lT0, dlT, nT, lD0, dlD, nD, v = tab
    x = np.clip((np.log10(np.maximum(T, 1.0)) - lT0)/dlT, 0, nT - 1.000001)
    y = np.clip((np.log10(np.maximum(rho, 1e-30)) - lD0)/dlD, 0, nD - 1.000001)
    i, j = x.astype(np.int32), y.astype(np.int32)
    fx, fy = x - i, y - j
    lk = ((1 - fx)*((1 - fy)*v[i, j] + fy*v[i, j + 1]) + fx*((1 - fy)*v[i + 1, j] + fy*v[i + 1, j + 1]))
    return 10.0**lk


def assemble(d, var):
    """-> (6, N, N, nr) [panel, k(eta), j(xi), i]; MeshBlocks panel-major (as rsgE_page/src/reduce.py)"""
    nmb = len(d['mb_logical'])
    per = nmb // 6
    n3, n2, n1 = d['mb_data'][var][0].shape
    N = n2*int(round(np.sqrt(per)))
    P = np.empty((6, N, N, n1), np.float32)
    for m in range(nmb):
        p = m // per
        lx1, lx2, lx3, lev = d['mb_logical'][m]
        P[p, lx3*n3:(lx3 + 1)*n3, lx2*n2:(lx2 + 1)*n2, :] = d['mb_data'][var][m]
    return P


def lookup(q, N):
    """unit vectors q (..., 3) -> (panel, k, j) of the cell containing them"""
    p = np.argmax(q @ FN.T, axis=-1)
    qa = np.einsum('...c,...c->...', q, FA[p])
    qb = np.einsum('...c,...c->...', q, FB[p])
    qn = np.einsum('...c,...c->...', q, FN[p])
    j = np.clip(((np.arctan(qa/qn)/(np.pi/4) + 1)*0.5*N).astype(int), 0, N - 1)
    k = np.clip(((np.arctan(qb/qn)/(np.pi/4) + 1)*0.5*N).astype(int), 0, N - 1)
    return p, k, j


def centroid(rl, rr):
    q = rl/rr
    return 0.25*(q*q + 1)/((q*q + q + 1)/3)*(rr + rl)


def face_defs(N):
    """cut faces of the wedge lon 0..90 x lat 0..45: (name, edge dirs (M+1,3), centre dirs (M,3), normal(s))"""
    out = []
    w = np.pi/(2*N)                                     # native cell width along the panel centre lines

    def mk(name, s0, s1, M, f, nrm):
        se = np.linspace(s0, s1, M + 1)
        sc = 0.5*(se[1:] + se[:-1])
        out.append((name, f(se), f(sc), nrm(sc)))
    z = lambda s: np.zeros_like(s)                      # noqa: E731
    mk('eq', 0, np.pi/2, N, lambda s: np.stack([np.cos(s), np.sin(s), z(s)], -1),
       lambda s: np.tile([0.0, 0.0, 1.0], (len(s), 1)))
    mk('m0', 0, np.pi/2, N, lambda s: np.stack([np.cos(s), z(s), np.sin(s)], -1),
       lambda s: np.tile([0.0, 1.0, 0.0], (len(s), 1)))
    mk('m90', 0, np.pi/2, N, lambda s: np.stack([z(s), np.cos(s), np.sin(s)], -1),
       lambda s: np.tile([1.0, 0.0, 0.0], (len(s), 1)))
    # layered view: full equatorial ring (lon -180..180) + meridional half-planes (lat 0..90) at lon = +-L
    mk('eqF', -np.pi, np.pi, 4*N, lambda s: np.stack([np.cos(s), np.sin(s), z(s)], -1),
       lambda s: np.tile([0.0, 0.0, 1.0], (len(s), 1)))
    for L in LAYER_L:
        for sg in (1, -1):
            lo = np.radians(sg*L)
            mk('mer%+05.1f' % (sg*L), 0, np.pi/2, N,
               lambda s, lo=lo: np.stack([np.cos(s)*np.cos(lo), np.cos(s)*np.sin(lo), np.sin(s)], -1),
               lambda s, lo=lo: np.tile([-np.sin(lo), np.cos(lo), 0.0], (len(s), 1)))
    return out, w


def face_cells(dc, nrm, N, w):
    """the two cells either side of each face sample (offset +-0.3 cell along the face normal)"""
    res = []
    for sg in (1, -1):
        q = dc + sg*0.3*w*nrm
        q /= np.linalg.norm(q, axis=-1, keepdims=True)
        res.append(lookup(q, N))
    return res


_EOS = {}


def eos_table(run, comp):
    if run not in _EOS:
        cf = os.path.join(VIS, run, 'eos_fe.npy')
        tb = eos_gen.Table.__new__(eos_gen.Table)
        tbl_args = dict(logd=(-21.0, -1.0), logt=(3.3, 6.8), dlog=0.02)
        if os.path.exists(cf):
            tb.fe = np.load(cf)
            tb.ny, tb.nx = tb.fe.shape
            tb.x0, tb.dx = -21.0, 20.0/(tb.nx - 1)
            tb.y0, tb.dy = 3.3, 3.5/(tb.ny - 1)
        else:
            tb = eos_gen.Table(**tbl_args, **comp)
            np.save(cf + '.%d.npy' % os.getpid(), tb.fe)
            os.replace(cf + '.%d.npy' % os.getpid(), cf)
        _EOS[run] = tb
    return _EOS[run]


def reduce_dump(run, fh, cfg):
    """compact per-dump data: per column (6,N,N) r_ph, idx, F_r(r_ph), T_ph; v_r rings on the 4 cut faces"""
    t0 = time.time()
    tab = extend_ldmin(read_ross(cfg['table']), cfg.get('ldmin'))
    tb = eos_table(run, cfg['comp'])
    a = bc.read_binary(fh)
    t = float(a['time'])
    rf = np.asarray(a['mb_x1f'][0], float)
    rho = assemble(a, 'dens')
    vr = assemble(a, 'velx')
    ei = assemble(a, 'eint')
    del a
    b = bc.read_binary(fh.replace('hydro_w', 'm1'))
    F1 = assemble(b, 'm1_f1')
    del b
    tread = time.time() - t0
    N, nr = rho.shape[1], rho.shape[3]
    dr = np.diff(rf)
    rc = centroid(rf[:-1], rf[1:])
    rph = np.empty((6, N, N))
    idx = np.empty((6, N, N), np.int32)
    Fph = np.empty((6, N, N))
    Tph = np.empty((6, N, N))
    ish = [int(np.argmin(np.abs(rc - x*cfg['rph0']))) for x in cfg['shells']]
    Tsh = np.empty((len(ish), 6, N, N))
    for p in range(6):
        r_ = rho[p].astype(np.float64)
        T = tb.temperature(r_, ei[p])[0].astype(np.float64)
        Tsh[:, p] = np.moveaxis(T[..., ish], -1, 0)
        kr = kappa(tab, r_, T)*r_
        tauf = np.zeros((N, N, nr + 1))
        tauf[..., :nr] = np.flip(np.cumsum(np.flip(kr*dr, -1), -1), -1)
        t23 = 2.0/3.0
        ge = tauf[..., :nr] >= t23
        i = nr - 1 - np.argmax(ge[..., ::-1], axis=-1)
        ta = np.take_along_axis(tauf, (i + 1)[..., None], -1)[..., 0]
        k = np.take_along_axis(kr, i[..., None], -1)[..., 0]
        rp = np.where(ge[..., 0], rf[i + 1] - (t23 - ta)/np.maximum(k, 1e-300), rf[0])
        rph[p], idx[p] = rp, np.where(ge[..., 0], i, 0)
        ic = np.clip(np.searchsorted(rc, rp) - 1, 0, nr - 2)
        ww = np.clip((rp - rc[ic])/(rc[ic + 1] - rc[ic]), 0.0, 1.0)

        def at(Q):
            q0 = np.take_along_axis(Q, ic[..., None], -1)[..., 0]
            q1 = np.take_along_axis(Q, (ic + 1)[..., None], -1)[..., 0]
            return (1 - ww)*q0 + ww*q1
        Fph[p] = at(F1[p].astype(np.float64))
        Tph[p] = 10**at(np.log10(T))
        del T, kr, tauf, ge
    out = dict(t=t, rf=rf, rph=rph, idx=idx, Fph=Fph, Tph=Tph.astype(np.float32))
    cc = np.tan((-1 + (np.arange(N) + 0.5)*2.0/N)*np.pi/4)          # cell solid-angle weights (csplot.shell_mean)
    X, Y = np.meshgrid(cc, cc, indexing='xy')
    wom = (1 + X**2)*(1 + Y**2)/(1 + X**2 + Y**2)**1.5
    for m, x in enumerate(cfg['shells']):               # native shells (layered view), nearest radial cell
        i = ish[m]
        out['sh_%03d' % round(100*x)] = (vr[..., i]/1e5).astype(np.float32)
        out['shr_%03d' % round(100*x)] = rc[i]
        Tm = float((Tsh[m]*wom[None]).sum()/(6*wom.sum()))
        out['shT_%03d' % round(100*x)] = (100*(Tsh[m]/Tm - 1)).astype(np.float32)   # dT_gas/<T_gas>_shell [%]
    faces, w = face_defs(N)
    for name, de, dc, nrm in faces:
        (p1, k1, j1), (p2, k2, j2) = face_cells(dc, nrm, N, w)
        out['vr_' + name] = (0.5*(vr[p1, k1, j1] + vr[p2, k2, j2])/1e5).astype(np.float32)   # (M, nr) km/s
    print('reduced %s t=%.3f d  R_ph med %.1f Rsun  F<=0 %.2e  read %.0f s total %.0f s' % (
        os.path.basename(fh), t/DAY, np.median(rph)/RSUN, np.mean(Fph <= 0), tread, time.time() - t0), flush=True)
    return out


def load_reduced(run, n, fh, key, cfg):
    d = os.path.join(VIS, run, 'red')
    os.makedirs(d, exist_ok=True)
    cf = os.path.join(d, 'r_%05d.npz' % n)
    rk = key + '|' + REDUCE_VERSION + '|' + json.dumps([cfg['table'], cfg['comp'], cfg.get('ldmin')])
    if os.path.exists(cf):
        try:
            z = np.load(cf)
            if str(z['key']) == rk:
                return {k: (z[k][()] if z[k].ndim == 0 else z[k]) for k in z.files if k != 'key'}
        except Exception:
            pass
    D = reduce_dump(run, fh, cfg)
    np.savez(cf + '.tmp.npz', key=rk, **D)
    os.replace(cf + '.tmp.npz', cf)
    return D


def teff(F):
    fl = 1e-6*max(float(np.median(np.abs(F))), 1e-30)
    return (np.maximum(F, fl)/SIGMA_SB)**0.25


# ---------------------------------------------------------------- geometry (global, seamless)
_MESH = {}


def cs_mesh(N):
    """unique mesh points of the 6 panels (seams shared): unit vectors P (np,3), quads Q (6N^2,4), cell centres C"""
    if N in _MESH:
        return _MESH[N]
    e = np.tan((-1 + np.arange(N + 1)*2.0/N)*np.pi/4)
    c = np.tan((-1 + (np.arange(N) + 0.5)*2.0/N)*np.pi/4)
    pts, cen = [], []
    for p in range(6):
        X, Y = np.meshgrid(e, e, indexing='xy')                       # [k, j]
        v = FA[p]*X[..., None] + FB[p]*Y[..., None] + FN[p]
        pts.append((v/np.linalg.norm(v, axis=-1, keepdims=True)).reshape(-1, 3))
        X, Y = np.meshgrid(c, c, indexing='xy')
        v = FA[p]*X[..., None] + FB[p]*Y[..., None] + FN[p]
        cen.append(v/np.linalg.norm(v, axis=-1, keepdims=True))
    allp = np.concatenate(pts)
    keyq = np.round(allp*2e6).astype(np.int64)
    _, first, inv = np.unique(keyq, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    P = allp[first]
    n1 = N + 1
    k, j = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    Q = []
    for p in range(6):
        base = p*n1*n1
        loc = np.stack([k*n1 + j, k*n1 + j + 1, (k + 1)*n1 + j + 1, (k + 1)*n1 + j], -1) + base
        Q.append(inv[loc].reshape(-1, 4))
    Q = np.concatenate(Q)
    C = np.concatenate(cen).reshape(6, N, N, 3)
    _MESH[N] = (P, Q, C)
    return _MESH[N]


def c2p(vals, Q, npt):
    s = np.zeros(npt)
    cnt = np.zeros(npt)
    for m in range(4):
        np.add.at(s, Q[:, m], vals)
        np.add.at(cnt, Q[:, m], 1.0)
    return s/cnt


def lut_from(cmapspec, lo, hi, n=256):
    import vtk
    import matplotlib
    name, a, b = cmapspec
    cm = matplotlib.colormaps[name]
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(n)
    lut.SetRange(lo, hi)
    for i, x in enumerate(np.linspace(a, b, n)):
        r, g, bb, _ = cm(x)
        lut.SetTableValue(i, r, g, bb, 1.0)
    lut.Build()
    return lut


def np2vtk_pts(a):
    import vtk
    from vtkmodules.util.numpy_support import numpy_to_vtk
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(a, dtype=np.float64), deep=True))
    return pts


def quads_poly(P, Q, cell=None, point=None):
    import vtk
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    pd = vtk.vtkPolyData()
    pd.SetPoints(np2vtk_pts(P))
    ca = vtk.vtkCellArray()
    nq = len(Q)
    off = numpy_to_vtkIdTypeArray(np.arange(0, 4*nq + 1, 4, dtype=np.int64), deep=True)
    con = numpy_to_vtkIdTypeArray(np.ascontiguousarray(Q.ravel(), dtype=np.int64), deep=True)
    ca.SetData(off, con)
    pd.SetPolys(ca)
    if cell is not None:
        arr = numpy_to_vtk(np.ascontiguousarray(cell, dtype=np.float64), deep=True)
        arr.SetName('c')
        pd.GetCellData().SetScalars(arr)
    if point is not None:
        arr = numpy_to_vtk(np.ascontiguousarray(point, dtype=np.float64), deep=True)
        arr.SetName('c')
        pd.GetPointData().SetScalars(arr)
    return pd


def wedge_f(U):
    """single view: < 0 inside the removed 90 x 90 deg wedge (lon 0..90, lat 0..90 = octant x, y, z > 0)"""
    return np.max(np.stack([-U[:, 0], -U[:, 1], -U[:, 2]]), axis=0)


def layer_f(U, L):
    """layered view: < 0 inside the wedge lat 0..90 x |lon| < L (deg)"""
    rxy = np.sqrt(U[:, 0]**2 + U[:, 1]**2)
    return np.maximum(-U[:, 2], np.cos(np.radians(L))*rxy - U[:, 0])


def clip(pd, fv, keep_neg):
    import vtk
    from vtkmodules.util.numpy_support import numpy_to_vtk
    f = numpy_to_vtk(np.ascontiguousarray(fv, dtype=np.float64), deep=True)
    f.SetName('f')
    pd.GetPointData().AddArray(f)
    pd.GetPointData().SetActiveScalars('f')
    cl = vtk.vtkClipPolyData()
    cl.SetInputData(pd)
    cl.SetValue(0.0)
    cl.SetInsideOut(1 if keep_neg else 0)
    cl.Update()
    out = cl.GetOutput()
    out.GetPointData().RemoveArray('f')
    if out.GetPointData().GetArray('c') is not None:
        out.GetPointData().SetActiveScalars('c')
    return out


def unit_pts(pd):
    from vtkmodules.util.numpy_support import vtk_to_numpy
    P = vtk_to_numpy(pd.GetPoints().GetData()).astype(np.float64)
    return P/np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-30)


def smooth_r(D, relief=1.0):
    """per-column r_ph [Rsun] (exact tau_R = 2/3 crossing per radial column), 2 light cell->point->cell passes, then
    the deviation from the median exaggerated x relief (USER 10-21: the real relief is ~1e-5 early on)"""
    N = D['rph'].shape[1]
    P, Q, C = cs_mesh(N)
    r = (D['rph']/RSUN).ravel()
    med = float(np.median(r))
    for _ in range(2):                                   # light geometry smoothing (render3d rph_smooth 2 analogue)
        rp = c2p(r, Q, len(P))
        r = rp[Q].mean(1)
    r = np.clip(med + relief*(r - med), 0.5*med, 2.0*med)   # guard: exaggerated relief kept within 0.5-2 x median
    return c2p(r, Q, len(P)), r.reshape(6, N, N)


def surface(D, view, relief):
    """photosphere polydata (smoothed r_ph geometry, T_eff [kK] at points), clipped exactly by the view's wedge"""
    N = D['rph'].shape[1]
    P, Q, C = cs_mesh(N)
    rp, rsm = smooth_r(D, relief)
    te = c2p((teff(D['Fph'])/1e3).ravel(), Q, len(P))   # T_eff at mesh points (USER 10-23: keep smooth, not per cell)
    pd = quads_poly(P*rp[:, None], Q, point=te)
    fv = wedge_f(P) if view == 'single' else layer_f(P, LAYER_L[0])
    return clip(pd, fv, False), rsm, rp


def shells(D, cfg):
    """layered view: native v_r shells at the 2 inner radii, each shown only between its own and the next wedge"""
    N = D['rph'].shape[1]
    P, Q, C = cs_mesh(N)
    res = []
    for m, x in enumerate(cfg['shells']):
        tg = '%03d' % round(100*x)
        r = float(D['shr_' + tg])/RSUN
        pd = quads_poly(P*r, Q, cell=D['shT_' + tg].ravel())
        pd = clip(pd, layer_f(P, LAYER_L[m]), True)           # inside the outer wedge
        pd = clip(pd, layer_f(unit_pts(pd), LAYER_L[m + 1]), False)   # outside the inner one
        res.append((pd, r))
    return res


def face_poly(de, V, ra, rb, rf, rbe=None):
    """quads (r, angle) on a cut face: radial cells clipped to [ra, rb] per angle sample; cell scalars V (M, nr).
    rbe (M+1,): outer radius at the sample EDGES (the photosphere mesh points on the cut line), so the face's outer
    edge coincides with the clipped surface edge"""
    M = len(ra)
    top = rb if rbe is None else np.maximum(rbe[:-1], rbe[1:])
    nmax = min(int(np.searchsorted(rf, top.max())) + 1, len(rf) - 1)
    Rl = np.maximum(rf[:nmax][None, :], ra[:, None])
    Rr = np.minimum(rf[1:nmax + 1][None, :], top[:, None])
    mi, ii = np.nonzero(Rr > Rl)
    r0 = Rl[mi, ii]
    if rbe is None:
        r1a = r1b = Rr[mi, ii]
    else:
        r1a = np.maximum(np.minimum(rf[ii + 1], rbe[mi]), r0)
        r1b = np.maximum(np.minimum(rf[ii + 1], rbe[mi + 1]), r0)
    a0, a1 = de[mi], de[mi + 1]
    pts = np.stack([a0*r0[:, None], a1*r0[:, None], a1*r1b[:, None], a0*r1a[:, None]], 1).reshape(-1, 3)
    return quads_poly(pts, np.arange(len(pts)).reshape(-1, 4), cell=V[mi, ii]), V[mi, ii]


_TREE = {}


def _tree(N):
    if N not in _TREE:
        from scipy.spatial import cKDTree
        _TREE[N] = cKDTree(cs_mesh(N)[0])
    return _TREE[N]


def faces(D, rsm, cfg, rpt=None):
    """cut-face polydata (cell scalars v_r km/s) of the view, out to the local (smoothed) photosphere"""
    N = D['rph'].shape[1]
    rf = D['rf']/RSUN
    fdefs, w = face_defs(N)
    view = cfg['view']
    if view == 'layer':
        rs = [float(D['shr_%03d' % round(100*x)])/RSUN for x in cfg['shells']]   # r2 > r3
    res = []
    for name, de, dc, nrm in fdefs:
        (p1, k1, j1), (p2, k2, j2) = face_cells(dc, nrm, N, w)
        redge = 0.5*(rsm[p1, k1, j1] + rsm[p2, k2, j2])
        rbe = None
        if rpt is not None:                              # photosphere radius at the nearest mesh point of each edge
            rbe = rpt[_tree(N).query(de)[1]]
        M = len(dc)
        if view == 'single':
            if name not in ('eq', 'm0', 'm90'):
                continue
            ra = np.zeros(M)
            rb = redge
        else:
            if name == 'eqF':
                # floor: every radius inside the outer wedge is removed above z = 0 (the 3 wedges nest), so the floor
                # runs from the centre to the photosphere for all |lon| < 90
                lon = np.degrees(np.abs(np.arctan2(dc[:, 1], dc[:, 0])))
                ra = np.zeros(M)
                rb = np.where(lon < LAYER_L[0], redge, 0.0)
            elif name.startswith('mer'):
                L = abs(float(name[3:]))
                m = LAYER_L.index(L)
                ra = np.full(M, (rs[0] if m == 0 else rs[1] if m == 1 else 0.0))
                rb = redge if m == 0 else np.full(M, rs[0] if m == 1 else rs[1])
            else:
                continue
        outer = (view == 'single') or name == 'eqF' or name.startswith('mer') and abs(float(name[3:])) == LAYER_L[0]
        if outer and rbe is not None:
            if name == 'eqF':                            # only inside |lon| < 90; elsewhere nothing
                lon_e = np.degrees(np.abs(np.arctan2(de[:, 1], de[:, 0])))
                rbe = np.where(lon_e <= LAYER_L[0] + 1e-6, rbe, 0.0)
            res.append(face_poly(de, D['vr_' + name], ra, rb, rf, rbe))
        else:
            res.append(face_poly(de, D['vr_' + name], ra, rb, rf))
    return res


def render(D, cfg, clim, out):
    import vtk
    rw = vtk.vtkRenderWindow()
    rw.SetOffScreenRendering(1)
    rw.SetSize(*SIZE)
    ren = vtk.vtkRenderer()
    ren.SetBackground(0, 0, 0)
    rw.AddRenderer(ren)
    ren.AutomaticLightCreationOff()
    surf, rsm, rpt = surface(D, cfg['view'], cfg['relief'])
    vlut = lut_from(('coolwarm', 0.0, 1.0), -clim['vr'], clim['vr'])
    tlut = lut_from(cfg['cmap_teff'], *clim['Teff'])
    fl = [(pd, vlut) for pd, _ in faces(D, rsm, cfg, rpt)]
    if cfg['view'] == 'layer':
        dlut = lut_from(('PuOr_r', 0.0, 1.0), -clim['dT'], clim['dT'])
        fl += [(pd, dlut) for pd, _ in shells(D, cfg)]
    for pd, lu in fl:
        m = vtk.vtkPolyDataMapper()
        m.SetInputData(pd)
        m.SetLookupTable(lu)
        m.SetScalarModeToUseCellData()
        m.UseLookupTableScalarRangeOn()
        a = vtk.vtkActor()
        a.SetMapper(m)
        pr = a.GetProperty()
        pr.SetAmbient(0.75); pr.SetDiffuse(0.3); pr.SetSpecular(0.0)
        if lu is not vlut:
            pr.BackfaceCullingOn()                       # shells: outer side only
        ren.AddActor(a)
    nf = vtk.vtkPolyDataNormals()
    nf.SetInputData(surf)
    nf.SplittingOff()
    nf.ConsistencyOn()
    nf.Update()
    m = vtk.vtkPolyDataMapper()
    m.SetInputData(nf.GetOutput())
    m.SetLookupTable(tlut)
    m.SetScalarModeToUsePointData()
    m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor()
    a.SetMapper(m)
    pr = a.GetProperty()
    pr.BackfaceCullingOn()                               # never show the inside of the photosphere through the cut
    pr.SetInterpolationToGouraud()
    pr.SetAmbient(0.3); pr.SetDiffuse(0.7); pr.SetSpecular(0.15); pr.SetSpecularPower(20)
    ren.AddActor(a)
    sp = vtk.vtkSphereSource()                           # inner boundary (grey)
    sp.SetRadius(D['rf'][0]/RSUN)
    sp.SetThetaResolution(48); sp.SetPhiResolution(24)
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(sp.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.33, 0.33, 0.33); a.GetProperty().SetAmbient(0.3)
    ren.AddActor(a)
    # camera: fixed, looking into the wedge
    az, el = np.radians(VIEWS[cfg['view']]['az']), np.radians(VIEWS[cfg['view']]['el'])
    eye = np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])
    R = clim['rview']
    cam = ren.GetActiveCamera()
    cam.SetViewAngle(30)
    vw = VIEWS[cfg['view']]
    dist = R/np.sin(np.radians(15))*1.02/vw.get('zoom', cfg['zoom'])
    cam.SetFocalPoint(0, 0, 0)
    cam.SetPosition(*(eye*dist))
    cam.SetViewUp(0, 0, 1)
    cam.SetWindowCenter(*vw.get('wcenter', cfg['wcenter']))
    ren.ResetCameraClippingRange()
    hl = vtk.vtkLight(); hl.SetLightTypeToHeadlight(); hl.SetIntensity(0.55); ren.AddLight(hl)
    sl = vtk.vtkLight(); sl.SetLightTypeToSceneLight(); sl.SetPositional(False)
    ld = eye + np.array([0.0, 0.0, 0.8]) + 0.6*np.array([-np.sin(az), np.cos(az), 0.0])
    sl.SetPosition(*(ld*dist)); sl.SetFocalPoint(0, 0, 0); sl.SetIntensity(0.55); ren.AddLight(sl)

    def sbar(lut, title, vertical, pos, size):
        sb = vtk.vtkScalarBarActor()
        sb.SetLookupTable(lut)
        sb.SetTitle(title)
        sb.SetNumberOfLabels(3)
        sb.SetLabelFormat('%.3g')
        sb.UnconstrainedFontSizeOn()
        for tp, fs in ((sb.GetTitleTextProperty(), round(24*FS)), (sb.GetLabelTextProperty(), round(22*FS))):
            tp.SetFontSize(fs); tp.SetColor(1, 1, 1); tp.BoldOff(); tp.ItalicOff(); tp.ShadowOff()
        sb.SetOrientationToVertical() if vertical else sb.SetOrientationToHorizontal()
        sb.SetPosition(*pos)
        sb.SetPosition2(*size)
        sb.SetBarRatio(0.35 if vertical else 0.30)
        ren.AddActor2D(sb)
    sbar(vlut, 'radial velocity (km/s)', False, (0.22, 0.01), (0.50, 0.11))
    sbar(tlut, 'T_eff (kK)', True, (0.84, 0.20), (0.11, 0.66))
    if cfg['view'] == 'layer':
        sbar(dlut, 'dT/<T> (%)', True, (0.02, 0.20), (0.11, 0.66))
    tx = vtk.vtkTextActor()
    tx.SetInput('%s\nt = %.2f d   median R_ph = %.0f R_sun   photosphere at true scale' % (
        cfg['title'], D['t']/DAY, np.median(D['rph'])/RSUN))
    tp = tx.GetTextProperty()
    tp.SetFontSize(round(28*FS)); tp.SetColor(1, 1, 1); tp.SetJustificationToCentered(); tp.SetVerticalJustificationToTop()
    tx.SetPosition(SIZE[0]/2, SIZE[1] - 8*FS)
    ren.AddActor2D(tx)
    rw.Render()
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(rw)
    w2i.ReadFrontBufferOff()
    w2i.Update()
    wr = vtk.vtkPNGWriter()
    wr.SetFileName(out)
    wr.SetInputConnection(w2i.GetOutputPort())
    wr.Write()
    rw.Finalize()


# ---------------------------------------------------------------- driver
def dump_list(cfg):
    out = []
    now = time.time()
    for fh in sorted(glob.glob(cfg['hydro'])):
        fm = fh.replace('hydro_w', 'm1')
        if not os.path.exists(fm):
            continue
        try:
            st = [os.stat(f) for f in (fh, fm)]
        except OSError:
            continue
        if now - max(s.st_mtime for s in st) < 180:     # maybe still being written
            continue
        n = int(fh.split('.')[-2])
        key = '|'.join('%s:%d:%d' % (os.path.realpath(f), s.st_size, int(s.st_mtime)) for f, s in zip((fh, fm), st))
        out.append((n, fh, key))
    if len(out) > MAXFRAMES:
        s = int(np.ceil(len(out)/MAXFRAMES))
        out = out[::s] if (len(out) - 1) % s == 0 else out[::s] + [out[-1]]
    return out


def red_one(args):
    run, n, fh, key = args
    D = load_reduced(run, n, fh, key, RUNS[run])
    med = float(np.median(D['rph']))
    lo, hi = np.percentile(D['rph'], [0.5, 99.5])
    return n, D['t'], med, float(lo), float(hi)


def auto_relief(st, cap):
    """one relief for the whole movie: the largest 0.5-99.5 % r_ph excursion over all dumps -> 4 % of the median
    (USER 10-22: x1000 became spiky once the real relief grew); 1-2-5 steps keep the frame cache stable"""
    x = max(max(hi - m, m - lo)/m for _, _, m, lo, hi in st)
    want = min(cap, 0.04/max(x, 1e-12))
    e = 10.0**np.floor(np.log10(want))
    return float(max(s*e for s in (1, 2, 5) if s*e <= want))


def sig2(x):
    if x == 0:
        return 0.0
    e = int(np.floor(np.log10(abs(x)))) - 2              # 3 significant digits
    return float(round(x, -e))


def colour_limits(run, cfg, n, fh, key):
    """1-99 % of the newest dump: T_eff map; |v_r| on the cut faces inside r_ph (symmetric)"""
    D = load_reduced(run, n, fh, key, cfg)
    te = teff(D['Fph'])/1e3
    lo, hi = np.percentile(te, [1, 99])
    rpt, rsm = smooth_r(D, cfg['relief'])
    vs = np.concatenate([v for _, v in faces(D, rsm, cfg, rpt)])
    p1, p99 = np.percentile(vs, [1, 99])
    if hi - lo < 0.02*hi:                                # degenerate (IC dump): +-5 % around the value
        mid = 0.5*(lo + hi)
        lo, hi = 0.95*mid, 1.05*mid
    out = dict(Teff=[sig2(lo), sig2(hi)], vr=max(sig2(max(abs(p1), abs(p99))), 0.1))
    if cfg['view'] == 'layer':
        d = np.concatenate([D['shT_%03d' % round(100*x)].ravel() for x in cfg['shells']])
        q1, q99 = np.percentile(d, [1, 99])
        out['dT'] = max(sig2(max(abs(q1), abs(q99))), 0.01)
    return out


def sticky_clim(path, new):
    """USER 10-22 "page refreshes really fast": keep the previous colour limits unless the newest dump leaves them
    by > 25 % (range or |v| scale), so a new dump renders only its own frame instead of re-rendering the movie"""
    try:
        old = json.load(open(path))
    except Exception:
        old = None
    def near(a, b):
        if isinstance(a, list):
            span = max(b[1] - b[0], 1e-30)
            return abs(a[0] - b[0]) <= 0.25*span and abs(a[1] - b[1]) <= 0.25*span
        return 0.8*b <= a <= 1.25*b
    def ok(k):                                           # the view radius may shrink a little, never grow (clipping)
        return 0.9*old[k] <= new[k] <= old[k] if k == 'rview' else near(new[k], old[k])
    if old and set(old) == set(new) and old['relief'] == new['relief'] \
            and all(ok(k) for k in new if k != 'relief'):
        return old
    json.dump(new, open(path, 'w'))
    return new


def work(args):
    run, n, fh, key, fkey, clim, tag, view = args
    cfg = dict(RUNS[run], view=view, relief=clim['relief'])
    fdir = os.path.join(VIS, run, 'test' if tag else 'frames' + VIEWS[cfg['view']]['suffix'])
    os.makedirs(fdir, exist_ok=True)
    t0 = time.time()
    D = load_reduced(run, n, fh, key, cfg)
    t1 = time.time()
    out = os.path.join(fdir, 'f_%05d%s%s.png' % (n, VIEWS[cfg['view']]['suffix'] if tag else '', tag))
    render(D, cfg, clim, out + '.tmp.png')
    os.replace(out + '.tmp.png', out)
    te = teff(D['Fph'])/1e3
    rec = dict(n=n, key=fkey, t=D['t'], rph_med=float(np.median(D['rph']))/RSUN, rph_max=float(D['rph'].max())/RSUN,
               teff_p=[float(x) for x in np.percentile(te, [1, 50, 99])], frac_Fle0=float(np.mean(D['Fph'] <= 0)),
               t_load=t1 - t0, t_render=time.time() - t1)
    if not tag:
        json.dump(rec, open(os.path.join(fdir, 'f_%05d.json' % n), 'w'))
    print('frame %05d t=%.2f d  R_ph med %.1f max %.1f Rsun  T_eff kK p1/50/99 %s  F<=0 %.2e  load %.1f s render %.1f s' % (
        n, D['t']/DAY, rec['rph_med'], rec['rph_max'], np.round(rec['teff_p'], 2), rec['frac_Fle0'], rec['t_load'],
        rec['t_render']), flush=True)
    return rec


def encode(run, cfg, recs, clim):
    sfx = VIEWS[cfg['view']]['suffix']
    fdir = os.path.join(VIS, run, 'frames' + sfx)
    od = os.path.join(VIS, run, 'out')
    os.makedirs(od, exist_ok=True)
    recs = sorted(recs, key=lambda q: q['t'])
    if len(recs) > 2:                                    # drop duplicate times (restart dumps written twice)
        dtm = float(np.median(np.diff([q['t'] for q in recs])))
        keep = [recs[0]]
        for q in recs[1:]:
            if q['t'] - keep[-1]['t'] > 0.1*dtm:
                keep.append(q)
        recs = keep
    tmp = os.path.join(VIS, run, 'seq' + sfx)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    for i, q in enumerate(recs):
        os.symlink(os.path.join(fdir, 'f_%05d.png' % q['n']), os.path.join(tmp, 's_%05d.png' % i))
    nf = len(recs)
    fps = int(round(1.3*max(8, min(24, round(nf/15.0)))))   # as render3d.py (x1.3, user 10-15)
    mp4 = os.path.join(od, cfg['movie'] + sfx + '.mp4')
    for crf in (22, 26, 30, 34, 38):
        subprocess.run([FFMPEG, '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', tmp + '/s_%05d.png',
                        '-c:v', 'libx264', '-preset', 'slow', '-crf', str(crf), '-pix_fmt', 'yuv420p',
                        '-profile:v', 'high', '-level', '4.0', '-movflags', '+faststart', mp4 + '.tmp.mp4'], check=True)
        if os.path.getsize(mp4 + '.tmp.mp4') < MAXBYTES:
            break
    os.replace(mp4 + '.tmp.mp4', mp4)
    shutil.rmtree(tmp)
    poster = os.path.join(od, cfg['movie'] + sfx + '_poster.png')
    shutil.copy(os.path.join(fdir, 'f_%05d.png' % recs[-1]['n']), poster)
    info = dict(run=run, nframes=nf, fps=fps, crf=crf, size=os.path.getsize(mp4), t0_d=recs[0]['t']/DAY,
                t1_d=recs[-1]['t']/DAY, clim=clim, cmap_teff=cfg['cmap_teff'], built=time.strftime('%Y-%m-%d %H:%M'))
    json.dump(info, open(os.path.join(od, cfg['movie'] + sfx + '.json'), 'w'), indent=1)
    print('encoded %s: %d frames %.2f-%.2f d, %d fps, crf %d, %.2f MB' % (mp4, nf, info['t0_d'], info['t1_d'], fps, crf,
                                                                          info['size']/1e6), flush=True)
    pub = cfg.get('publish')
    if pub:
        os.makedirs(pub, exist_ok=True)
        for src in (mp4, poster):
            dst = os.path.join(pub, os.path.basename(src))
            shutil.copy(src, dst + '.tmp')
            os.replace(dst + '.tmp', dst)
        print('copied to', pub, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', choices=sorted(RUNS))
    ap.add_argument('--only', default='')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--no-encode', action='store_true')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--tag', default='')
    ap.add_argument('--set', default='{}')
    ap.add_argument('--view', default='single', choices=sorted(VIEWS))
    a = ap.parse_args()
    run = a.run
    RUNS[run].update(json.loads(a.set))
    RUNS[run]['view'] = a.view
    cfg = RUNS[run]
    sfx = VIEWS[a.view]['suffix']
    os.makedirs(os.path.join(VIS, run, 'frames' + sfx), exist_ok=True)
    dl = dump_list(cfg)
    if not dl:
        print('no complete dumps'); return
    import multiprocessing as mp
    if a.check:
        have = {int(os.path.basename(f)[2:7]) for f in glob.glob(os.path.join(VIS, run, 'frames' + sfx, 'f_*.json'))}
        print(len([q for q in dl if q[0] not in have])); return
    ctx = mp.get_context('spawn')
    # 1) reductions (cached), 2) colour ranges + view radius, 3) frames
    with ctx.Pool(max(1, min(a.workers, int(os.environ.get('RED_CAP', 4))))) as pool:   # reductions: ~3 GB each, cap the parallelism (RED_CAP env, 10-22)
        st = pool.map(red_one, [(run, n, fh, key) for n, fh, key in dl], chunksize=1)
    cfg['relief'] = 1.0   # USER 10-27: never exaggerate (auto_relief had squashed the real +-20-40 % relief to 4 %)
    n_new, fh_new, key_new = dl[-1]
    clim = colour_limits(run, cfg, n_new, fh_new, key_new)
    rview = max(m + cfg['relief']*(hi - m) for _, _, m, lo, hi in st)/RSUN*1.06
    clim['rview'] = float(np.ceil(rview/10.0)*10.0)
    clim['relief'] = cfg['relief']                      # spawned workers re-import RUNS: carry it in clim
    clim = sticky_clim(os.path.join(VIS, run, 'clim%s.json' % sfx), clim)
    h = hashlib.md5((json.dumps({k: v for k, v in cfg.items() if k != 'publish'}, sort_keys=True) + json.dumps(clim)
                     + RENDER_VERSION + repr(SIZE)).encode()).hexdigest()[:10]
    print('%s: %d dumps, clim %s (cfg %s)' % (run, len(dl), json.dumps(clim), h), flush=True)
    if a.only:
        want = {int(x) for x in a.only.split(',')}
        dl = [q for q in dl if q[0] in want]
    todo, recs = [], []
    for n, fh, key in dl:
        fkey = key + '|' + h
        js = os.path.join(VIS, run, 'frames' + sfx, 'f_%05d.json' % n)
        if not a.tag and os.path.exists(js) and os.path.exists(js[:-5] + '.png'):
            q = json.load(open(js))
            if q.get('key') == fkey:
                recs.append(q)
                continue
        todo.append((run, n, fh, key, fkey, clim, a.tag, a.view))
    print('%d cached frames, %d to render' % (len(recs), len(todo)), flush=True)
    if todo:
        with ctx.Pool(max(1, a.workers)) as pool:
            recs += pool.map(work, todo, chunksize=1)
    if a.no_encode or a.tag or a.only or not recs:
        return
    encode(run, cfg, recs, clim)


if __name__ == '__main__':
    main()
