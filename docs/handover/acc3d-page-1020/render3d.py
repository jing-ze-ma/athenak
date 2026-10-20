"""ACC3D page 10-20 (bundle copy of page/mov3d/render3d.py): 3-D isosurface frames (pvbatch or python + vtk>=9, off-screen).
Title text from env SITE / RUN_LABEL. Pure VTK (no paraview.simple), so a pip vtk wheel with OSMesa/EGL works too.
usage: pvbatch render3d.py <red_dir> <frame_dir> <cfg.json> [i/n]   renders frame_dir/f_NNNNN.png for every red_NNNNN.npz
without a frame (FORCE=1: all; 'i/n' = this worker renders every n-th file starting at i).
Nested density isosurfaces (cfg levels, log10 rho code) of the r <= d_L1 region, coloured by log10 T_gas (EOS-table inversion),
fixed colour bar + camera; grey inner boundary (r_in = 6.3 Rsun) band."""
import os
import sys
import glob
import json
import numpy as np
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
import vtkmodules.vtkRenderingFreeType  # noqa: F401
from vtkmodules.util import numpy_support as ns
from vtkmodules.vtkCommonCore import vtkPoints, vtkLookupTable
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid
from vtkmodules.vtkFiltersCore import vtkContourFilter, vtkPolyDataNormals
from vtkmodules.vtkFiltersGeometry import vtkStructuredGridGeometryFilter
from vtkmodules.vtkRenderingCore import (vtkRenderer, vtkRenderWindow, vtkPolyDataMapper, vtkActor,
                                         vtkWindowToImageFilter, vtkTextActor, vtkLight)
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
from vtkmodules.vtkIOImage import vtkPNGWriter

RED, FR, CFG = sys.argv[1], sys.argv[2], sys.argv[3]
WI, WN = (int(v) for v in sys.argv[4].split('/')) if len(sys.argv) > 4 else (0, 1)
C = json.load(open(CFG))
os.makedirs(FR, exist_ok=True)
FORCE = os.environ.get('FORCE', '') == '1'
DL1 = 16.2266
RIN = 6.3
W, H = C['size']


def lut(lo, hi):
    # inferno-like sequential (dark purple -> red -> orange -> pale yellow), perceptually ordered
    st = np.array([[0.00, 0.001, 0.000, 0.014], [0.15, 0.159, 0.044, 0.329], [0.30, 0.397, 0.083, 0.433],
                   [0.45, 0.623, 0.165, 0.388], [0.60, 0.831, 0.283, 0.259], [0.75, 0.962, 0.490, 0.084],
                   [0.88, 0.982, 0.730, 0.200], [1.00, 0.988, 0.998, 0.645]])
    lt = vtkLookupTable()
    lt.SetNumberOfTableValues(256)
    lt.SetTableRange(lo, hi)
    for i, s in enumerate(np.linspace(0, 1, 256)):
        c = [np.interp(s, st[:, 0], st[:, k]) for k in (1, 2, 3)]
        lt.SetTableValue(i, c[0], c[1], c[2], 1.0)
    lt.SetBelowRangeColor(*lt.GetTableValue(0)); lt.SetAboveRangeColor(*lt.GetTableValue(255))
    lt.Build()
    return lt


LT = lut(*C['logT'])
GEO = {}


def grid(d):
    """structured grid of cell centres, phi wrapped (seam closed), one extra radial layer at r = d_L1 with floor density
    so the surfaces close on the d_L1 sphere."""
    r = np.concatenate([d['r'], [DL1]]).astype(np.float64)
    th = np.asarray(d['th'], np.float64); ph = np.asarray(d['ph'], np.float64)
    ph = np.concatenate([ph, [ph[0] + 2 * np.pi]])
    key = (len(r), len(th), len(ph))
    if key not in GEO:
        P, T, R = np.meshgrid(ph, th, r, indexing='ij')           # (phi, theta, r): r fastest = VTK i
        xyz = np.stack([R * np.sin(T) * np.cos(P), R * np.sin(T) * np.sin(P), R * np.cos(T)], axis=-1).reshape(-1, 3)
        pts = vtkPoints(); pts.SetData(ns.numpy_to_vtk(np.ascontiguousarray(xyz, np.float32), deep=1))
        GEO[key] = pts
    def pad(q, fill):
        q = np.asarray(q, np.float32)
        q = np.concatenate([q, q[:1]], axis=0)                     # phi wrap
        if fill is None:
            q = np.concatenate([q, q[:, :, -1:]], axis=2)
        else:
            q = np.concatenate([q, np.full(q[:, :, -1:].shape, fill, np.float32)], axis=2)
        return np.ascontiguousarray(q.reshape(-1))
    g = vtkStructuredGrid()
    g.SetDimensions(len(r), len(th), len(ph))
    g.SetPoints(GEO[key])
    a = ns.numpy_to_vtk(pad(d['lrho'], -9.0), deep=1); a.SetName('lrho')
    b = ns.numpy_to_vtk(pad(d['lT'], None), deep=1); b.SetName('lT')
    g.GetPointData().AddArray(a); g.GetPointData().AddArray(b)
    g.GetPointData().SetActiveScalars('lrho')
    return g, th, ph


def iso(g, level, op):
    cf = vtkContourFilter()
    cf.SetInputData(g)
    cf.SetValue(0, level)
    cf.SetInputArrayToProcess(0, 0, 0, 0, 'lrho')
    cf.ComputeScalarsOff()
    cf.ComputeNormalsOff()
    nf = vtkPolyDataNormals()
    nf.SetInputConnection(cf.GetOutputPort())
    nf.SplittingOff(); nf.ConsistencyOn()
    m = vtkPolyDataMapper()
    m.SetInputConnection(nf.GetOutputPort())
    m.SetScalarModeToUsePointFieldData(); m.SelectColorArray('lT')
    m.SetLookupTable(LT); m.SetScalarRange(LT.GetTableRange()); m.UseLookupTableScalarRangeOn()
    m.ScalarVisibilityOn(); m.InterpolateScalarsBeforeMappingOn()
    a = vtkActor(); a.SetMapper(m)
    p = a.GetProperty()
    p.SetOpacity(op); p.SetAmbient(0.35); p.SetDiffuse(0.75); p.SetSpecular(0.15 if op > 0.9 else 0.05); p.SetSpecularPower(20)
    p.BackfaceCullingOff()
    return a


def inner_band(th, ph):
    tt = np.linspace(th[0] - 0.5 * (th[1] - th[0]), th[-1] + 0.5 * (th[-1] - th[-2]), 30)
    pp = np.linspace(0, 2 * np.pi, 257)
    P, T = np.meshgrid(pp, tt, indexing='ij')
    R = RIN * 0.999
    xyz = np.stack([R * np.sin(T) * np.cos(P), R * np.sin(T) * np.sin(P), R * np.cos(T)], axis=-1).reshape(-1, 3)
    g = vtkStructuredGrid(); g.SetDimensions(len(tt), len(pp), 1)
    pts = vtkPoints(); pts.SetData(ns.numpy_to_vtk(np.ascontiguousarray(xyz, np.float32), deep=1)); g.SetPoints(pts)
    f = vtkStructuredGridGeometryFilter(); f.SetInputData(g); f.Update()
    m = vtkPolyDataMapper(); m.SetInputData(f.GetOutput()); m.ScalarVisibilityOff()
    a = vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.55, 0.55, 0.58); a.GetProperty().SetAmbient(0.3)
    return a


def ring(radius, rgb, z=0.0):
    """thin reference circle in the orbital plane (d_L1 sphere outline)."""
    from vtkmodules.vtkFiltersSources import vtkRegularPolygonSource
    s = vtkRegularPolygonSource(); s.SetNumberOfSides(256); s.SetRadius(radius); s.SetCenter(0, 0, z); s.GeneratePolygonOff()
    m = vtkPolyDataMapper(); m.SetInputConnection(s.GetOutputPort())
    a = vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*rgb); a.GetProperty().SetOpacity(0.5); a.GetProperty().SetLineWidth(1.5)
    return a


def text(s, x, y, size, just='left', rgb=(1, 1, 1)):
    ta = vtkTextActor(); ta.SetInput(s)
    tp = ta.GetTextProperty(); tp.SetFontSize(size); tp.SetColor(*rgb)
    {'left': tp.SetJustificationToLeft, 'center': tp.SetJustificationToCentered, 'right': tp.SetJustificationToRight}[just]()
    ta.SetDisplayPosition(x, y)
    return ta


REDS = sorted(glob.glob(os.path.join(RED, 'red_*.npz')))[WI::WN]
for f in REDS:
    idx = f.split('_')[-1].split('.')[0]
    out = os.path.join(FR, 'f_%s.png' % idx)
    if os.path.exists(out) and not FORCE:
        continue
    d = np.load(f)
    g, th, ph = grid(d)
    ren = vtkRenderer()
    ren.SetBackground(0.02, 0.02, 0.05)
    ren.SetUseDepthPeeling(1); ren.SetMaximumNumberOfPeels(12); ren.SetOcclusionRatio(0.0)
    ren.AddActor(inner_band(th, ph))
    for L, op in zip(C['levels'], C['opacity']):
        ren.AddActor(iso(g, L, op))
    ren.AddActor(ring(DL1, (0.6, 0.8, 1.0)))
    cam = ren.GetActiveCamera()
    fx, fy, fz = C['focus']
    el, az = np.radians(C['cam_elev']), np.radians(C['cam_az'])
    D = C['cam_dist']
    cam.SetFocalPoint(fx, fy, fz)
    cam.SetPosition(fx + D * np.cos(el) * np.cos(az), fy + D * np.cos(el) * np.sin(az), fz + D * np.sin(el))
    cam.SetViewUp(0, 0, 1); cam.SetViewAngle(C['view_angle'])
    ren.ResetCameraClippingRange()
    t = float(d['t'])
    ts = t - C['t_stream']
    ren.AddActor2D(text('RY Per accretor, 3-D: %s (%s)' % (os.environ.get('RUN_LABEL', 'production'), os.environ.get('SITE', 'Raven')), 20, H - 40, 22))
    ren.AddActor2D(text('t = %.3f code = %.3f P_orb   (%s)' % (t, t / C['porb'],
                        'static phase, stream off' if ts < -1e-9 else 'stream on for %.3f P_orb' % (ts / C['porb'])), 20, H - 70, 20))
    ren.AddActor2D(text('isosurfaces log10 rho [code, 5.11e-8 g/cc] = ' + ', '.join('%g' % L for L in C['levels'])
                        + ' (opacity ' + ', '.join('%g' % o for o in C['opacity']) + ');  r <= d_L1 = 16.23 Rsun (blue ring)',
                        20, 38, 16, rgb=(0.85, 0.85, 0.85)))
    ren.AddActor2D(text('colour: gas temperature from the run\'s general-EOS table (rho, e_int); grey = inner boundary r = 6.3 Rsun;'
                        ' donor along +x (towards the lower right)', 20, 14, 16, rgb=(0.85, 0.85, 0.85)))
    sb = vtkScalarBarActor()
    sb.SetLookupTable(LT); sb.SetTitle('log10 T_gas [K]'); sb.SetNumberOfLabels(7); sb.SetLabelFormat('%.1f')
    sb.SetPosition(0.905, 0.18); sb.SetWidth(0.07); sb.SetHeight(0.64)
    sb.GetTitleTextProperty().SetFontSize(16); sb.GetLabelTextProperty().SetFontSize(15)
    sb.SetUnconstrainedFontSize(True); sb.SetDrawBelowRangeSwatch(False); sb.SetDrawAboveRangeSwatch(False)
    ren.AddActor2D(sb)
    rw = vtkRenderWindow()
    rw.SetOffScreenRendering(1); rw.SetAlphaBitPlanes(1); rw.SetMultiSamples(0)
    rw.SetSize(W, H)
    rw.AddRenderer(ren)
    rw.Render()
    wi = vtkWindowToImageFilter(); wi.SetInput(rw); wi.Update()
    pw = vtkPNGWriter(); pw.SetFileName(out + '.tmp.png'); pw.SetInputConnection(wi.GetOutputPort()); pw.Write()
    os.replace(out + '.tmp.png', out)
    print('frame', idx, 't %.4f' % t, flush=True)
