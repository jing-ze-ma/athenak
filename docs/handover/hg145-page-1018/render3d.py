"""HGPAGE 10-18: 3-D movie frames of the hg145 wedge (pvbatch, pure VTK, EGL off-screen).
usage: pvbatch render3d.py <red_dir> <frame_dir> <cfg.json>   renders frame_dir/f_NNNNN.png for every red_NNNNN.npz without a frame
(FORCE=1: all). Photosphere R_ph(theta,phi) coloured by T_eff; the wedge side faces (phi = 0, pi/3; theta = pi/3, 2pi/3) from
r_in to the local R_ph coloured by v_r."""
import os
import sys
import glob
import json
import numpy as np
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
import vtkmodules.vtkRenderingFreeType  # noqa: F401
from vtkmodules.vtkCommonCore import vtkPoints, vtkFloatArray, vtkLookupTable
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid
from vtkmodules.vtkFiltersGeometry import vtkStructuredGridGeometryFilter
from vtkmodules.vtkRenderingCore import (vtkRenderer, vtkRenderWindow, vtkPolyDataMapper, vtkActor,
                                         vtkWindowToImageFilter, vtkTextActor)
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
from vtkmodules.vtkIOImage import vtkPNGWriter

RED, FR, CFG = sys.argv[1], sys.argv[2], sys.argv[3]
C = json.load(open(CFG))
os.makedirs(FR, exist_ok=True)
RS = 6.957e10
FORCE = os.environ.get('FORCE', '') == '1'


def lut(name, lo, hi):
    lt = vtkLookupTable()
    lt.SetNumberOfTableValues(256)
    lt.SetTableRange(lo, hi)
    x = np.linspace(0, 1, 256)
    if name == 'coolwarm':   # diverging blue - white - red (Moreland-like)
        c0, cm, c1 = np.array([59, 76, 192]) / 255., np.array([0.93, 0.93, 0.93]), np.array([180, 4, 38]) / 255.
        for i, s in enumerate(x):
            c = c0 + (cm - c0) * (s / 0.5) if s < 0.5 else cm + (c1 - cm) * ((s - 0.5) / 0.5)
            lt.SetTableValue(i, c[0], c[1], c[2], 1.0)
    else:                    # warm sequential for T_eff: dark red -> orange -> pale yellow
        stops = np.array([[0.0, 0.18, 0.02, 0.02], [0.35, 0.70, 0.15, 0.05], [0.7, 0.98, 0.55, 0.12], [1.0, 1.0, 0.95, 0.75]])
        for i, s in enumerate(x):
            c = [np.interp(s, stops[:, 0], stops[:, k]) for k in (1, 2, 3)]
            lt.SetTableValue(i, c[0], c[1], c[2], 1.0)
    lt.Build()
    return lt


def sgrid(X, Y, Z, S):
    """structured surface from 2-D arrays (n0, n1)."""
    n0, n1 = X.shape
    p = vtkPoints()
    p.SetNumberOfPoints(n0 * n1)
    a = vtkFloatArray()
    a.SetName('s')
    a.SetNumberOfTuples(n0 * n1)
    k = 0
    for j in range(n0):
        for i in range(n1):
            p.SetPoint(k, X[j, i], Y[j, i], Z[j, i])
            a.SetValue(k, S[j, i])
            k += 1
    g = vtkStructuredGrid()
    g.SetDimensions(n1, n0, 1)
    g.SetPoints(p)
    g.GetPointData().SetScalars(a)
    f = vtkStructuredGridGeometryFilter()
    f.SetInputData(g)
    f.Update()
    return f.GetOutput()


def actor(poly, lt):
    m = vtkPolyDataMapper()
    m.SetInputData(poly)
    m.SetLookupTable(lt)
    m.SetScalarRange(lt.GetTableRange())
    m.SetScalarModeToUsePointData()
    m.InterpolateScalarsBeforeMappingOn()
    a = vtkActor()
    a.SetMapper(m)
    a.GetProperty().SetAmbient(0.25)
    a.GetProperty().SetDiffuse(0.8)
    return a


def sph(r, th, ph):
    return r * np.sin(th) * np.cos(ph), r * np.sin(th) * np.sin(ph), r * np.cos(th)


def face(rr, Rcap, ang_fixed, angs, vr, kind):
    """side face over (angle, r) with r clipped at the local photosphere radius Rcap(angle)."""
    R2 = np.minimum(rr[None, :], Rcap[:, None])
    A2 = np.repeat(angs[:, None], len(rr), axis=1)
    if kind == 'phi':      # phi fixed, angs = theta
        X, Y, Z = sph(R2, A2, ang_fixed)
    else:                  # theta fixed, angs = phi
        X, Y, Z = sph(R2, ang_fixed, A2)
    return sgrid(X, Y, Z, vr)


# colour scales FIXED over the whole movie: 'auto' = robust percentiles over ALL reduced dumps (|v_r| on the side faces
# below the local R_ph: 99.5 %; T_eff: 0.5-99.5 %). When the resolved scale moves by > 5 % (run grew), all frames are re-rendered.
REDS = sorted(glob.glob(os.path.join(RED, 'red_*.npz')))
if C['vr'] == 'auto' or C['teff'] == 'auto':
    av, at = [], []
    for f in REDS:
        d = np.load(f)
        for k in ('f_p0_vr', 'f_p1_vr', 'f_t0_vr', 'f_t1_vr'):
            av.append(np.abs(np.asarray(d[k], np.float32)).ravel())
        at.append(np.asarray(d['Teff'], np.float32).ravel() / 1e3)
    av = np.concatenate(av); at = np.concatenate(at)
    VR = float(np.nanpercentile(av, 99.5)) if C['vr'] == 'auto' else float(C['vr'])
    TE = [float(np.nanpercentile(at, 0.5)), float(np.nanpercentile(at, 99.5))] if C['teff'] == 'auto' else C['teff']
    if TE[1] - TE[0] < 0.02:
        mid = 0.5 * (TE[0] + TE[1]); TE = [mid - 0.01, mid + 0.01]
else:
    VR, TE = float(C['vr']), C['teff']
SC = os.path.join(FR, 'scale.json')
try:
    old = json.load(open(SC))
    if abs(old['vr'] / VR - 1) > 0.05 or abs((old['teff'][1] - old['teff'][0]) / (TE[1] - TE[0]) - 1) > 0.05 \
            or abs(old['teff'][0] - TE[0]) > 0.05 * (TE[1] - TE[0]):
        FORCE = True
    else:
        VR, TE = old['vr'], old['teff']     # keep the scale the existing frames use
except Exception:
    FORCE = True
json.dump({'vr': VR, 'teff': TE}, open(SC, 'w'))
print('scale v_r +-%.3g km/s, T_eff %.4g-%.4g kK, re-render all: %s' % (VR, TE[0], TE[1], FORCE), flush=True)
lt_v = lut('coolwarm', -VR, VR)
lt_T = lut('warm', TE[0], TE[1])
W, H = C.get('size', [1000, 760])

for f in REDS:
    idx = f.split('_')[-1].split('.')[0]
    out = os.path.join(FR, 'f_%s.png' % idx)
    if os.path.exists(out) and not FORCE:
        continue
    d = np.load(f)
    r = d['r'] / RS
    ns = int(d['ns'])
    rr = r[:ns]
    th, ph = d['th'], d['ph']
    Rph = d['Rph'] / RS
    med = np.nanmedian(Rph)
    Rph = np.where(np.isfinite(Rph), Rph, med)
    Te = d['Teff'] / 1e3
    Te = np.where(np.isfinite(Te), Te, np.nanmedian(Te))
    # extend the cap to the wedge edges (faces): use cell-centre angles padded to the face angles
    thp = np.concatenate([[th[0] - 0.5 * (th[1] - th[0])], th, [th[-1] + 0.5 * (th[-1] - th[-2])]])
    php = np.concatenate([[ph[0] - 0.5 * (ph[1] - ph[0])], ph, [ph[-1] + 0.5 * (ph[-1] - ph[-2])]])
    Rp = np.pad(Rph, 1, mode='edge')       # (k, j) = (phi, theta)
    Tp = np.pad(Te, 1, mode='edge')
    PH, TH = np.meshgrid(php, thp, indexing='ij')
    X, Y, Z = sph(Rp, TH, PH)
    ren = vtkRenderer()
    ren.SetBackground(0, 0, 0)
    ren.AddActor(actor(sgrid(X, Y, Z, Tp), lt_T))
    vr = {k: np.asarray(d[k], float) for k in ('f_p0_vr', 'f_p1_vr', 'f_t0_vr', 'f_t1_vr')}
    ren.AddActor(actor(face(rr, Rp[0, :], php[0], thp, np.pad(vr['f_p0_vr'], ((1, 1), (0, 0)), mode='edge'), 'phi'), lt_v))
    ren.AddActor(actor(face(rr, Rp[-1, :], php[-1], thp, np.pad(vr['f_p1_vr'], ((1, 1), (0, 0)), mode='edge'), 'phi'), lt_v))
    ren.AddActor(actor(face(rr, Rp[:, 0], thp[0], php, np.pad(vr['f_t0_vr'], ((1, 1), (0, 0)), mode='edge'), 'theta'), lt_v))
    ren.AddActor(actor(face(rr, Rp[:, -1], thp[-1], php, np.pad(vr['f_t1_vr'], ((1, 1), (0, 0)), mode='edge'), 'theta'), lt_v))
    cam = ren.GetActiveCamera()
    fx, fy, fz = sph(C['focus_r'], np.pi / 2, np.pi / 6)
    cam.SetFocalPoint(fx, fy, fz)
    px, py, pz = sph(C['cam_r'], np.radians(C['cam_th']), np.radians(C['cam_ph']))
    cam.SetPosition(px, py, pz)
    cam.SetViewUp(0, 0, 1)
    cam.SetViewAngle(C.get('view_angle', 30))
    ren.ResetCameraClippingRange()
    t = float(d['t']) / 86400
    ta = vtkTextActor()
    ta.SetInput('He giant 2.75 Msun He star (production, %s)\nt = %.2f d     median R_ph = %.1f R_sun' % (os.environ.get("SITE", "Raven"), t, med))
    tp = ta.GetTextProperty()
    tp.SetFontSize(26); tp.SetJustificationToCentered(); tp.SetColor(1, 1, 1)
    ta.SetDisplayPosition(W // 2, H - 80)
    ren.AddActor2D(ta)
    sb = vtkScalarBarActor()
    sb.SetLookupTable(lt_T); sb.SetTitle('T_eff (kK)'); sb.SetNumberOfLabels(3); sb.SetLabelFormat('%.2f')
    sb.SetPosition(0.86, 0.2); sb.SetWidth(0.08); sb.SetHeight(0.6)
    sb.GetTitleTextProperty().SetFontSize(18); sb.GetLabelTextProperty().SetFontSize(16)
    sb.SetUnconstrainedFontSize(True)
    ren.AddActor2D(sb)
    sv = vtkScalarBarActor()
    sv.SetLookupTable(lt_v); sv.SetTitle('radial velocity (km/s)'); sv.SetOrientationToHorizontal(); sv.SetNumberOfLabels(3)
    sv.SetLabelFormat('%.1f')
    sv.SetPosition(0.22, 0.02); sv.SetWidth(0.5); sv.SetHeight(0.11)
    sv.GetTitleTextProperty().SetFontSize(18); sv.GetLabelTextProperty().SetFontSize(16)
    sv.SetUnconstrainedFontSize(True)
    ren.AddActor2D(sv)
    rw = vtkRenderWindow()
    rw.SetOffScreenRendering(1)
    rw.SetSize(W, H)
    rw.AddRenderer(ren)
    rw.Render()
    wi = vtkWindowToImageFilter()
    wi.SetInput(rw)
    wi.Update()
    pw = vtkPNGWriter()
    pw.SetFileName(out)
    pw.SetInputConnection(wi.GetOutputPort())
    pw.Write()
    print('frame', idx, 't %.2f d' % t, flush=True)
