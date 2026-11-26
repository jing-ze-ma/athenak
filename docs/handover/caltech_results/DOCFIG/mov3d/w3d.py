"""w3d.py STAR [--only a-b] [--workers W] [--variant notch|face]: 3-D wedge movie frames, style of vis3d/render3d.py v7 (functions imported),
colours of the USER's 11-25 pick (T_eff cmocean ice, v_r cmasher fusion), 1920 x 1080.  Existing dumps only, read only."""
import os, sys, json, time, argparse, re, glob
os.environ['VTK_DEFAULT_OPENGL_WINDOW'] = 'vtkEGLRenderWindow'
import numpy as np
G = '/resnick/groups/carnegie_poc/jingze'; D0 = G + '/mov3d_1125'
sys.path.insert(0, G + '/vis3d'); sys.dont_write_bytecode = True
import render3d as r3
for k in [m for m in sys.modules if m == 'bin_convert']: del sys.modules[k]
import types; sys.modules.setdefault('h5py', types.ModuleType('h5py'))   # reader imports h5py but does not need it here
sys.path.insert(0, G + '/rsgE_page/vispy_4de8168d'); import bin_convert as bc; r3.bc = bc
from matplotlib.colors import ListedColormap
ICE, FUS = (ListedColormap(np.load(D0 + '/lut_%s.npy' % n)) for n in ('ice', 'fusion'))
RS, DAY = r3.RSUN, 86400.0
RUNS = dict(bsg=dict(title='M20 BSG', hydro=G + '/prod_c1_1124/bsg/bin/bsg3d.hydro_w.%05d.bin', rad='vet', var='vet_f1', grid=G + '/prod_c1_1124/bsg/bsg3d.x1grid.txt',
    table=G + '/fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt', mu=0.62, gamma=5/3., rin=0.0, rview=56.0, az=30, el=30, zoom=1.17))
TB2 = G + '/fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
def _c(title, d, pre, mu, rin, rview, **k):
    q = dict(title=title, hydro=G + '/' + d + '/bin/' + pre + '.hydro_w.%05d.bin', rad='vet', var='vet_f1', grid=G + '/' + d + '/' + pre + '.x1grid.txt', table=TB2, mu=mu, gamma=5/3., rin=rin, rview=rview, az=30, el=30, zoom=1.0, fit=True, last=None, te=None); q.update(k); return q
RUNS['bsg'].update(last=114, te=None)
RUNS['m35mms'] = _c('M35 MMS', 'prod_c1_1124/m35mms', 'm35mms', 0.6215, 0.86*9.554e11/6.957e10, 14.2, last=480, note='frames interpolated')
RUNS['m35zams'] = _c('M35 ZAMS', 'prod_c1_1124/m35zams', 'm35zams', 0.6220, 0.93*5.4972e11/6.957e10, 8.0, last=480, note='frames interpolated')
RUNS['t9l62'] = _c('M56 LBV (radiative IC)', 'prod_c4_1124/t9l62', 't9l62', 0.674, 0.42*6.5326e12/6.957e10, 560.0, last=249, te=[4.9, 10.0], amb=0.80, dif=0.30, spec=0.0)
RUNS['t9mlt'] = _c('M56 LBV', 'prod_mlt_1125/t9mlt', 't9l62', 0.674, 0.62*3.9282e12/6.957e10, 64.0, last=178, te=[9.0, 29.5], amb=0.80, dif=0.30, spec=0.0)
def lastdump(s):                                               # newest dump complete for > 10 min (the run may still be writing)
    c = RUNS[s]
    if c['last'] is not None: return c['last']
    P = []
    for f in sorted(glob.glob(c['hydro'].replace('%05d', '*'))):
        g = f.replace('hydro_w', c['rad'])
        if os.path.exists(g) and time.time() - max(os.path.getmtime(f), os.path.getmtime(g)) > 600: P.append(int(f.split('.')[-2]))
    return max(P)
SIZE = (1920, 1080)
for _s, _L, _R in (('bsg', 10., 52.5), ('m35mms', 1., 13.65), ('m35zams', 0.1, 7.90), ('t9l62', 100., 470.), ('t9mlt', 10., 55.)): RUNS[_s].update(bar=_L, rref=_R, key=_s)
def basis(cam):
    ep, fp, up = (np.array(v, float) for v in (cam.position, cam.focal_point, cam.up)); fwd = (fp-ep)/np.linalg.norm(fp-ep); rt = np.cross(fwd, up); rt /= np.linalg.norm(rt); return ep, fwd, rt, np.cross(rt, fwd)
def project(cam, q):                                           # pixel (x from left, y from bottom) of world points [R_sun]
    ep, fwd, rt, uu = basis(cam); d = q-ep; dep = d@fwd; th = np.tan(np.radians(cam.view_angle)/2); asp = SIZE[0]/SIZE[1]; wc = np.array(cam.GetWindowCenter())
    return np.stack([((d@rt)/dep/(th*asp)-wc[0]+1)/2*SIZE[0], ((d@uu)/dep/th-wc[1]+1)/2*SIZE[1]], -1), dep
def addbar(pv, pl, D, c, caps):
    cam = pl.camera; ep, fwd, rt, uu = basis(cam); th = np.tan(np.radians(cam.view_angle)/2); asp = SIZE[0]/SIZE[1]; wc = np.array(cam.GetWindowCenter()); Q = []
    for k0, k1, j0, j1 in caps:                                # depth of the point of the shown photosphere (radius rref, fixed over the movie) nearest to the camera
        PH, TH = np.meshgrid(D['phf'][k0:k1+1], D['thf'][j0:j1+1], indexing='ij'); Q.append((r3.sph(np.full(PH.shape, c['rref']), TH, PH).reshape(-1, 3)-ep)@fwd)
    dep = float(np.concatenate(Q).min()); ln = c['bar']/(dep*th*asp)*SIZE[0]/2; _m, py0 = json.load(open(D0 + '/barpos.json'))[c['key']]; _m += 163.0 if os.environ.get('HERO') else 0.0; px0 = _m-ln/2; X = (2*px0/SIZE[0]-1+wc[0])*th*asp*dep; Y = (2*py0/SIZE[1]-1+wc[1])*th*dep; P0 = ep+fwd*dep+rt*X+uu*Y; P1 = P0+rt*c['bar']
    ln = c['bar']/(dep*th*asp)*SIZE[0]/2
    for a, b in ((P0, P1),) + tuple((p-uu*8*dep*th*2/SIZE[1], p+uu*8*dep*th*2/SIZE[1]) for p in (P0, P1)): pl.add_mesh(pv.Line(a, b), color='white', line_width=5, lighting=False)
    lb = pl.add_text('%g R_sun' % c['bar'], position=(int(round(px0+ln/2)), int(py0+14)), font_size=20 if os.environ.get('HERO') else 15, color='white'); lb.prop.SetJustificationToCentered(); pl.renderer.ResetCameraClippingRange()
    return ln, dep
def assemble(d, names):                                            # float64 version of render3d.assemble
    lg = np.array(d['mb_logical'])[:, :3]; n3, n2, n1 = d['mb_data'][names[0]][0].shape
    out = {k: np.empty((n3*(lg[:, 2].max()+1), n2*(lg[:, 1].max()+1), n1), np.float64) for k in names}; thf, phf = np.empty(n2*(lg[:, 1].max()+1)+1), np.empty(n3*(lg[:, 2].max()+1)+1)
    for b in range(d['n_mbs']):
        j, k = lg[b, 1], lg[b, 2]
        for nm in names: out[nm][k*n3:(k+1)*n3, j*n2:(j+1)*n2] = d['mb_data'][nm][b]
        thf[j*n2:(j+1)*n2+1] = d['mb_x2f'][b]; phf[k*n3:(k+1)*n3+1] = d['mb_x3f'][b]
    return out, np.array(d['mb_x1f'][0], float), thf, phf
def reduce(s, n):
    f = D0 + '/red/%s_%05d.npz' % (s, n)
    if os.path.exists(f):
        z = np.load(f); return {k: (z[k][()] if z[k].ndim == 0 else z[k]) for k in z.files}
    c = RUNS[s]; tab = r3.read_ross(c['table']); fh = c['hydro'] % n; a = bc.read_binary(fh); t = float(a['time']); A, rf, thf, phf = assemble(a, ['dens', 'velx', 'eint']); del a
    gr = np.array([x.split()[:4] for x in open(c['grid']) if not x.startswith('#') and len(x.split()) == 4], float); nr = A['dens'].shape[2]; rff = gr[3:3+nr+1, 1]
    if len(rff) == nr + 1: rf = rff                               # stretched grid: true faces from x1grid.txt (the dump header assumes uniform)
    T = (c['gamma']-1)*c['mu']*r3.MU_K*A['eint']/A['dens']; rho, vr = A['dens'], A['velx']; del A
    F1 = assemble(bc.read_binary(fh.replace('hydro_w', c['rad'])), [c['var']])[0][c['var']]
    nph, nth, nr = rho.shape; dr = np.diff(rf); rc = r3.centroid(rf[:-1], rf[1:]); kr = r3.kappa(tab, rho, T)*rho; tauf = np.zeros((nph, nth, nr+1)); tauf[..., :nr] = np.flip(np.cumsum(np.flip(kr*dr, 2), 2), 2)
    ge = tauf[..., :nr] >= 2/3.; i = nr-1-np.argmax(ge[..., ::-1], -1); ta = np.take_along_axis(tauf, (i+1)[..., None], -1)[..., 0]; k = np.take_along_axis(kr, i[..., None], -1)[..., 0]
    rph = np.where(ge[..., 0], rf[i+1]-(2/3.-ta)/np.maximum(k, 1e-300), rf[0]); idx = np.where(ge[..., 0], i, 0)
    ic = np.clip(np.searchsorted(rc, rph)-1, 0, nr-2); w = np.clip((rph-rc[ic])/(rc[ic+1]-rc[ic]), 0, 1); at = lambda Q: (1-w)*np.take_along_axis(Q, ic[..., None], -1)[..., 0]+w*np.take_along_axis(Q, (ic+1)[..., None], -1)[..., 0]
    o = dict(t=t, rf=rf, thf=thf, phf=phf, rph=rph, idx=idx, Fph=at(F1)); km, jm = nph//2, nth//2
    for q in sorted({0, km-1, nph-1}): o['fk_%d' % q] = (vr[q]/1e5).astype(np.float32)
    for q in sorted({0, jm, nth-1}): o['fj_%d' % q] = (vr[:, q]/1e5).astype(np.float32)
    np.savez(f + '.%d.tmp.npz' % os.getpid(), **o); os.replace(f + '.%d.tmp.npz' % os.getpid(), f); return o
def trim(D, rin):                                                 # cut faces from the envelope base outward only
    i0 = int(np.searchsorted(D['rf'], rin*RS)); D = dict(D); D['rf'] = D['rf'][i0:]; D['idx'] = D['idx']-i0
    for k in list(D):
        if k[:3] in ('fk_', 'fj_'): D[k] = D[k][:, i0:]
    return D
def render(D, c, L, out, variant):
    import pyvista as pv
    pv.OFF_SCREEN = True; D = trim(D, c['rin']); cfg = dict(cap='Teff', cap_points=True, teff_smooth=0, rph_smooth=2); nph, nth = D['rph'].shape; km, jm = nph//2, nth//2
    pl = pv.Plotter(off_screen=True, window_size=SIZE, lighting='none'); pl.background_color = 'black'
    pl.add_light(pv.Light(light_type='headlight', intensity=0.55)); pl.add_light(pv.Light(position=(1.0, -0.3, 1.2), focal_point=(0, 0, 0), light_type='scene light', intensity=0.55))
    if variant == 'notch':
        fs = [r3.face_phi(pv, D, 0, 0, nth, 0), r3.face_phi(pv, D, nph-1, jm, nth, 1), r3.face_th(pv, D, nth-1, 0, nph, 1), r3.face_th(pv, D, 0, 0, km, 0), r3.face_phi(pv, D, km-1, 0, jm, 1), r3.face_th(pv, D, jm, km, nph, 0)]
        caps = [(0, nph, jm, nth), (0, km, 0, jm)]; inn = (km, nph, 0, jm)
    else:
        fs = [r3.face_phi(pv, D, nph-1, 0, nth, 1), r3.face_th(pv, D, 0, 0, nph, 0)]; caps = [(0, nph, 0, nth)]; inn = None
    fa = [pl.add_mesh(f, scalars='f', cmap=FUS, clim=(-L['vr'], L['vr']), show_scalar_bar=False, ambient=0.75, diffuse=0.3, specular=0.0) for f in fs if f.n_cells]
    ca = [pl.add_mesh(r3.cap(pv, D, cfg, *q), scalars='c', cmap=ICE, clim=L['te'], show_scalar_bar=False, ambient=c.get('amb', 0.3), diffuse=c.get('dif', 0.7), specular=c.get('spec', 0.15), specular_power=20) for q in caps]
    if inn: pl.add_mesh(r3.inner(pv, D, *inn), color='#555555', ambient=0.3)
    sb = dict(color='white', title_font_size=34, label_font_size=30, n_labels=3, fmt='%.3g', unconstrained_font_size=True)
    if not os.environ.get('HERO'): pl.add_scalar_bar('radial velocity (km/s)', mapper=fa[0].mapper, vertical=False, position_x=0.30, position_y=0.03, width=0.36, height=0.08, **sb)
    if not os.environ.get('HERO'): pl.add_scalar_bar('T_eff (kK)', mapper=ca[0].mapper, vertical=True, position_x=0.72, position_y=0.22, width=0.05, height=0.62, **sb)
    
    if c.get('note') and not os.environ.get('HERO'): pl.add_text(c['note'], position=(40, 40), font_size=13, color='white')
    if os.environ.get('HERO'): pl.add_text('t = %.1f d' % (D['t']/DAY), position=(36, 1080-70), font_size=20, color='white')
    tt = None if os.environ.get('HERO') else pl.add_text('%s\nt = %.2f d     median R_ph = %s R_sun' % (c['title'], D['t']/DAY, ('%.2f' if c['rview'] < 20 else '%.0f') % np.median(D['rph']/RS)), position=(960, 1040), font_size=15, color='white')
    if tt is not None: tt.prop.SetJustificationToCentered(); tt.prop.SetVerticalJustificationToTop()
    phc = 0.5*(D['phf'][0]+D['phf'][-1]); d0 = np.array([np.cos(phc), np.sin(phc), 0.]); up = np.array([0., 0., 1.]); side = np.cross(up, d0); az, el = np.radians(c['az']), np.radians(c['el']); eye = (d0*np.cos(az)+side*np.sin(az))*np.cos(el)+up*np.sin(el)
    TH, PH = np.meshgrid(np.linspace(D['thf'][0], D['thf'][-1], 9), np.linspace(D['phf'][0], D['phf'][-1], 9)); box = np.concatenate([r3.sph(np.full(TH.shape, rr), TH, PH).reshape(-1, 3) for rr in (D['rf'][0]/RS, c['rview'])]); foc = 0.5*(box.min(0)+box.max(0)); ext = np.linalg.norm(box.max(0)-box.min(0))
    pl.camera.position = tuple(foc+eye*ext*2); pl.camera.focal_point = tuple(foc); pl.camera.up = (0, 0, 1); pl.camera.view_angle = 30
    
    lead = None
    if c.get('fit'):                                           # wedge fills the frame: perspective projection of the wedge points -> view angle and window centre
        ep = foc+eye*ext*3.0; fwd = -eye; rt = np.cross(fwd, up); rt /= np.linalg.norm(rt); uu = np.cross(rt, fwd); q = box-ep; dep = q@fwd; X, Y = (q@rt)/dep, (q@uu)/dep; asp = SIZE[0]/SIZE[1]
        th = max((Y.max()-Y.min())/2/c.get('fh', 0.72), (X.max()-X.min())/2/(c.get('fw', 0.62)*asp)); va = 2*np.degrees(np.arctan(th)); cx, cy = (X.max()+X.min())/2/(th*asp), (Y.max()+Y.min())/2/th
        pl.camera.position = tuple(ep); pl.camera.focal_point = tuple(ep+fwd); pl.camera.up = (0, 0, 1); pl.camera.view_angle = va; pl.reset_camera_clipping_range(); pl.camera.SetWindowCenter(cx-c.get('tx', -0.06), cy-c.get('ty', 0.03)); lead = 1
    if lead is None: pl.reset_camera(bounds=(box[:, 0].min(), box[:, 0].max(), box[:, 1].min(), box[:, 1].max(), box[:, 2].min(), box[:, 2].max())); pl.camera.zoom(c['zoom']); pl.camera.SetWindowCenter(-0.085 if os.environ.get('HERO') else 0.0, 0.06)
    ln, dep = addbar(pv, pl, D, c, caps)
    if os.environ.get('CHK'):                                  # check data: projected outline of the true cap mesh (pixel box) and the bar length
        Q = []
        for q in caps:
            g = r3.cap(pv, D, cfg, *q); Q.append(np.asarray(g.points))
        Q = np.concatenate(Q); XY, _ = project(pl.camera, Q)
        print('CHK', c['title'], 'bar %.4g R_sun = %.2f px predicted; cap box px x %.1f %.1f y %.1f %.1f' % (c['bar'], ln, XY[:, 0].min(), XY[:, 0].max(), XY[:, 1].min(), XY[:, 1].max()), flush=True)
    pl.screenshot(out); pl.close()
def limits(s, D):
    te = r3.teff_map(D, dict(teff_smooth=0))/1e3; fv = np.concatenate([np.abs(D[k][..., :]).ravel() for k in D if k[:3] in ('fk_', 'fj_')])
    vis = []
    for k in D:
        if k[:3] in ('fk_', 'fj_'):
            i = D['idx'][int(k[3:])] if k[:2] == 'fk' else D['idx'][:, int(k[3:])]
            m = np.arange(D[k].shape[1])[None, :] <= i[:, None]; vis.append(np.abs(D[k][m]))
    return dict(te=[float(x) for x in np.percentile(te, [1, 99])], vr=float(np.percentile(np.concatenate(vis), 99)), Fle0=float(np.mean(D['Fph'] <= 0)))
def work(a):
    s, n, L, v, tag = a; t0 = time.time(); D = reduce(s, n); t1 = time.time(); f = D0 + ('/hero/' if os.environ.get('HERO') else '/frames/') + '%s%s_%05d.png' % (s, tag, n); render(D, RUNS[s], L, f, v); t2 = time.time(); print('frame %s %05d t=%.2f d reduce %.1f s render %.1f s R_ph med %.2f F<=0 %.1e' % (s, n, D['t']/DAY, t1-t0, t2-t1, np.median(D['rph'])/RS, np.mean(D['Fph'] <= 0)), flush=True)
def redw(q):
    t0 = time.time(); D = reduce(*q); print('red', q, '%.1f s' % (time.time()-t0), 'R_ph med %.2f p99.9 %.2f' % (np.median(D['rph'])/RS, np.percentile(D['rph'], 99.9)/RS), flush=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('star'); ap.add_argument('--only', default='last'); ap.add_argument('--workers', type=int, default=1); ap.add_argument('--reduce', action='store_true'); ap.add_argument('--variant', default='notch'); ap.add_argument('--tag', default=''); a = ap.parse_args()
    s = a.star; lf = D0 + '/lim_%s.json' % s; (x, y) = (lastdump(s),)*2 if a.only == 'last' else ((0, lastdump(s)) if a.only == 'all' else [int(q) for q in a.only.split('-')])
    if a.reduce:
        import multiprocessing as mp
        with mp.get_context('spawn').Pool(a.workers) as p: p.map(redw, [(s, n) for n in range(x, y+1)], chunksize=1)
        sys.exit()
    if not os.path.exists(lf): _l = limits(s, trim(reduce(s, lastdump(s)), RUNS[s]['rin']) if s == 'm35zams' else reduce(s, lastdump(s))); _l['te'] = RUNS[s]['te'] or _l['te']; json.dump(_l, open(lf, 'w'))
    L = json.load(open(lf)); print('limits', L); jobs = [(s, n, L, a.variant, a.tag) for n in range(x, y+1)]
    if a.workers == 1: [work(j) for j in jobs]
    else:
        import multiprocessing as mp
        with mp.get_context('spawn').Pool(a.workers) as p: p.map(work, jobs, chunksize=1)
