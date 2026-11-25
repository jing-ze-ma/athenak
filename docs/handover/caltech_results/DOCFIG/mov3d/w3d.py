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
RUNS = dict(bsg=dict(title='Blue supergiant', hydro=G + '/prod_c1_1124/bsg/bin/bsg3d.hydro_w.%05d.bin', rad='vet', var='vet_f1', grid=G + '/prod_c1_1124/bsg/bsg3d.x1grid.txt',
    table=G + '/fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt', mu=0.62, gamma=5/3., rin=0.0, rview=56.0, az=30, el=30, zoom=1.17))
SIZE = (1920, 1080)
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
    ca = [pl.add_mesh(r3.cap(pv, D, cfg, *q), scalars='c', cmap=ICE, clim=L['te'], show_scalar_bar=False, ambient=0.3, diffuse=0.7, specular=0.15, specular_power=20) for q in caps]
    if inn: pl.add_mesh(r3.inner(pv, D, *inn), color='#555555', ambient=0.3)
    sb = dict(color='white', title_font_size=34, label_font_size=30, n_labels=3, fmt='%.3g', unconstrained_font_size=True)
    pl.add_scalar_bar('radial velocity (km/s)', mapper=fa[0].mapper, vertical=False, position_x=0.30, position_y=0.03, width=0.36, height=0.08, **sb)
    pl.add_scalar_bar('T_eff (kK)', mapper=ca[0].mapper, vertical=True, position_x=0.72, position_y=0.22, width=0.05, height=0.62, **sb)
    pl.add_text('%s\nt = %.2f d     median R_ph = %.0f R_sun' % (c['title'], D['t']/DAY, np.median(D['rph']/RS)), position='upper_edge', font_size=15, color='white')
    phc = 0.5*(D['phf'][0]+D['phf'][-1]); d0 = np.array([np.cos(phc), np.sin(phc), 0.]); up = np.array([0., 0., 1.]); side = np.cross(up, d0); az, el = np.radians(c['az']), np.radians(c['el']); eye = (d0*np.cos(az)+side*np.sin(az))*np.cos(el)+up*np.sin(el)
    TH, PH = np.meshgrid(np.linspace(D['thf'][0], D['thf'][-1], 9), np.linspace(D['phf'][0], D['phf'][-1], 9)); box = np.concatenate([r3.sph(np.full(TH.shape, rr), TH, PH).reshape(-1, 3) for rr in (D['rf'][0]/RS, c['rview'])]); foc = 0.5*(box.min(0)+box.max(0)); ext = np.linalg.norm(box.max(0)-box.min(0))
    pl.camera.position = tuple(foc+eye*ext*2); pl.camera.focal_point = tuple(foc); pl.camera.up = (0, 0, 1); pl.camera.view_angle = 30
    pl.reset_camera(bounds=(box[:, 0].min(), box[:, 0].max(), box[:, 1].min(), box[:, 1].max(), box[:, 2].min(), box[:, 2].max())); pl.camera.zoom(c['zoom']); pl.camera.SetWindowCenter(0.0, 0.06); pl.screenshot(out); pl.close()
def limits(s, D):
    te = r3.teff_map(D, dict(teff_smooth=0))/1e3; fv = np.concatenate([np.abs(D[k][..., :]).ravel() for k in D if k[:3] in ('fk_', 'fj_')])
    vis = []
    for k in D:
        if k[:3] in ('fk_', 'fj_'):
            i = D['idx'][int(k[3:])] if k[:2] == 'fk' else D['idx'][:, int(k[3:])]
            m = np.arange(D[k].shape[1])[None, :] <= i[:, None]; vis.append(np.abs(D[k][m]))
    return dict(te=[float(x) for x in np.percentile(te, [1, 99])], vr=float(np.percentile(np.concatenate(vis), 99)), Fle0=float(np.mean(D['Fph'] <= 0)))
def work(a):
    s, n, L, v, tag = a; t0 = time.time(); D = reduce(s, n); t1 = time.time(); f = D0 + '/frames/%s%s_%05d.png' % (s, tag, n); render(D, RUNS[s], L, f, v); t2 = time.time(); print('frame %s %05d t=%.2f d reduce %.1f s render %.1f s R_ph med %.2f F<=0 %.1e' % (s, n, D['t']/DAY, t1-t0, t2-t1, np.median(D['rph'])/RS, np.mean(D['Fph'] <= 0)), flush=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('star'); ap.add_argument('--only', default='114-114'); ap.add_argument('--workers', type=int, default=1); ap.add_argument('--variant', default='notch'); ap.add_argument('--tag', default=''); a = ap.parse_args()
    s = a.star; lf = D0 + '/lim_%s.json' % s; x, y = (int(q) for q in a.only.split('-'))
    if not os.path.exists(lf): json.dump(limits(s, reduce(s, 114)), open(lf, 'w'))
    L = json.load(open(lf)); print('limits', L); jobs = [(s, n, L, a.variant, a.tag) for n in range(x, y+1)]
    if a.workers == 1: [work(j) for j in jobs]
    else:
        import multiprocessing as mp
        with mp.get_context('spawn').Pool(a.workers) as p: p.map(work, jobs, chunksize=1)
