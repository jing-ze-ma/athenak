"""e1div_prof.py RUNDIR OUTDIR [NPROC]  (DR, e1div-1014; python3 + numpy (+ scipy optional), bin_convert.py from ../e1caltech-1013)
Per hydro_w/m1 dump pair in RUNDIR/bin: radial profiles over the (theta, phi) shell -- sin(theta)-weighted mean, min, max -- of
rho, v_r, v_t (=sqrt(v_theta^2+v_phi^2)), eint, T_gas (general-EOS table eos_dump.txt.gz, nearest-row inversion), E, T_rad, F1;
plus count of cells with rho > 2e-16 (he_sponge_dmax) per radius.  -> OUTDIR/prof_<cycle>.npz (~0.2 MB each).
For dumps with t in SLICE_T (+- 6e4 s): the full (phi, theta) maps of rho, v_r, T, E, F1 at the radii IR, and the meridional
slice (all r, theta) at the phi index of the largest |v_r| above 1.0 R_ph -> OUTDIR/slice_<cycle>.npz (~10 MB each, float32)."""
import glob, gzip, os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e1caltech-1013'))
import bin_convert as bc  # noqa: E402
D, OUT = sys.argv[1], sys.argv[2]; NP = int(sys.argv[3]) if len(sys.argv) > 3 else 4
os.makedirs(OUT, exist_ok=True)
AR, RPH = 7.5657332503e-15, 2.70140777e13
SLICE_T = [9.0e6, 9.5e6, 9.9e6, 1.03e7, 1.08e7, 1.12e7, 1.16e7, 1.2e7, 1.22e7, 1.24e7]
IR = [341, 344, 347, 350, 360, 400, 450, 475, 490, 500, 508, 520]   # active radial indices (0..535)
with gzip.open(os.path.join(HERE, 'eos_dump.txt.gz'), 'rt') as fh:
    hdr = [fh.readline() for _ in range(6)]
    tab = np.loadtxt(fh, comments='#')
nx, ny, xmin, dx, ymin, dy = [float(v) for v in hdr[2].split()[1:]]; nx, ny = int(nx), int(ny)
le = tab[:, 0].reshape(ny, nx); lT = ymin + dy*np.arange(ny)
lemono = np.maximum.accumulate(le, axis=0)
def temp(rho, e):
    ci = np.clip(np.rint((np.log10(np.maximum(rho, 1e-21)) - xmin)/dx).astype(int), 0, nx - 1)
    ler = np.log10(np.maximum(e/rho, 1.0)); out = np.empty_like(rho)
    for c in np.unique(ci):
        w = ci == c; out[w] = 10**np.interp(ler[w], lemono[:, c], lT)
    return out
def one(hf):
    m1f = hf.replace('.hydro_w.', '.m1.')
    if not os.path.exists(m1f): return None
    h = bc.read_binary_as_athdf(hf); m = bc.read_binary_as_athdf(m1f)
    c = int(h['NumCycles']); t = float(h['Time'])
    po = os.path.join(OUT, 'prof_%06d.npz' % c)
    th = np.asarray(h['x2v']); w = np.sin(th)[None, :, None]; w = w/np.sin(th).sum()/h['dens'].shape[0]
    rho = np.asarray(h['dens'], float); e = np.asarray(h['eint'], float); vr = np.asarray(h['velx'], float)
    vt = np.sqrt(np.asarray(h['vely'], float)**2 + np.asarray(h['velz'], float)**2)
    E = np.asarray(m['m1_e'], float); F1 = np.asarray(m['m1_f1'], float); T = temp(rho, e); Tr = (np.maximum(E, 1e-30)/AR)**0.25
    o = dict(time=t, cycle=c, r=np.asarray(h['x1v'], float))
    for n, a in (('rho', rho), ('vr', vr), ('vt', vt), ('e', e), ('T', T), ('E', E), ('Tr', Tr), ('F1', F1)):
        o[n + '_mean'] = (a*w).sum(axis=(0, 1)); o[n + '_min'] = a.min(axis=(0, 1)); o[n + '_max'] = a.max(axis=(0, 1))
    o['n_above_dmax'] = (rho > 2e-16).sum(axis=(0, 1))
    np.savez(po, **o)
    if any(abs(t - s) < 6e4 for s in SLICE_T):
        x = o['r']/RPH; sel = x >= 1.0
        k = np.unravel_index(np.argmax(np.abs(vr[:, :, sel])), vr[:, :, sel].shape)[0]
        s = dict(time=t, cycle=c, r=o['r'], theta=th, ir=np.array(IR), kslice=k)
        for n, a in (('rho', rho), ('vr', vr), ('T', T), ('E', E), ('F1', F1), ('vt', vt)):
            s[n + '_at_ir'] = a[:, :, IR].astype(np.float32); s[n + '_merid'] = a[k].astype(np.float32)
        np.savez_compressed(os.path.join(OUT, 'slice_%06d.npz' % c), **s)
    return c, t
if __name__ == '__main__':
    fs = sorted(glob.glob(os.path.join(D, 'bin', '*.hydro_w.*.bin')))
    with Pool(NP) as p:
        for r in p.imap_unordered(one, fs):
            if r: print('done cycle %d t %.4e' % r, flush=True)
