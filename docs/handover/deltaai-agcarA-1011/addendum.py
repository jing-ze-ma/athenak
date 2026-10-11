# addendum.py <runs/JOBID>: does ic_rad relax perturb the envelope gas? (TASK-...-icrad-restart, item 5)
import sys, os, glob, warnings
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, '/u/jma20/ATHENAK/builds/wt/hes_gpu/vis/python')
import bin_convert as bc
S = sys.argv[1]; RPH = 388.3 * 6.957e10
BANDS = [('<0.8', 0, 0.8), ('0.95-1.05', 0.95, 1.05), ('1.05-1.5', 1.05, 1.5), ('>1.5', 1.5, 1e9)]
x1v = np.loadtxt('/work/nvme/bivj/jma20/agcarA_1011/agcar3d.x1grid.txt', max_rows=486)[3:483, 2]; r = x1v / RPH

def dumps(run):
    out = {}
    for f in glob.glob(f'{S}/{run}/bin/agcar3d.hydro_w.*.bin'):
        b = bc.read_binary(f); out.setdefault(b['cycle'], b)
    return out

def smean(b, v):
    num = 0.; den = 0.
    for m in range(b['n_mbs']):
        g = b['mb_geometry'][m]; a = np.array(b['mb_data'][v][m]); n2 = a.shape[1]
        w = np.sin(g[2] + (np.arange(n2) + .5) * (g[3] - g[2]) / n2)[None, :, None]
        num = num + (a * w).sum(axis=(0, 1)); den = den + w.sum() * a.shape[0]
    return num / den

def vmax(b, v='velx'):
    return max(np.abs(np.array(b['mb_data'][v][m])).max() for m in range(b['n_mbs']))

runs = ['B_def', 'B_relax', 'B_fsrelax']
D = {n: dumps(n) for n in runs}
print('cycle-0 dumps (written after the IC step); bitwise identical arrays vs B_def:')
for n in runs[1:]:
    for v in ('dens', 'velx', 'vely', 'velz', 'eint'):
        a = np.concatenate([np.ravel(x) for x in D[n][0]['mb_data'][v]]); b0 = np.concatenate([np.ravel(x) for x in D['B_def'][0]['mb_data'][v]])
        print(f'  {n:10s} {v:5s} identical={np.array_equal(a, b0)}  max|rel| cellwise={np.max(np.abs(a-b0)/np.maximum(np.abs(b0),1e-300)):.3e}')
print('\ncycle 0: max rel change of shell means vs B_def, per band')
for n in runs[1:]:
    for v in ('eint', 'dens'):
        m0 = smean(D['B_def'][0], v); m1 = smean(D[n][0], v); rel = np.abs(m1 - m0) / np.abs(m0)
        print(f'  {n:10s} {v:5s} ' + ' | '.join(f'{bn}: {rel[(r>=lo)&(r<hi)].max():.3e}' for bn, lo, hi in BANDS))
print('\nmax |velx| (cm/s) per dump')
for n in runs:
    print(f'  {n:10s} ' + '  '.join(f'c{c}: {vmax(D[n][c]):.4e}' for c in sorted(D[n])))
print('\nKE_int, KEr_int (hst cols 11, 12) per cycle')
H = {n: np.loadtxt(f'{S}/{n}/agcar3d.user.hst') for n in runs}
print('cycle ' + ' '.join(f'{n+" KE":>14s} {n+" KEr":>14s}' for n in runs))
for c in list(range(0, 11)) + [15, 20, 30, 40, 50, 60]:
    print(f'{c:5d} ' + ' '.join(f'{H[n][c,10]:14.6e} {H[n][c,11]:14.6e}' for n in runs))
for n in runs:
    k = H[n][:61, 10]; print(f'  {n:10s} KE_int max {k.max():.4e} at cycle {k.argmax()}, KEr_int max {H[n][:61,11].max():.4e} at cycle {H[n][:61,11].argmax()}')
print('\ncycle 0: signed shell-mean eint ratio relax/default at sample radii, and eint ratio at cycle 50, 60')
m0 = smean(D['B_def'][0], 'eint'); m1 = smean(D['B_relax'][0], 'eint')
for rr in (0.5, 0.8, 0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.2, 1.5, 2.0, 2.5, 2.9):
    i = np.argmin(np.abs(r - rr)); s = f'  r/Rph {r[i]:.3f}: c0 {m1[i]/m0[i]:.4f}'
    for c in (50, 60):
        if c in D['B_def'] and c in D['B_relax']:
            s += f'  c{c} {smean(D["B_relax"][c],"eint")[i]/smean(D["B_def"][c],"eint")[i]:.4f}'
    print(s)
