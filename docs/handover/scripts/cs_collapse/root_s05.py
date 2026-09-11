"""si/s05 post-mortem root-cell analysis, following root2_task1..4 / twin2_analysis.
usage: root_s05.py <cyclediag_dir> <k> <j> <abort_cycle>
Finds the root cell (first cell above the deep interior to exceed 5000 K), then dumps
its 200-cycle history, the drain history at i and i+-4, and the per-cycle block counts."""
import sys
import os
import glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_fast_shim import read_rows, read_full
import numpy as np

DIR, K, J, CAB = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
fs = sorted(glob.glob(DIR+'/*.dat'))
cyc = [int(os.path.basename(f).split('.')[2]) for f in fs]
p = dict(zip(cyc, fs))
print('cyclediag files: %d, cycles %d .. %d' % (len(fs), min(cyc), max(cyc)))
print('column (k,j) = (%d,%d), abort cycle %d' % (K, J, CAB))

# ---- root cell: first cell with i > 45 to exceed 5000 K, earliest cycle ----
IMIN = 45
root_i, root_c = None, None
for c in sorted(cyc):
    r = read_rows(p[c], ['rt_T'], K, J)
    T = r['rt_T']
    hot = np.where(T[IMIN:130] > 5000.0)[0]
    if len(hot):
        root_i = IMIN + int(hot[0])
        root_c = c
        break
print('\n=== ROOT CELL: first cell at i > %d above 5000 K ===' % IMIN)
if root_i is None:
    print('  none in this column -- scanning the whole block instead')
    for c in sorted(cyc):
        r = read_full(p[c], ['rt_T'])
        T = r['rt_T'][:, :, IMIN:130]
        if T.max() > 5000.0:
            kk, jj, ii = np.unravel_index(np.argmax(T > 5000.0), T.shape)
            print('  cycle %d: (k,j,i) = (%d,%d,%d)  T = %.4e' %
                  (c, kk, jj, ii+IMIN, T[kk, jj, ii]))
            root_i, root_c = ii+IMIN, c
            K, J = int(kk), int(jj)
            break
else:
    print('  i = %d, first at cycle %d (%d cycles before the abort)' %
          (root_i, root_c, CAB-root_c))

I = root_i
V = ['u_ener', 'rt_T', 'w_dens', 'w_velx', 'w_eint', 'rt_src', 'rt_de', 'rt_Fb', 'rt_Ft',
     'rt_Em',
     'rt_clip', 'cond_f1', 'cond_dtcell', 'cond_kappa']
sel = sorted(set([c for c in cyc
                  if c <= root_c+20 and (root_c-c) % 5 == 0 and c >= root_c-200]
                 + [c for c in cyc if root_c-25 <= c <= min(root_c+40, max(cyc))]))
print('\n=== TASK 4a: root-cell history, (k,j,i) = (%d,%d,%d) ===' % (K, J, I))
print('%8s %9s %11s %11s %12s %12s %13s %12s %12s %12s %12s %12s %7s %12s %11s' % (
    'cycle', 'dt', 'T[K]', 'rho', 'v_r', 'e_int', 'delta_u', 'rt_src', 'rt_de', 'rt_Fb',
    'rt_Ft', 'rt_Em',
    'rt_clip', 'div(cond_f1)', 'cond_dtcell'))
for c in sel:
    r = read_rows(p[c], V, K, J)
    du = float('nan')
    if c-1 in p:
        du = r['u_ener'][I] - read_rows(p[c-1], ['u_ener'], K, J)['u_ener'][I]
    divf = r['cond_f1'][I+1] - r['cond_f1'][I]
    print('%8d %9.4f %11.4e %11.4e %12.4e %12.4e %13.5e %12.4e'
          ' %12.4e %12.4e %12.4e %12.4e %7.3f %12.4e %11.4e' % (
              c, r['dt'], r['rt_T'][I], r['w_dens'][I], r['w_velx'][I], r['w_eint'][I],
              du,
              r['rt_src'][I], r['rt_de'][I], r['rt_Fb'][I], r['rt_Ft'][I], r['rt_Em'][I],
              r['rt_clip'][I], divf, r['cond_dtcell'][I]))

print('\n=== TASK 4b: drain history at i-4, i, i+4 (rho, v_r) ===')
print('%8s %9s | %11s %12s | %11s %12s | %11s %12s' % (
    'cycle', 'dt', 'rho[i-4]', 'v_r[i-4]', 'rho[i]', 'v_r[i]', 'rho[i+4]', 'v_r[i+4]'))
for c in sel:
    r = read_rows(p[c], ['w_dens', 'w_velx'], K, J)
    print('%8d %9.4f | %11.4e %12.4e | %11.4e %12.4e | %11.4e %12.4e' % (
        c, r['dt'], r['w_dens'][I-4], r['w_velx'][I-4], r['w_dens'][I], r['w_velx'][I],
        r['w_dens'][I+4], r['w_velx'][I+4]))

print('\n=== TASK 4c: per-cycle counts over the whole meshblock ===')
print('active cells per block = 16*16*128 = %d' % (16*16*128))
print('%8s %13s %10s %12s %11s %11s %12s %12s' % (
    'cycle', 'time', 'n(T=200K)', 'n(rho<1e-9)', 'n(rt_clip)', 'n(T>1e4)', 'max T[K]',
    'min rho'))
step = max(1, (max(cyc)-min(cyc))//120)
selb = sorted(set([c for c in cyc if (c-min(cyc)) %
              step == 0] + [c for c in cyc if c >= root_c-30]))
for c in selb:
    r = read_full(p[c], ['rt_T', 'w_dens', 'rt_clip'])
    s = (slice(r['ks'], r['ke']+1), slice(r['js'], r['je']+1), slice(r['is'], r['ie']+1))
    T = r['rt_T'][s]
    rho = r['w_dens'][s]
    cl = r['rt_clip'][s]
    print('%8d %13.4f %10d %12d %11d %11d %12.4e %12.4e' % (
        c, r['time'], int(np.sum(np.abs(T-200.) < 0.5)), int(np.sum(rho < 1e-9)),
        int(np.sum(cl != 0)), int(np.sum(T > 1e4)), float(T.max()), float(rho.min())))
