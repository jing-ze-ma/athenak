# usage: binv.py <hydro_w bin file> <r_int> [gamma]: max |v| interior / top, per-cell CFL dt location
import sys, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/sprhd_0926/wt/vis/python')
import bin_convert as bc
f = bc.read_binary(sys.argv[1]); rint = float(sys.argv[2]); gam = float(sys.argv[3]) if len(sys.argv) > 3 else 5/3
out = []
for b in range(f['n_mbs']):
    g = f['mb_geometry'][b]
    nx1, nx2, nx3 = f['nx1_out_mb'], f['nx2_out_mb'], f['nx3_out_mb']
    r = g[0] + (np.arange(nx1) + 0.5)*(g[1] - g[0])/nx1
    th = g[2] + (np.arange(nx2) + 0.5)*(g[3] - g[2])/nx2
    dr = (g[1]-g[0])/nx1; dth = (g[3]-g[2])/nx2; dph = (g[5]-g[4])/nx3
    d = f['mb_data']['dens'][b]; e = f['mb_data']['eint'][b]
    v1 = f['mb_data']['velx'][b]; v2 = f['mb_data']['vely'][b]; v3 = f['mb_data']['velz'][b]
    cs = np.sqrt(gam*(gam-1)*e/d)
    R = r[None, None, :]; TH = th[None, :, None]
    dtc = np.minimum(np.minimum(dr/(abs(v1)+cs), R*dth/(abs(v2)+cs)), R*np.sin(TH)*dph/(abs(v3)+cs))
    out.append((np.broadcast_to(R, d.shape).ravel(), np.sqrt(v1**2+v2**2+v3**2).ravel(), cs.ravel(), dtc.ravel()))
R = np.concatenate([o[0] for o in out]); V = np.concatenate([o[1] for o in out])
C = np.concatenate([o[2] for o in out]); D = np.concatenate([o[3] for o in out])
mi = R <= rint
i1 = np.argmax(V*mi); i2 = np.argmax(V*(~mi)); i3 = np.argmin(D)
print('t=%.5e interior max|v|=%.3e (Mach %.3e at r=%.6e) | top max|v|=%.3e (Mach %.3e at r=%.6e) | '
      'min dt_cell/cfl=%.4e at r=%.6e (interior-only min %.4e)'
      % (f['time'], V[i1], V[i1]/C[i1], R[i1], V[i2], V[i2]/C[i2], R[i2], D[i3], R[i3], D[mi].min()))
