# usage: conv.py <run dir> [dump list]: horizontal-mean convection diagnostics per dump (sph wedge)
# v_r' rms (radial minus horizontal mean), v_lat rms, F_conv/F (gas enthalpy (gamma 5/3 approx) +
# 4/3 E_r fluctuation fluxes), at radii bands; v_MLT(He box) = 1.86e4 cm/s at r = R0 - 1.02e8.
import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/sprhd_0926/wt/vis/python')
import bin_convert as bc
d = sys.argv[1]; gam = 5/3; R0 = 3.240543e10
fs = sorted(glob.glob(d + '/bin/*.hydro_w.*.bin'))
sel = [int(s) for s in sys.argv[2].split(',')] if len(sys.argv) > 2 else range(0, len(fs), 4)
def glob3(f, key):
    nb = f['n_mbs']; geo = np.array(f['mb_geometry'])
    x1 = sorted(set(geo[:, 0])); x2 = sorted(set(geo[:, 2])); x3 = sorted(set(geo[:, 4]))
    n1, n2, n3 = f['nx1_out_mb'], f['nx2_out_mb'], f['nx3_out_mb']
    A = np.zeros((len(x3)*n3, len(x2)*n2, len(x1)*n1))
    for b in range(nb):
        g = geo[b]; i = x1.index(g[0]); j = x2.index(g[2]); k = x3.index(g[4])
        A[k*n3:(k+1)*n3, j*n2:(j+1)*n2, i*n1:(i+1)*n1] = f['mb_data'][key][b]
    r0, r1 = min(geo[:, 0]), max(geo[:, 1]); N = len(x1)*n1
    return A, r0 + (np.arange(N)+0.5)*(r1-r0)/N
zb = [-1.6e8, -1.3e8, -1.02e8, -0.7e8, -0.4e8, 0.0]
print('bands z [cm]:', zb, ' (R0 = %.6e)' % R0)
for s in sel:
    fh = bc.read_binary(fs[s]); fm = bc.read_binary(fs[s].replace('hydro_w', 'm1'))
    rho, r = glob3(fh, 'dens'); v1, _ = glob3(fh, 'velx'); v2, _ = glob3(fh, 'vely')
    v3, _ = glob3(fh, 'velz'); e, _ = glob3(fh, 'eint'); Er, _ = glob3(fm, 'm1_e'); F1, _ = glob3(fm, 'm1_f1')
    hm = lambda a: a.mean(axis=(0, 1))
    v1p = v1 - hm(rho*v1)/hm(rho)
    hg = gam*e; hr = 4/3*Er
    Fg = hm((hg - hm(hg))*v1p); Fr = hm((hr - hm(hr))*v1p)
    Ftot = hm(F1) + hm(hg*v1)
    vr = np.sqrt(hm(rho*v1p**2)/hm(rho)); vl = np.sqrt(hm(rho*(v2**2+v3**2))/hm(rho))
    vmr = np.sqrt(hm(v1)**2)
    sp = e/hm(e) - gam*rho/hm(rho); sp = sp - hm(sp)
    cs_ = hm(sp*v1p)/np.sqrt(hm(sp**2)*hm(v1p**2) + 1e-300)
    z = r - R0; out = []
    for z0 in zb:
        i = np.argmin(abs(z - z0))
        out.append('%4.2f|%4.2f|%8.1e|%8.1e|%5.2f' % (vr[i]/1.86e4, vl[i]/1.86e4, Fg[i]/Ftot[i], Fr[i]/Ftot[i], cs_[i]))
    i = np.argmax(Fg/Ftot)
    print('t=%8.0f  [vr1/vMLT|vlat/vMLT|Fg/F|Frc/F|corr(s1,vr1)] ' % fh['time'] + ' '.join(out) +
          '  max Fg/F %.2e at z=%.3e; <F1>/F0 mid %.4f' % (Fg[i]/Ftot[i], z[i], hm(F1)[len(r)//2]/2.475e15*(R0/r[len(r)//2])**2))
