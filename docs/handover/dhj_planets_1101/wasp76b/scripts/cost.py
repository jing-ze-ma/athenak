"""cost.py (dhj_planets_1101): cost per rotation of a planet deck relative to the WASP-121b
hydro fresh deck (w121hfresh_1101: C32 x nx1 112, measured dt 13.79 s, 0.0277 s/cycle,
3.7 min/rot on 1 viper MI300A node), and radiative / dynamical times on the IC column.

dt: the radial CFL estimate 0.3 min(dr/(c_s + 1 km/s)) on the IC column of each planet; the
planet's dt is taken as 13.79 s x (its estimate / the W121 estimate) (both runs are radial-
CFL limited near the jet base; the horizontal C32 limit is ~10x larger).
Cost per cycle scales with nx1 (same C32 angular mesh, same ck cadence).
t_rad(p) = c_p p / (4 sigma g T^3) (thermal relaxation of a layer of column mass p/g);
c_p = c_v Gamma_1/chi_rho from the 3x EOS dump (incl. H2 dissociation / ionisation).
usage: cost.py PROFILE GRAV AP FACES P_ROT[s] [NAME]
"""
import os
import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import planet_setup as ps  # noqa: E402
import dhj_remap as dr_  # noqa: E402
import grid112 as g112  # noqa: E402

BAR, SIG = 1e6, 5.670374419e-5
W_DT, W_SPC, W_MINROT, W_P = 13.79, 0.0277, 3.7, 110153.5


class CP:
    def __init__(self, fn):
        g = np.loadtxt(fn, comments='#')
        lr, lt = np.unique(g[:, 0]), np.unique(g[:, 1])
        rho, T = 10**g[:, 0], 10**g[:, 1]
        p, chr_, cht, cv = g[:, 3], g[:, 4], g[:, 5], g[:, 6]
        g1 = chr_ + p*cht**2/(rho*T*cv)
        cp = cv*g1/chr_
        self.f = RGI((lr, lt), np.log10(cp).reshape(len(lt), len(lr)).T, bounds_error=False,
                     fill_value=None)

    def __call__(self, rho, T):
        return 10**self.f(np.stack([np.log10(rho), np.log10(T)], -1))


def column(prof, grav, ap, faces, eos):
    col = ps.Column(prof, eos, ptop=1e-9)
    gm = grav*ap**2
    e = faces
    rc, dr, p, H = g112.cells(col, gm, ap, e)
    cs = np.interp(np.log10(p*BAR), col.lp, col.cs)
    dt = 0.3*np.min(dr/(cs + 1e5))
    return col, gm, rc, dr, p, dt


def trad(col, gm, ap, cpf, eos, pbar):
    T = col.T_at(pbar)
    rho = eos.rho(pbar*BAR, T)
    r = g112.r_of_p_ext(col, gm, ap, pbar)
    g = gm/r**2
    rho = float(np.ravel(rho)[0]); T = float(T)
    cp = float(np.ravel(cpf(np.array([rho]), np.array([T])))[0])
    return cp*pbar*BAR/(4*SIG*g*T**3), T, g


def main():
    prof, grav, ap, ffaces, prot = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), \
        sys.argv[4], float(sys.argv[5])
    name = sys.argv[6] if len(sys.argv) > 6 else 'planet'
    eos = ps.EOSTable(g112.EOSF)
    cpf = CP(os.path.join(os.path.dirname(g112.EOSF), 'eos_grid.txt'))
    fw = dr_.edges(g112.W_N, g112.W_AP, g112.W_X1MAX, np.array(g112.W_C, float))
    colw, gmw, _, _, _, dtw = column(g112.W_IC, g112.W_G, g112.W_AP, fw, eos)
    fp = np.loadtxt(ffaces)
    n = len(fp) - 1
    eosp = ps.EOSTable(g112.EOSF_P)            # planet EOS dump (env PLANET_EOSD)
    cpfp = CP(os.path.join(g112.EOSD_P, 'eos_grid.txt'))
    col, gm, _, _, _, dtp = column(prof, grav, ap, fp, eosp)
    dt = W_DT*dtp/dtw
    cyc_w, cyc = W_P/W_DT, prot/dt
    ratio = cyc*n/(cyc_w*112)
    print('radial dt estimate: W121 %.2f s, %s %.2f s -> expected dt %.2f s (W121 measured %.2f)'
          % (dtw, name, dtp, dt, W_DT))
    print('cycles/rot: W121 %.0f, %s %.0f;  nx1 %d (x %.3f per cycle)' % (cyc_w, name, cyc, n,
                                                                        n/112))
    print('cost per rotation relative to W121 hfresh: %.2f  -> %.1f min/rot on 1 viper node '
          '(W121 3.7)' % (ratio, ratio*W_MINROT))
    print('cost per unit simulated time relative to W121: %.2f' % (ratio*W_P/prot))
    print('300 rot: %.1f h on 1 viper node;  rotations for the W121 300-rot physical time '
          '(3.30e7 s): %.0f' % (300*ratio*W_MINROT/60, 300*W_P/prot))
    print('\n t_rad = c_p p/(4 sigma g T^3) on the IC columns [s] (and in rotations):')
    print('  p [bar]   W121: T      t_rad      rot  |  %s: T      t_rad      rot   ratio'
          % name)
    for pb in (1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100):
        tw, Tw, _ = trad(colw, gmw, g112.W_AP, cpf, eos, pb)
        tp, Tp, _ = trad(col, gm, ap, cpfp, eosp, pb)
        print('  %7.0e  %8.0f  %9.3e  %7.3f  | %8.0f  %9.3e  %7.3f  %.2f'
              % (pb, Tw, tw, tw/W_P, Tp, tp, tp/prot, tp/tw))


if __name__ == '__main__':
    main()
