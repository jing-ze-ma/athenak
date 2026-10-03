# flake8: noqa
"""Angle-averaged and per-ray profiles of pp2 dumps (64-70 ks) for the radial grid redesign.
-> prof70.npz.  H_p uses the TOTAL pressure p_gas + E/3 (as RESULTS.md sect. 3)."""
import bin_convert as bc
import ana3d as A3
import sys
import glob
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_mltpp_1002')
GM = 4.18143e26
R = A3.R
D = '/viper/ptmp2/jinma/he_mltpp_1002/race/pp2/'
rows = [
    l.split() for l in open(
        D +
        'hepresn.x1grid.txt') if l.strip() and l[0] != '#' and len(
            l.split()) >= 4]
xf = np.array([float(r[1]) for r in rows])
xv = np.array([float(r[2]) for r in rows])
IS, IE = 3, 472
rv = xv[IS:IE + 1]
re_ = np.append(xf[IS:IE + 1], xf[IE] + float(rows[IE][3]))
dr = np.diff(re_)


def one(fh):
    hw = bc.read_binary(fh)
    m1 = bc.read_binary(fh.replace('hydro_w', 'm1'))
    rho = A3.glob3(hw, 'dens')
    ei = A3.glob3(hw, 'eint')
    E = A3.glob3(m1, 'm1_e')
    v1 = A3.glob3(hw, 'velx')
    v2 = A3.glob3(hw, 'vely')
    v3 = A3.glob3(hw, 'velz')
    N3, N2, _ = rho.shape
    thf = np.linspace(hw['x2min'], hw['x2max'], N2 + 1)
    w = (np.cos(thf[:-1]) - np.cos(thf[1:]))[None, :, None] * np.ones((N3, 1, 1))
    w /= w.sum()
    T = A3.temperature(rho, ei)
    lt = np.log10(np.clip(T, 10**3.8, 10**6.5))
    lr = np.clip(np.log10(rho), -16, -2)
    pg = rho * 10**A3.HC.SLP.ev(lt, lr)
    kap = 10**A3.HC.SK.ev(lt, np.log10(rho))
    g = GM / rv**2
    P = pg + E / 3.0
    rm = (rho * w).sum(axis=(0, 1))
    Pm = (P * w).sum(axis=(0, 1))
    Ht = Pm / (rm * g)
    Hray = P / (rho * g)
    Hr5 = np.percentile(Hray.reshape(-1, rv.size), 5, axis=0)
    Hrmin = Hray.min(axis=(0, 1))
    vrm = (v1 * w).sum(axis=(0, 1))
    vr_rms = np.sqrt((((v1 - vrm)**2) * w).sum(axis=(0, 1)))
    vh_rms = np.sqrt(((v2**2 + v3**2) * w).sum(axis=(0, 1)))
    drho = np.sqrt((((rho - rm)**2) * w).sum(axis=(0, 1))) / rm
    rmax = rho.max(axis=(0, 1))
    rmin = rho.min(axis=(0, 1))
    ffl = ((rho < 3e-14) * w).sum(axis=(0, 1))
    dtau = (kap * rho * dr[None, None, :])[..., ::-1]
    tau = np.cumsum(dtau, axis=-1)
    rph = rv[::-1][np.argmax(tau > 2.0 / 3.0, axis=-1)]
    return dict(
        t=hw['time'],
        rm=rm,
        Ht=Ht,
        Hr5=Hr5,
        Hrmin=Hrmin,
        vr=vr_rms,
        vh=vh_rms,
        drho=drho,
        rmax=rmax,
        rmin=rmin,
        ffl=ffl,
        rph=rph,
        Em=(
            E *
            w).sum(
            axis=(
                0,
                1)),
        pgm=(
            pg *
            w).sum(
            axis=(
                0,
                1)))


if __name__ == '__main__':
    fs = sorted(glob.glob(D + 'bin/hepresn.hydro_w.*.bin'))[-12:]
    with Pool(6) as p:
        res = p.map(one, fs)
    out = {'rv': rv, 're': re_}
    for k in res[0]:
        out[k] = np.array([r[k] for r in res])
    np.savez('/viper/ptmp2/jinma/he_mltpp_1002/remap128/grid/prof70.npz', **out)
    print('t', out['t'])
