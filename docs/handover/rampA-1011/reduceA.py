"""rampA (AG Car A, adapted from the AG Car B shake-down reducer): per-dump shell profiles + a few slices (run on Raven, nice).
usage: python3 reduce.py RUNDIR OUTDIR NPROC
Shell means are sin(theta)-weighted over the wedge.  eint = volumetric gas internal energy;
T, P_gas from the general-EOS table dump (inverted per log rho row, bilinear); kappa_R from the
TOPS table (linear, clamped).  Fluxes in units of L (full sphere, 4 pi r^2 <.>):
Lrad_lab = <F1>, Lrad_co = <F1 - v_r 4E/3> (P_rr ~ E/3), Lent_g = <(e+P) v_r>,
Lent_r = <4E/3 v_r>, Lkin = <rho v^2/2 v_r>, Lpot = <rho v_r> (-GM/r)."""
import glob
import os
import sys
import numpy as np
from multiprocessing import Pool
from scipy.ndimage import map_coordinates
from scipy.interpolate import RectBivariateSpline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc  # noqa: E402

D, OUT, NP = sys.argv[1], sys.argv[2], int(sys.argv[3])
HERE = os.path.dirname(os.path.abspath(__file__))
RS, MSUN, YR = 6.957e10, 1.98847e33, 3.15576e7
L0, RPH = 3.411709e39, 2.70140777e13
GM, C, AR, SIG = 7.29918420e27, 2.99792458e10, 7.5657332503e-15, 5.670374e-5
RFE = 245.2*RS

# ---- EOS inversion table: (log rho, log e/rho) -> log T, log p/rho
EOSD = os.environ.get('EOSDUMP', os.path.join(HERE, 'eos_dump.txt'))
with open(EOSD) as fh:
    hdr = [fh.readline() for _ in range(6)]
nx, ny, xmin, dx, ymin, dy = [float(v) for v in hdr[2].split()[1:]]
nx, ny = int(nx), int(ny)
d = np.loadtxt(EOSD, comments='#')
le = d[:, 0].reshape(ny, nx)   # [T, rho]
lp = d[:, 1].reshape(ny, nx)
lT = ymin + dy*np.arange(ny)
NE = 3000
le_lo, le_hi = le.min(), le.max()
leg = np.linspace(le_lo, le_hi, NE)
invT = np.empty((nx, NE))
invP = np.empty((nx, NE))
for i in range(nx):
    col = np.maximum.accumulate(le[:, i])
    invT[i] = np.interp(leg, col, lT)
    invP[i] = np.interp(leg, col, lp[:, i])

ross = os.environ['ROSS']   # rosseland_ext2_gs98_x0.36_z0.02.txt
nT, nD, lT0, dlT, lD0, dlD = 217, 421, 2.6, 0.025, -21.0, 0.05
kv = np.loadtxt(ross, comments='#').reshape(nT, nD)


def eos(rho, eint):
    lr = np.log10(np.maximum(rho, 1e-21))
    ler = np.log10(np.maximum(eint/rho, 1.0))
    ci = np.clip((lr - xmin)/dx, 0, nx - 1)
    cj = np.clip((ler - le_lo)/(leg[1] - leg[0]), 0, NE - 1)
    T = 10.0**map_coordinates(invT, [ci, cj], order=1)
    P = rho*10.0**map_coordinates(invP, [ci, cj], order=1)
    return T, P


def kap(rho, T):
    ci = np.clip((np.log10(T) - lT0)/dlT, 0, nT - 1)
    cj = np.clip((np.log10(np.maximum(rho, 1e-30)) - lD0)/dlD, 0, nD - 1)
    return 10.0**map_coordinates(kv, [ci, cj], order=1)


g = np.array([[float(v) for v in ln.split()] for ln in
              open(glob.glob(os.path.join(D, '*.x1grid.txt'))[0])
              if not ln.startswith('#') and len(ln.split()) == 4])
r = g[3:483, 2]
dr = g[3:483, 3]
SLICE_T = [1.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6, 7.0e6, 8.0e6]


def one(k):
    hf = sorted(glob.glob(os.path.join(D, 'bin', '*.hydro_w.*.bin')))[k]
    cf = os.path.join(OUT, 'cache', os.path.basename(hf) + '.npz')
    if os.path.exists(cf):
        z = np.load(cf)
        return {a: z[a] for a in z.files}
    h = bc.read_binary_as_athdf(hf)
    m = bc.read_binary_as_athdf(hf.replace('.hydro_w.', '.m1.'))
    th = h['x2v']
    w = np.sin(th)[None, :, None]
    w = w/np.sin(th).sum()/h['dens'].shape[0]
    rho = h['dens'].astype(np.float64)
    e = h['eint'].astype(np.float64)
    vr, vt, vp = (h[v].astype(np.float64) for v in ('velx', 'vely', 'velz'))
    E = m['m1_e'].astype(np.float64)
    F1 = m['m1_f1'].astype(np.float64)
    T, P = eos(rho, e)
    ka = kap(rho, T)
    v2 = vr**2 + vt**2 + vp**2

    def sm(a):
        return (a*w).sum(axis=(0, 1))
    A = 4*np.pi*r**2/L0
    o = {}
    o['rho'] = sm(rho)
    o['rho2'] = sm(rho**2)
    o['T'] = sm(T)
    o['Pg'] = sm(P)
    o['e'] = sm(e)
    o['E'] = sm(E)
    o['vr'] = sm(vr)
    o['vr2'] = sm(vr**2)
    o['vh2'] = sm(vt**2 + vp**2)
    o['rvr2'] = sm(rho*vr**2)
    o['cs2'] = sm((P + E/3.0)/rho)
    o['mflux'] = sm(rho*vr)
    o['Lrad_lab'] = A*sm(F1)
    Fco = F1 - vr*4.0*E/3.0
    o['Lrad_co'] = A*sm(Fco)
    o['Lent_g'] = A*sm((e + P)*vr)
    o['eP'] = sm(e + P)
    o['Lconv_g'] = o['Lent_g'] - A*o['eP']*o['vr']      # fluctuation (resolved convective) enthalpy flux
    o['Lconv_r'] = A*sm(4.0*E/3.0*vr) - A*4.0/3.0*o['E']*o['vr']
    o['ke'] = sm(0.5*rho*v2)
    o['ke_r'] = sm(0.5*rho*vr**2)
    o['Lent_r'] = A*sm(4.0*E/3.0*vr)
    o['Lkin'] = A*sm(0.5*rho*v2*vr)
    o['Lpot'] = A*o['mflux']*(-GM/r)
    o['kFco'] = sm(ka*Fco)
    o['kap_m'] = sm(ka)
    o['krho'] = sm(ka*rho)
    o['Tr'] = sm((E/AR)**0.25)
    ct = np.sqrt((5.0/3.0*P + 4.0/9.0*E)/rho)
    dth = th[1] - th[0]
    o['sig1'] = ((np.abs(vr) + ct)/dr).max(axis=(0, 1))
    o['sig23'] = ((np.maximum(np.abs(vt), np.abs(vp)) + ct)/(r*dth)).max(axis=(0, 1))
    o['vmax'] = np.sqrt(v2).max(axis=(0, 1))
    # per-column photosphere (tau 2/3 from the top along r)
    tau = np.cumsum((ka*rho*dr)[:, :, ::-1], axis=2)[:, :, ::-1]
    above = tau >= 2.0/3.0
    ii = above.shape[2] - 1 - np.argmax(above[:, :, ::-1], axis=2)
    rph_col = r[ii]/RS
    o['rph_col'] = np.array([rph_col.mean(), rph_col.min(), rph_col.max(), rph_col.std()])
    o['time'] = np.array([h['Time']])
    o['cycle'] = np.array([h['NumCycles']])
    t = h['Time']
    if any(abs(t - s) < 6e3 for s in SLICE_T) or k == 0:
        j = rho.shape[1]//2
        ife = int(np.argmin(abs(r - RFE)))
        i09 = int(np.argmin(abs(r - 0.9*RPH)))
        sl = {'vr': vr[:, j, :], 'drho': rho[:, j, :]/o['rho'] - 1, 'dT': T[:, j, :]/o['T'] - 1,
              'Fco': Fco[:, j, :]*A, 'rph_col': rph_col,
              'vr_fe': vr[:, :, ife], 'drho_fe': rho[:, :, ife]/o['rho'][ife] - 1,
              'vr_09': vr[:, :, i09], 'drho_09': rho[:, :, i09]/o['rho'][i09] - 1,
              'x2v': th, 'x3v': h['x3v'], 'r': r, 'time': t}
        np.savez_compressed(os.path.join(OUT, 'slice_%05d.npz' % k),
                            **{a: np.asarray(b, dtype=np.float32) for a, b in sl.items()})
    print('done', k, t, flush=True)
    np.savez(cf, **o)
    return o


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    fs = sorted(glob.glob(os.path.join(D, 'bin', '*.hydro_w.*.bin')))
    n = len([f for f in fs if os.path.exists(f.replace('.hydro_w.', '.m1.'))])
    if len(sys.argv) < 5 or sys.argv[4] != 'all':
        n -= 1   # the newest dump may still be written
    with Pool(NP) as pool:
        res = pool.map(one, range(n))
    keys = res[0].keys()
    np.savez_compressed(os.path.join(OUT, 'prof.npz'), r=r, dr=dr,
                        **{a: np.array([q[a] for q in res]) for a in keys})
    print('ALL DONE', n)
