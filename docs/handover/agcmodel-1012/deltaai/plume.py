"""plume.py RUNDIR OUT K1,K2,...: plume structure of the 3-D AG Car A dumps.
Per dump and radius: rms v_r, v_h, upflow area fraction f_up, <v_r>_up / <v_r>_down, enthalpy-flux split (up/down),
kinetic flux, corr(dT, v_r), dominant lateral wavelength of v_r (2-D FFT over theta-phi, power-weighted), in cells and in H_p;
mean-profile superadiabaticity nabla - nabla_ad (shell means, total P, StarPP.thermo).  Writes OUT/plume.txt and plume.npz."""
import sys, os, glob, types
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc
import make_ic_mlt_star as M
from make_ic_agcar_mltpp import StarPP
from multiprocessing import Pool
from scipy.ndimage import map_coordinates

D, OUT = sys.argv[1], sys.argv[2]
KS = [int(k) for k in sys.argv[3].split(',')]
RPH, RS, L0, GM = 2.70140777e13, 6.957e10, 3.411709e39, 7.29918420e27
AR = 7.5657332503e-15
os.makedirs(OUT, exist_ok=True)
# EOS inversion (as reduceA.py)
EF = os.environ.get('EOSDUMP', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eos_dump.txt'))
with open(EF) as fh:
    hdr = [fh.readline() for _ in range(6)]
nx, ny, xmin, dx, ymin, dy = [float(v) for v in hdr[2].split()[1:]]
nx, ny = int(nx), int(ny)
dd = np.loadtxt(EF, comments='#')
le = dd[:, 0].reshape(ny, nx); lp = dd[:, 1].reshape(ny, nx)
lT = ymin + dy*np.arange(ny)
NE = 3000
leg = np.linspace(le.min(), le.max(), NE)
invT = np.empty((nx, NE)); invP = np.empty((nx, NE))
for i in range(nx):
    col = np.maximum.accumulate(le[:, i])
    invT[i] = np.interp(leg, col, lT); invP[i] = np.interp(leg, col, lp[:, i])


def eos(rho, eint):
    lr = np.log10(np.maximum(rho, 1e-21)); ler = np.log10(np.maximum(eint/rho, 1.0))
    ci = np.clip((lr - xmin)/dx, 0, nx - 1); cj = np.clip((ler - leg[0])/(leg[1] - leg[0]), 0, NE - 1)
    return 10.0**map_coordinates(invT, [ci, cj], order=1), rho*10.0**map_coordinates(invP, [ci, cj], order=1)


a = types.SimpleNamespace(gm=GM, lum=L0, rstar=RPH, teff=9000.0,
                          eos='table:' + EF,
                          ross=os.environ['ROSS'],
                          ld_ext=-21.0, alpha=1.5, mltpp=0, f1=1e-4, f2=1e-3)
st = StarPP(a)
XS = [0.36, 0.40, 0.45, 0.50, 0.55, 0.63, 0.70, 0.80, 0.90, 0.95]


def one(k):
    hf = sorted(glob.glob(os.path.join(D, 'bin', '*.hydro_w.*.bin')))[k]
    h = bc.read_binary_as_athdf(hf)
    mf = hf.replace('.hydro_w.', '.m1.')
    m = bc.read_binary_as_athdf(mf) if os.path.exists(mf) else None
    r = h['x1v']; th = h['x2v']
    rho = h['dens'].astype(np.float64); e = h['eint'].astype(np.float64)
    vr = h['velx'].astype(np.float64); vt = h['vely'].astype(np.float64); vp = h['velz'].astype(np.float64)
    T, P = eos(rho, e)
    E = m['m1_e'].astype(np.float64) if m is not None else AR*T**4
    w = np.sin(th)[None, :, None]; w = w/w.sum()/rho.shape[0]
    sm = lambda q: (q*w).sum(axis=(0, 1))
    A = 4*np.pi*r**2/L0
    rows = []
    # mean-profile superadiabaticity
    Tm, Pm, rhom = sm(T), sm(P + E/3.0), sm(rho)
    nab = np.gradient(np.log(Tm), r)/np.gradient(np.log(Pm), r)
    nad = np.array([st.thermo(Tm[i], rhom[i])[6] for i in range(r.size)])
    Hp = Pm/(rhom*GM/r**2)
    dth = th[1] - th[0]
    for x in XS:
        i = int(np.argmin(abs(r/RPH - x)))
        v = vr[:, :, i]; dT = T[:, :, i]/sm(T)[i] - 1
        ww = w[:, :, 0]*rho.shape[0]
        ww = np.broadcast_to(ww, v.shape)/v.shape[0]
        up = v > 0
        fup = ww[up].sum()
        vup = (v*ww)[up].sum()/max(fup, 1e-30); vdn = (v*ww)[~up].sum()/max(1 - fup, 1e-30)
        hflux = (e[:, :, i] + P[:, :, i] + 4*E[:, :, i]/3)*v
        hm = (hflux*ww).sum(); vm = (v*ww).sum(); hbar = ((e[:, :, i] + P[:, :, i] + 4*E[:, :, i]/3)*ww).sum()
        Lc = A[i]*(hm - hbar*vm)
        Lc_up = A[i]*((hflux - hbar*v)*ww)[up].sum(); Lc_dn = A[i]*((hflux - hbar*v)*ww)[~up].sum()
        Lk = A[i]*(0.5*rho[:, :, i]*(vr[:, :, i]**2 + vt[:, :, i]**2 + vp[:, :, i]**2)*v*ww).sum()
        vrms = np.sqrt((v**2*ww).sum()); vh = np.sqrt(((vt[:, :, i]**2 + vp[:, :, i]**2)*ww).sum())
        c = np.corrcoef(dT.ravel(), v.ravel())[0, 1]
        Pk = np.abs(np.fft.fft2(v - v.mean()))**2
        k3 = np.fft.fftfreq(v.shape[0]); k2 = np.fft.fftfreq(v.shape[1])
        kk = np.sqrt(k3[:, None]**2 + k2[None, :]**2)
        kbar = (kk*Pk).sum()/Pk.sum()     # cycles per cell
        lam_cells = 1.0/max(kbar, 1e-30)
        lam_Hp = lam_cells*r[i]*dth/Hp[i]
        rows.append((x, vrms/1e5, vh/1e5, fup, vup/1e5, vdn/1e5, Lc, Lc_up, Lc_dn, Lk, c, lam_cells, lam_Hp,
                     nab[i] - nad[i], Hp[i]/(r[i]*dth)))
    return h['Time'], np.array(rows), r, nab - nad


if __name__ == '__main__':
    with Pool(min(len(KS), int(os.environ.get('NPROC', '8')))) as p:
        res = p.map(one, KS)
    with open(os.path.join(OUT, 'plume.txt'), 'w') as fh:
        hd = (' x    vrms  vh    f_up  v_up   v_dn   Lconv  Lc_up  Lc_dn  Lkin   cor(dT,vr) lam[cells] lam/Hp  (nab-nad)_mean Hp/dx_lat')
        for (t, rows, r, sad), k in zip(res, KS):
            fh.write('# dump %d t %.4e\n#%s\n' % (k, t, hd))
            for q in rows:
                fh.write('%5.2f %5.2f %5.2f %5.3f %6.2f %6.2f %6.3f %6.3f %6.3f %6.3f %6.3f %7.1f %6.2f %+8.4f %6.1f\n' % tuple(q))
    np.savez(os.path.join(OUT, 'plume.npz'), t=np.array([q[0] for q in res]), rows=np.array([q[1] for q in res]),
             r=res[0][2], sad=np.array([q[3] for q in res]))
    print(open(os.path.join(OUT, 'plume.txt')).read())
