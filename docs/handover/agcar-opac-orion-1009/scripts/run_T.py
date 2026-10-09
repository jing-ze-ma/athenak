"""One temperature node: LTE Planck (absorption only) and Rosseland (absorption +
electron + Rayleigh scattering) means on the log rho grid.
usage: python run_T.py logT R [tag]      -> work/out/T<logT>_R<R>[_tag].npz
R = nu/dnu of the log-uniform opacity-sampling grid (300 .. 4e5 cm^-1).
Lines: Kurucz gfall (H..Zn, all stages present), Voigt (Doppler + 2 km/s micro-
turbulence, radiative damping), opacity sampling; the EXACT line Planck mean (sum of the
line strengths, no grid) is computed alongside."""
import glob, os, sys, time
import numpy as np
sys.path.insert(0, '/orion/ptmp/jinma/agcar_opac_1009/scripts')
import lowrho_lib as L

lT = float(sys.argv[1])
R = float(sys.argv[2])
tag = sys.argv[3] if len(sys.argv) > 3 else ''
WCAP = float(os.environ.get('WCAP', '0.01'))        # line wing cut, fraction of nu
EPS = float(os.environ.get('EPS', '1e-10'))         # skip a line in the grid below this
T = 10**lT
DLR = float(os.environ.get('DLR', '0.25'))
LR = np.round(np.arange(-21.0, -11.99, DLR), 2)
t0 = time.time()

at = L.Atoms()
UT = {k: at.U(*k, T) for k in at.ions}
# ionisation over the rho grid (molecular depletion applied below where FastChem says so)
chem = np.load(L.B + 'work/chem.npz')
ich = np.argmin(abs(chem['lT'] - lT)) if abs(chem['lT'] - lT).min() < 1e-6 else None
MOLS = {'C1O1': ('CO', {6: 1, 8: 1}, 28.010), 'C1N1': ('CN', {6: 1, 7: 1}, 26.018),
        'H1O1': ('OH', {1: 1, 8: 1}, 17.007), 'O1Si1': ('SiO', {14: 1, 8: 1}, 44.085),
        'H2O1': ('H2O', {1: 2, 8: 1}, 18.015), 'O1Ti1': ('TiO', {22: 1, 8: 1}, 63.866),
        'H2': (None, {1: 2}, 2.016)}


def npg_depleted(j):
    npg = dict(L.NPG)
    mol = {}
    if ich is None:
        return npg, mol
    for m, (nm, st, mw) in MOLS.items():
        jc = int(np.argmin(abs(chem['lr'] - LR[j])))
        assert abs(chem['lr'][jc] - LR[j]) < 1e-6
        nm_g = chem['n_' + m][ich, jc]
        mol[m] = nm_g
        for z, k in st.items():
            npg[z] -= k*nm_g
    return npg, mol


SAH = []
for j, lr in enumerate(LR):
    npg, mol = npg_depleted(j)
    ne, n, nHm = L.saha(at, T, 10**lr, UT, npg)
    SAH.append((ne, n, nHm, mol))
fmax = {k: max(SAH[j][1][k]/10**LR[j] for j in range(len(LR))) for k in at.ions}
ions = [k for k in at.ions if fmax[k] > 1e-12*L.NPG[k[0]]]

# frequency grid
nu0, nu1 = 300.0*L.c, 4.0e5*L.c
dln = 1.0/R
nnu = int(np.log(nu1/nu0)/dln) + 1
lnu0 = np.log(nu0)
nu = np.exp(lnu0 + dln*np.arange(nnu))
stim = -np.expm1(-L.h*nu/(L.kB*T))
Bn = L.planck_nu(nu, T)
dB = L.dBdT_nu(nu, T)
wP = Bn*nu*dln/(L.sigSB*T**4/np.pi)               # sum(wP) ~ 1
wR = dB*nu*dln
wR /= 4*L.sigSB*T**3/np.pi
print(f'logT {lT} R {R:.0e} nnu {nnu} ions {len(ions)} wP {wP.sum():.6f} '
      f'wR {wR.sum():.6f}', flush=True)

SIG = np.zeros((len(ions), nnu), np.float32)       # lines, cm^2 per particle of ion
SBF = np.zeros((len(ions), nnu), np.float32)       # bound-free x stim, per particle
PEX = np.zeros(len(ions))                          # exact line Planck sum per particle
PEXg = np.zeros(len(ions))                         # same, lines inside the grid only
nskip = 0
for ii, (z, s) in enumerate(ions):
    U = UT[(z, s)]
    sb = L.bf_ion(at, z, s, nu, T, U)*stim
    SBF[ii] = sb.astype(np.float32)
    ls = L.line_strengths(at, z, s, T, U)
    if ls is None:
        continue
    lnu, S, dnuD, a = ls
    Bl = L.planck_nu(lnu, T)/(L.sigSB*T**4/np.pi)
    PEX[ii] = np.sum(S*Bl)
    ing = (lnu > nu0) & (lnu < nu1)
    PEXg[ii] = np.sum((S*Bl)[ing])
    k0 = S/(np.sqrt(np.pi)*dnuD)*fmax[(z, s)]     # peak opacity per gram at max population
    use = ing & (k0 > EPS)
    nskip += np.sum(ing & ~use)
    vG = np.sqrt(np.log(np.maximum(k0/EPS, 1.0)))
    vL = np.sqrt(k0*a/(np.sqrt(np.pi)*EPS))
    cap = WCAP*lnu/dnuD
    nwin = np.minimum(np.maximum(np.maximum(vG, vL), 4.0), cap)
    sig = np.zeros(nnu)
    L.add_lines(sig, lnu0, dln, lnu[use], S[use], dnuD[use], a[use], nwin[use])
    SIG[ii] += sig.astype(np.float32)
print(f'lines/bf done {time.time()-t0:.0f}s, grid-skipped lines {nskip}', flush=True)

# molecules: DACE 1e-8 bar cross sections (cm^2/g of species), nearest T node <= 2900 K
# held above (OH to 4900 K)
MOLK = {}
if ich is not None:
    for m, (nm, st, mw) in MOLS.items():
        if nm is None:
            continue
        fs = glob.glob(f'/orion/ptmp/jinma/rsg_wind_1008/dace/raw/{nm}/Out_*_n800.bin')
        Ta = np.array(sorted(int(os.path.basename(f).split('_')[3]) for f in fs))
        Tu = Ta[np.argmin(abs(Ta - T))]
        f = glob.glob(f'/orion/ptmp/jinma/rsg_wind_1008/dace/raw/{nm}/Out_*_{Tu:05d}_n800.bin')[0]
        kd = np.fromfile(f, '<f4').astype(np.float64)
        wn = np.arange(len(kd))*0.01
        # exact Planck per gram of species on the native 0.01 cm^-1 grid
        nun = np.maximum(wn, 1e-3)*L.c
        pex = np.sum(kd*L.planck_nu(nun, T)*0.01*L.c)/(L.sigSB*T**4/np.pi)
        idx = np.round(nu/L.c/0.01).astype(np.int64)
        ks = np.where(idx < len(kd), kd[np.minimum(idx, len(kd) - 1)], 0.0)
        MOLK[m] = (ks*mw*L.amu, pex*mw*L.amu, Tu)       # per molecule
    print('molecule T nodes', {m: v[2] for m, v in MOLK.items()}, flush=True)

sH, sHe, sH2 = L.rayleigh(nu)
hbf = L.hminus_bf(nu)*stim
hff = L.hminus_ff(nu, T)
elem_of = np.array([z for z, s in ions])
out = {k: np.zeros(len(LR)) for k in
       ['kP', 'kP_num', 'kR', 'kR_abs', 'kes', 'ne', 'kP_lines', 'kP_bf', 'kP_ff',
        'kP_hm', 'kP_mol', 'kR_cont', 'kP_lines_outgrid']}
ZL = sorted(set(elem_of))
out['kP_lines_el'] = np.zeros((len(LR), len(ZL)))
for j, lr in enumerate(LR):
    rho = 10**lr
    ne, n, nHm, mol = SAH[j]
    w = np.array([n[k]/rho for k in ions])                   # per gram
    w32 = w.astype(np.float32)
    kbf = (w32 @ SBF).astype(np.float64)
    kabs = (w32 @ SIG).astype(np.float64) + kbf
    q2 = sum(s*s*n[(z, s)] for (z, s) in at.ions if s > 0)
    kff = L.ff_hydrogenic(nu, T, ne, q2)/rho
    khm = (nHm*hbf + n[(1, 0)]*ne*L.kB*T*hff)/rho
    kmol = np.zeros(nnu)
    pmol = 0.0
    for m, (ks, pex, Tu) in MOLK.items():
        kmol += mol[m]*ks
        pmol += mol[m]*pex
    kcont_abs = kff + khm
    kabs_tot = kabs + kcont_abs + kmol
    nH2 = mol.get('H2', 0.0)*rho
    kes = ne*L.sigT/rho
    ksca = kes + (n[(1, 0)]*sH + n[(2, 0)]*sHe + nH2*sH2)/rho
    # Planck: exact lines + numerical continuum (bf, ff, H-, molecules exact)
    lines_ex = np.sum(w*PEX)
    lines_grid = np.sum(w*PEXg)
    out['kP_num'][j] = np.sum(kabs_tot*wP) + (lines_ex - lines_grid)
    out['kP_lines'][j] = lines_ex
    out['kP_lines_outgrid'][j] = lines_ex - lines_grid
    out['kP_ff'][j] = np.sum(kff*wP)
    out['kP_hm'][j] = np.sum(khm*wP)
    out['kP_mol'][j] = pmol
    out['kes'][j] = kes
    out['ne'][j] = ne
    for iz, z in enumerate(ZL):
        out['kP_lines_el'][j, iz] = np.sum((w*PEX)[elem_of == z])
    ktot = kabs_tot + ksca
    out['kR'][j] = 1.0/np.sum(wR/ktot)
    out['kR_abs'][j] = 1.0/np.sum(wR/kabs_tot)
    out['kR_cont'][j] = 1.0/np.sum(wR/(kbf + kcont_abs + ksca))
    SAH[j] = (ne, n, nHm, mol, w)
PBF = (SBF.astype(np.float64) @ wP)                 # bf Planck per particle
for j in range(len(LR)):
    out['kP_bf'][j] = np.sum(SAH[j][4]*PBF)
out['kP'] = out['kP_lines'] + out['kP_bf'] + out['kP_ff'] + out['kP_hm'] + out['kP_mol']
os.makedirs(L.B + 'work/out', exist_ok=True)
fn = L.B + f'work/out/T{lT:.3f}_R{R:.0e}{tag}.npz'
np.savez(fn, lT=lT, lr=LR, R=R, ZL=np.array(ZL), wPsum=wP.sum(), wRsum=wR.sum(), **out)
print(f'wrote {fn} {time.time()-t0:.0f}s', flush=True)
for j in range(0, len(LR), 4):
    print(f'lr {LR[j]:6.2f} kP {out["kP"][j]:.3e} num {out["kP_num"][j]:.3e} '
          f'lines {out["kP_lines"][j]:.3e} bf {out["kP_bf"][j]:.2e} ff {out["kP_ff"][j]:.2e}'
          f' hm {out["kP_hm"][j]:.2e} mol {out["kP_mol"][j]:.2e} | kR {out["kR"][j]:.3e}'
          f' kes {out["kes"][j]:.3e} kRc {out["kR_cont"][j]:.3e}')
