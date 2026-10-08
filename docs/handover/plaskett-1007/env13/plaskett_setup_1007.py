#!/usr/bin/env python3
"""plaskett_setup_1007.py -- physical input numbers for the corrected Plaskett accretor setup
(AthenaK pgen ry_per_accretor, 2-D r-phi equatorial wedge, inner = envelope).

Code units (as in the pgen): length Rsun, velocity km/s, time Rsun/(km/s) = 6.957e5 s,
density in units of the stream's PEAK density at r_out (problem/rho_stream = 1).
Constants are the pgen's (kG, kMsun, kRsun, kKB, kMH).

Sections (numbering follows the task):
  1  sound speeds, t_don, Omega, Ryu+2025 A, B, C at L1, L1 stream widths
  2  3-D ballistic FAN from the L1 plane to the sphere r = r_out
  3  physical peak density at r_out for Mdot = 1e-4 Msun/yr -> env_rho_ph; Ryu Mdot check
  4  equilibrium photosphere r_ph(phi); envelope-top radii (psi = -15 c_ph^2)
  5  stream penetration depth into the n = 3 envelope -> env_wb_depth (wb_rmax)
  6  radial grid (polynomial + plateau stretch, src/coordinates/grid_stretch.hpp)
  7  ambient sanity (env_amb_rho, env_cs_amb, dfloor)
  8  measurement radius r_meas
  9  run-length keys
usage: python3 plaskett_setup_1007.py   [env NX1=500 to change the radial cell count]
"""
import os
import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.integrate import solve_ivp

np.set_printoptions(linewidth=120)
SEP = '-'*88

# ------------------------------------------------------------------ constants (pgen's)
G, MSUN, RSUN = 6.674e-8, 1.989e33, 6.957e10
KB, MH = 1.381e-16, 1.6726e-24
KVEL = 1.0e5
KTIME = RSUN/KVEL                       # s
GU = G*MSUN/(RSUN*KVEL**2)              # G Msun in (km/s)^2 Rsun
YR = 3.156e7

# ------------------------------------------------------------------ decided inputs
MA, MD = 16.0, 18.0                     # Msun (accretor / donor)
PDAY = 3.69                             # d
RACC_KEY = 9.00129260997162             # problem/r_acc (a face of the grid; see sec. 6)
R_IN = 6.3
MDOT = 1.0e-4                           # Msun/yr (Wade+2026)
T_GAIN, T_DON, MU = 28300.0, 26200.0, 0.62   # K (Wade+2026 MESA track, Case-A onset)
R_DON = 12.8                            # Rsun (= its Roche lobe)
KAPPA_E = 0.34
NPOLY = 3.0
SPIN = 1.0                              # spin 1: Phi_wb = Roche potential (Delta = 0)
# ambient / floor (current input, checked in sec. 7)
AMB_RHO, CS_AMB, DFLOOR = 1.0e-6, 300.0, 1.0e-7
NX1 = int(os.environ.get('NX1', 500))
DTH_WEDGE = 1.5769322 - 1.5646604       # theta extent of the 4-cell wedge (rad)

def cs_iso(T):
    """isothermal sound speed sqrt(kT/(mu m_H)) in km/s"""
    return np.sqrt(KB*T/(MU*MH))/KVEL

# ================================================================== 1
print(SEP); print('1. sound speeds, binary, L1 curvature (Ryu+2025 A, B, C), L1 widths')
c_ph = cs_iso(T_GAIN)
c_don = cs_iso(T_DON)
asep = (G*(MA + MD)*MSUN*(PDAY*86400.0)**2/(4*np.pi**2))**(1/3)/RSUN
om = np.sqrt(GU*(MA + MD)/asep**3)      # code units (km/s / Rsun)
porb = 2*np.pi/om
GMA, GMD = GU*MA, GU*MD
XCM = asep*MD/(MA + MD)                 # CM distance from the accretor (code's x_cm)
mu = MD/(MA + MD)
aom = asep*om                           # velocity unit of the ballistic integrator
q = MD/MA                               # Ryu's q = donor/accretor
print(f'c_ph (gainer Teff {T_GAIN:.0f} K, mu {MU}) = {c_ph:.4f} km/s   [problem/env_cs_ph]')
print(f'c_s,don (Teff {T_DON:.0f} K)            = {c_don:.4f} km/s   [problem/env_cs_stream;'
      f' problem/t_don = {T_DON:.1f} with mu_don {MU}]')
print(f'a (Kepler, {MA + MD:.0f} Msun, P {PDAY} d) = {asep:.4f} Rsun; Omega = {om:.6f} code '
      f'= {om/KTIME:.6e} /s; P_orb = {porb:.6f} code; a Omega = {aom:.3f} km/s')

def roche(r, phi):
    """the pgen's RochePot (equatorial plane, accretor-centred r, phi; donor at phi = 0)"""
    rd = np.sqrt(np.maximum(r*r + asep*asep - 2*asep*r*np.cos(phi), 1e-30))
    return -GMA/r - GMD/rd - 0.5*om**2*(r*r + XCM*XCM - 2*XCM*r*np.cos(phi))

def roche3(x, y, z):
    """3-D Roche potential, accretor at origin, donor at (asep, 0, 0)"""
    r1 = np.sqrt(x*x + y*y + z*z); r2 = np.sqrt((x - asep)**2 + y*y + z*z)
    return -GMA/r1 - GMD/r2 - 0.5*om**2*((x - XCM)**2 + y*y)

d_l1 = brentq(lambda x: (roche(x + 1e-6, 0) - roche(x - 1e-6, 0))/2e-6, 5, 30)
r_out = round(0.85*d_l1, 4)            # = mesh/x1max as written in the input (the pgen's r_out)
print(f'd_L1 = {d_l1:.5f} Rsun; r_out = 0.85 d_L1 = {0.85*d_l1:.6f} -> x1max key {r_out}')
# curvature at L1: Phi ~ (1/2) Omega^2 (A x^2 + B y^2 + C z^2)
S = GMA/d_l1**3 + GMD/(asep - d_l1)**3
A = -(2*S + om**2)/om**2
B = (S - om**2)/om**2
C = S/om**2
A_fit = -16.8 + 7.53*np.tanh(0.67*np.log10(q))**2
print(f'A = {A:.4f} (Ryu fit {A_fit:.3f}), B = {B:.4f} (= -(A+3)/2 = {-(A+3)/2:.4f}), '
      f'C = {C:.4f} (= B+1)')
sy_an = c_don/(om*np.sqrt(B))
sz_an = c_don/(om*np.sqrt(C))
FY, FZ = 0.932, 1.0                     # Ryu q=1 isothermal: FWHM 0.932x in y; z: none given
sy = FY*sy_an
sz = FZ*sz_an
print(f'L1 isothermal profile rho_L exp(-Omega^2 (B y^2 + C z^2)/(2 c^2)):')
print(f'  sigma_y = c/(Omega sqrt B) = {sy_an:.4f} Rsun analytic, x{FY} -> {sy:.4f} (Ryu-corrected)')
print(f'  sigma_z = c/(Omega sqrt C) = {sz_an:.4f} Rsun analytic; kept x{FZ} (Ryu: vertical HSE '
      f'holds, no z correction quoted) -> {sz:.4f}')
print(f'  old default stream_width c_s/Omega: {c_don/om:.4f} Rsun (with c_s,don {c_don:.2f})')

# ================================================================== 2
print(SEP); print('2. ballistic FAN (3-D restricted three-body, rotating frame, a = Omega = 1)')
xa, xd = -mu, 1 - mu

def gx(x):
    return -(1 - mu)*(x - xa)/abs(x - xa)**3 - mu*(x - xd)/abs(x - xd)**3 + x

xl1 = brentq(gx, xa + 1e-3, xd - 1e-3)
eps = c_don/aom

def rhs(s):
    x, y, z, vx, vy, vz = s
    r1 = np.sqrt((x - xa)**2 + y*y + z*z); r2 = np.sqrt((x - xd)**2 + y*y + z*z)
    f1, f2 = (1 - mu)/r1**3, mu/r2**3
    return np.array([vx, vy, vz,
                     -f1*(x - xa) - f2*(x - xd) + x + 2*vy,
                     -f1*y - f2*y + y - 2*vx,
                     -f1*z - f2*z])

# launch grid on the L1 plane: y in +-4.5 sigma_y, z in 0..4.5 sigma_z (mirror symmetric)
NY, NZ = 121, 31
yy = np.linspace(-4.5*sy, 4.5*sy, NY)
zz = np.linspace(0.0, 4.5*sz, NZ)
Y, Z = np.meshgrid(yy, zz, indexing='ij')
dy, dz = yy[1] - yy[0], zz[1] - zz[0]
# mass-flux weight rho_L(y,z) c dy dz (rho_L c = 1 here; normalised later); z<0 mirror -> x2
W = np.exp(-0.5*(Y/sy)**2 - 0.5*(Z/sz)**2)*dy*dz*np.where(Z > 0, 2.0, 1.0)
y0, z0, w0 = Y.ravel()/asep, Z.ravel()/asep, W.ravel()
n = y0.size
s = np.zeros((6, n))
s[0] = xl1 - 1.0e-4                    # the pgen's launch: x_L1 - 1e-4, v = -eps x_hat
s[1], s[2], s[3] = y0, z0, -eps
rc = r_out/asep
H = 2.0e-5
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fan_cache.npz')
key = np.array([NY, NZ, sy, sz, eps, rc, H])
if os.path.exists(CACHE) and np.allclose(np.load(CACHE)['key'], key, rtol=1e-14, atol=0):
    cz_ = np.load(CACHE); X_OUT, T_OUT, done, back = cz_['X'], cz_['T'], cz_['d'], cz_['b']
else:
    done = np.zeros(n, bool); back = np.zeros(n, bool)
    X_OUT = np.zeros((6, n)); T_OUT = np.zeros(n)
    rprev = np.sqrt((s[0] - xa)**2 + s[1]**2 + s[2]**2)
    t = 0.0
    while t < 6.0 and not np.all(done | back):
        sp = s.copy()
        k1 = rhs(s); k2 = rhs(s + 0.5*H*k1); k3 = rhs(s + 0.5*H*k2); k4 = rhs(s + H*k3)
        s = s + H*(k1 + 2*k2 + 2*k3 + k4)/6
        t += H
        r = np.sqrt((s[0] - xa)**2 + s[1]**2 + s[2]**2)
        hit = (~done) & (~back) & (rprev > rc) & (r <= rc)
        if hit.any():
            wl = (rprev[hit] - rc)/(rprev[hit] - r[hit])
            X_OUT[:, hit] = sp[:, hit] + wl*(s[:, hit] - sp[:, hit])
            T_OUT[hit] = t - H + wl*H
            done |= hit
        back |= (~done) & (s[0] > xl1 + 2e-3)  # fell back into the donor's lobe
        rprev = r
    np.savez(CACHE, key=key, X=X_OUT, T=T_OUT, d=done, b=back)
print(f'{n} trajectories (y x z = {NY} x {NZ}, z >= 0 mirrored); reached r_out: {done.sum()}, '
      f'fell back to the donor: {back.sum()}, neither (t > 6): {(~done & ~back).sum()}')
print(f'  flux fraction that reaches r_out: {w0[done].sum()/w0.sum():.6f}')
# accretor-centred state at r_out
dx, yo, zo = X_OUT[0] - xa, X_OUT[1], X_OUT[2]
vx, vy, vz = X_OUT[3], X_OUT[4], X_OUT[5]
Rc = np.hypot(dx, yo)
PHI = np.arctan2(yo, dx)
LAT = np.arctan2(zo, Rc)
r3 = np.sqrt(dx*dx + yo*yo + zo*zo)
VR = (vx*dx + vy*yo + vz*zo)/r3*aom
VP = (-vx*yo + vy*dx)/Rc*aom
VZ = vz*aom
m = done
w = w0[m]/w0[m].sum()
ic = np.argmin(np.abs(Y.ravel()) + np.abs(Z.ravel()))   # central trajectory y = z = 0
phi_c, vr_c, vp_c = PHI[ic], VR[ic], VP[ic]
phi_m = np.sum(w*PHI[m]); phi_sd = np.sqrt(np.sum(w*(PHI[m] - phi_m)**2))
lat_sd = np.sqrt(np.sum(w*LAT[m]**2))
vr_m, vp_m = np.sum(w*VR[m]), np.sum(w*VP[m])
vr_sd = np.sqrt(np.sum(w*(VR[m] - vr_m)**2)); vp_sd = np.sqrt(np.sum(w*(VP[m] - vp_m)**2))
print(f'central trajectory at r_out: phi {np.degrees(phi_c):.3f} deg, v_r {vr_c:.2f}, '
      f'v_phi {vp_c:.2f} km/s, |v| {np.hypot(vr_c, vp_c):.2f}, t {T_OUT[ic]/om*1:.4f} code')
print(f'flux-weighted at r_out:      phi {np.degrees(phi_m):.3f} deg (shift '
      f'{np.degrees(phi_m - phi_c):+.3f} deg), v_r {vr_m:.2f}, v_phi {vp_m:.2f} km/s')
print(f'  flux-weighted std: phi {np.degrees(phi_sd):.3f} deg = {phi_sd*r_out:.4f} Rsun; '
      f'v_r {vr_sd:.2f}, v_phi {vp_sd:.2f} km/s; latitude {np.degrees(lat_sd):.3f} deg = '
      f'{lat_sd*r_out:.4f} Rsun')
# Gaussian fit to the flux-weighted phi distribution (projected over z: the 2-D code's profile)
nb = 121
edges = np.linspace(phi_m - 5*phi_sd, phi_m + 5*phi_sd, nb + 1)
hist, _ = np.histogram(PHI[m], bins=edges, weights=w)
cen = 0.5*(edges[1:] + edges[:-1]); hist = hist/(edges[1] - edges[0])
gauss = lambda x, a0, x0, s0: a0*np.exp(-0.5*((x - x0)/s0)**2)
pfit, _ = curve_fit(gauss, cen, hist, p0=[hist.max(), phi_m, phi_sd])
sig_phi = abs(pfit[2])
# direct FWHM of the histogram (interpolated half-max crossings)
hm = 0.5*hist.max(); above = np.where(hist >= hm)[0]
il, ir = above[0], above[-1]
xl_ = np.interp(hm, [hist[il-1], hist[il]], [cen[il-1], cen[il]])
xr_ = np.interp(hm, [hist[ir+1], hist[ir]], [cen[ir+1], cen[ir]])
fwhm_dir = xr_ - xl_
print(f'Gaussian fit to the flux-weighted phi distribution: centre {np.degrees(pfit[1]):.3f} deg, '
      f'sigma_phi {sig_phi:.5f} rad = {np.degrees(sig_phi):.3f} deg; FWHM fit '
      f'{np.degrees(2.3548*sig_phi):.3f} deg, histogram {np.degrees(fwhm_dir):.3f} deg')
width_fan = sig_phi*r_out
sig_phi_fan = sig_phi
print(f'   fan: sigma_phi r_out = {width_fan:.4f} Rsun (Gaussian fit), {phi_sd*r_out:.4f} (std)')
# SHELL CROSSING: phi(y) on the midplane row must be monotonic for the fan to describe a
# stream; a fold means pressureless trajectories have crossed (caustic, infinite density)
P2 = PHI.reshape(NY, NZ); L2 = LAT.reshape(NY, NZ); D2 = done.reshape(NY, NZ)
yr = yy[D2[:, 0]]; pr = P2[D2[:, 0], 0]
dpr = np.diff(pr)
ext = np.where(np.sign(dpr[1:]) != np.sign(dpr[:-1]))[0] + 1   # folds of phi(y)
print('shell crossing in-plane: phi(y) on the midplane row has folds (caustics) at ' +
      ', '.join(f'y = {yr[k]/sy:+.2f} sigma_y (phi {np.degrees(pr[k]):.3f} deg)' for k in ext) +
      ('' if len(ext) else 'none') +
      f'; the core between them is single-valued; flux beyond the folds: '
      f'{w0[(Z.ravel() >= 0) & done & ((Y.ravel() < yr[ext[0]]) | (Y.ravel() > yr[ext[-1]]))].sum()/w0[done].sum():.3f}'
      if len(ext) else 'shell crossing in-plane: none')
zr = L2[NY//2, 1]/(zz[1]/r_out)
print(f'shell crossing vertical: z(r_out)/z(L1) on the central column = {zr:+.3f} '
      f'(< 0: the trajectories crossed the midplane; launch-plane frequencies sqrt(B) = '
      f'{np.sqrt(B):.2f}, sqrt(C) = {np.sqrt(C):.2f} Omega, flight time '
      f'{T_OUT[ic]:.3f}/Omega -> phases {np.sqrt(B)*T_OUT[ic]:.2f}, {np.sqrt(C)*T_OUT[ic]:.2f} rad '
      f'> pi/2)')
# PRESSURE-SUPPORTED (ADOPTED) widths at r_out: the isothermal stream stays in transverse
# hydrostatic balance (Ryu+2025: vertical HSE holds, in-plane the stream expands at 0.2-0.3 c,
# it does not focus), so sigma = c / sqrt(Phi_nn) with Phi_nn the potential curvature
# across the stream at the central trajectory's r_out crossing (as sigma_y, sigma_z at L1).
xo, yo_ = (X_OUT[0, ic] - xa)*asep, X_OUT[1, ic]*asep
vhat = np.array([X_OUT[3, ic], X_OUT[4, ic]]); vhat /= np.hypot(*vhat)
nhat = np.array([-vhat[1], vhat[0]])                    # in-plane normal to the stream
hh = 1e-3
def phi_nn(x0, y0, nv):
    f0 = roche3(x0, y0, 0.0)
    return (roche3(x0 + hh*nv[0], y0 + hh*nv[1], 0) - 2*f0
            + roche3(x0 - hh*nv[0], y0 - hh*nv[1], 0))/hh**2
Pnn = phi_nn(xo, yo_, nhat)
Pzz = (roche3(xo, yo_, hh) - 2*roche3(xo, yo_, 0) + roche3(xo, yo_, -hh))/hh**2
sperp = c_don/np.sqrt(Pnn)
sz_hse = c_don/np.sqrt(Pzz)
cosa = abs(vr_c)/np.hypot(vr_c, vp_c)
width = sperp/cosa                                      # arc width on the r_out sphere
sig_phi = width/r_out
print(f'pressure-supported widths at r_out (central crossing, stream angle from radial '
      f'{np.degrees(np.arccos(cosa)):.1f} deg):')
print(f'  in-plane: Phi_nn = {Pnn/om**2:.3f} Omega^2 -> sigma_perp = c_don/sqrt(Phi_nn) = '
      f'{sperp:.4f} Rsun; along the r_out arc sigma_perp/cos = {width:.4f} Rsun '
      f'(L1: {sy:.4f} corrected, {sy_an:.4f} analytic)')
print(f'  vertical: Phi_zz = {Pzz/om**2:.3f} Omega^2 -> sigma_z = {sz_hse:.4f} Rsun (L1 {sz:.4f})')
print(f'=> problem/stream_width = sigma_phi r_out = {width:.4f} Rsun ADOPTED (sigma_phi '
      f'{sig_phi:.5f} rad = {np.degrees(sig_phi):.3f} deg, FWHM {2.3548*width:.4f} Rsun); the '
      f'ballistic fan value {width_fan:.4f} is a caustic artefact')
# theta (z) width of the fan, for the record
hz, ez = np.histogram(np.concatenate([LAT[m], -LAT[m]]), bins=81,
                      weights=np.concatenate([w, w]))
czz = 0.5*(ez[1:] + ez[:-1])
pz, _ = curve_fit(gauss, czz, hz, p0=[hz.max(), 0.0, lat_sd + 1e-6])
sz_bal = abs(pz[2])*r_out
print(f'theta (z) width at r_out: ballistic fan sigma_z {sz_bal:.4f} Rsun (fit; std '
      f'{lat_sd*r_out:.4f}) vs pressure-supported {sz_hse:.4f}')
# velocity spread across the stream: within +-1 sigma (adopted) of the flux-weighted centre
inn = m & (np.abs(PHI - phi_m) < sig_phi)
wi = w0[inn]/w0[inn].sum()
print(f'velocity across the stream (|phi - <phi>| < 1 adopted sigma): v_r {VR[inn].min():.1f}..'
      f'{VR[inn].max():.1f}, v_phi {VP[inn].min():.1f}..{VP[inn].max():.1f} km/s; '
      f'all flux (2.5-97.5%): v_r {np.percentile(VR[m], 2.5):.1f}..{np.percentile(VR[m], 97.5):.1f},'
      f' v_phi {np.percentile(VP[m], 2.5):.1f}..{np.percentile(VP[m], 97.5):.1f}')
row = (Z.ravel() == 0) & m
core = row & (np.abs(Y.ravel()) <= sy)
print(f'  midplane launch |y| <= sigma_y: v_r {VR[core].min():.1f}..{VR[core].max():.1f}, v_phi '
      f'{VP[core].min():.1f}..{VP[core].max():.1f}, |v| {np.hypot(VR[core], VP[core]).min():.1f}..'
      f'{np.hypot(VR[core], VP[core]).max():.1f} km/s; v_z rms {np.sqrt(np.sum(w*VZ[m]**2)):.2f}')
IY0 = NY//2
dphidy = (P2[IY0+1, 0] - P2[IY0-1, 0])/(2*dy); dlatdz = L2[IY0, 1]/dz
Jc = abs(dphidy*dlatdz)
print(f'pressureless (Jacobian) area ratio r_out/L1 at the centre: '
      f'r_out^2 |d(phi,lat)/d(y,z)| = {r_out**2*Jc:.4f}')

# ================================================================== 3
print(SEP); print('3. physical density scale for Mdot = 1e-4 Msun/yr')
mdot = MDOT*MSUN/YR                     # g/s
vabs_c = np.hypot(vr_c, vp_c)
def rho_peak_rout(s_perp, s_z, v):
    """Mdot = rho_peak |v| 2 pi sigma_perp sigma_z  (= rho_peak |v_r| 2 pi sigma_arc sigma_z)"""
    return mdot/(v*KVEL*2*np.pi*(s_perp*RSUN)*(s_z*RSUN))
gacc = G*MA*MSUN/(9.0*RSUN)**2          # g at R_acc = 9.0 Rsun
P_ph = 2/3*gacc/KAPPA_E
rho_ph = P_ph/(c_ph*KVEL)**2
print(f'g = G M_acc/R_acc^2 = {gacc:.4e} cm/s^2; P_ph = (2/3) g/kappa_e = {P_ph:.4e}; '
      f'rho_ph = P_ph/c_ph^2 = {rho_ph:.4e} g/cc (+-2x: formula)')
print(f'Mdot = {mdot:.4e} g/s; central |v| at r_out {vabs_c:.2f}, v_r {vr_c:.2f} km/s')
cases = [('ADOPTED: pressure-supported sigma_perp x sigma_z at r_out', sperp, sz_hse),
         ('L1 widths kept (Ryu-corrected sigma_y, sigma_z)', sy, sz),
         ('L1 analytic widths', sy_an, sz_an),
         ('ballistic fan (caustic, upper bound on rho)', width_fan*cosa, sz_bal)]
for lab, a_, b_ in cases:
    rp = rho_peak_rout(a_, b_, vabs_c)
    print(f'  [{lab}] {a_:.4f} x {b_:.4f} Rsun: rho_peak(r_out) = {rp:.4e} g/cc -> '
          f'env_rho_ph = {rho_ph/rp:.4f}')
rho_pk = rho_peak_rout(sperp, sz_hse, vabs_c)
env_rho_ph = rho_ph/rho_pk
print(f'=> problem/env_rho_ph = rho_ph/rho_peak(r_out) = {env_rho_ph:.4f}   (viper 0.3; simple '
      f'estimate 0.051); rho_stream = 1 code <-> {rho_pk:.4e} g/cc')
# Ryu Mdot check: analytic 2 pi rho_L c sigma_y sigma_z, times the Coriolis correction
fR = 0.721 - 0.149*np.tanh(0.522*np.log10(q))**2
rhoL_an = mdot/(2*np.pi*c_don*KVEL*sy_an*sz_an*RSUN**2)
rhoL = rhoL_an/fR
print(f'Ryu Mdot/Mdot_analytic (isothermal, q = {q:.4f}) = {fR:.4f}')
print(f'rho_L implied by Mdot 1e-4: 2 pi rho_L c sy sz (analytic) = Mdot -> {rhoL_an:.4e}; '
      f'with the {fR:.3f} correction {rhoL:.4e} g/cc (sim peak 0.935x -> {0.935*rhoL:.4e})')
gdon = G*MD*MSUN/(R_DON*RSUN)**2
rho_ph_don = 2/3*gdon/KAPPA_E/(c_don*KVEL)**2
Hdon = (c_don*KVEL)**2/gdon
# isothermal Bernoulli/HSE: rho_L = rho_ph,don exp(-(Phi_L - Phi_ph)/c^2 - 1/2)
ovf = -(c_don*KVEL)**2*(np.log(rho_ph_don/rhoL) - 0.5)   # Phi_ph - Phi_L (> 0: overfill)
print(f'donor photosphere (same formula, R 12.8, Teff {T_DON:.0f}): g {gdon:.4e}, rho_ph,don = '
      f'{rho_ph_don:.4e} g/cc; rho_L/rho_ph,don = {rhoL/rho_ph_don:.1f} -> the donor photosphere '
      f'lies {ovf/(c_don*KVEL)**2:.2f} c^2 above Phi_L1, i.e. ~{ovf/gdon/RSUN:.3f} Rsun = '
      f'{ovf/gdon/Hdon:.1f} H_p (H_p,don {Hdon/RSUN:.4f} Rsun), relative overfill '
      f'{ovf/gdon/RSUN/R_DON:.1e}')
print(f'  i.e. at Mdot 1e-4 the L1 gas is {rhoL/rho_ph_don:.0f}x photospheric density: optically thick '
      f'(tau across sigma_y ~ kappa_e rho_L sigma_y = {KAPPA_E*rhoL*sy*RSUN:.0f}); isothermal at Teff '
      f'is then only an approximation')

# ================================================================== 4
print(SEP); print('4. equilibrium photosphere (code RochePot, spin 1) through (R_acc, 90 deg)')
phis = roche(RACC_KEY, np.pi/2)
def r_eq(phi, dpsi=0.0, lo=5.0):
    """r where Phi(r, phi) = Phi_s - dpsi (psi = Phi_s - Phi = dpsi)"""
    # upper bracket: inside the gainer's lobe (the equipotential is closed, r < d_L1)
    hi = d_l1 - 1e-3 if abs(np.cos(phi) - 1) < 1e-12 else 12.0
    return brentq(lambda r: roche(r, phi) - (phis - dpsi), lo, hi)
phd = np.arange(0, 181, 15)
rph = np.array([r_eq(np.radians(p)) for p in phd])
gph = np.array([(roche(r + 1e-5, np.radians(p)) - roche(r - 1e-5, np.radians(p)))/2e-5
                for r, p in zip(rph, phd)])
Hp = c_ph**2/gph
print(' phi[deg]  r_ph[Rsun]  g_r[code]  H_p=c_ph^2/g [Rsun]')
for p, r, g_, h in zip(phd, rph, gph, Hp):
    print(f'  {p:5.0f}   {r:9.5f}  {g_:9.1f}   {h:.5f}')
rph_max = rph.max()
pf = np.linspace(0, np.pi, 721); rf_ = np.array([r_eq(p) for p in pf])
rph_max = max(rph_max, rf_.max())
print(f'max r_ph = {rph_max:.5f} Rsun at phi {np.degrees(pf[np.argmax(rf_)]):.2f} deg '
      f'(R_acc {RACC_KEY:.5f}); H_p range {Hp.min():.5f}..{Hp.max():.5f} Rsun')
for p in (0, 90, 180):
    rt = r_eq(np.radians(p), -15*c_ph**2, lo=RACC_KEY - 0.5)
    print(f'  envelope top psi = -15 c_ph^2 at phi {p:3d}: r = {rt:.5f} Rsun')
# the pgen's default r_top: first r = R_acc exp(n dlg), dlg = ln(x1max/x1min)/nx1, where the
# c_ph isothermal atmosphere ABOVE (R_acc, phi = 90) has dropped by e^15 -- phi = 90 ONLY
dlg = np.log(r_out/R_IN)/NX1
nn = 1
while -(roche(RACC_KEY*np.exp(nn*dlg), np.pi/2) - phis)/c_ph**2 >= -15: nn += 1
rtop_def = RACC_KEY*np.exp(nn*dlg)
rtop0 = r_eq(0.0, -15*c_ph**2, lo=RACC_KEY - 0.5)
trunc = pf[rf_ > rtop_def]
print(f'pgen default r_top (phi 90 search) = {rtop_def:.5f} < r_ph(phi) for phi in '
      f'{np.degrees(trunc[trunc < np.pi/2]).min():.0f}..{np.degrees(trunc[trunc < np.pi/2]).max():.1f} '
      f'and {np.degrees(trunc[trunc > np.pi/2]).min():.1f}..{np.degrees(trunc[trunc > np.pi/2]).max():.0f}'
      f' deg (+ mirror): those cells (r > r_top) are initialised as AMBIENT, i.e. the photosphere '
      f'caps facing L1 and L2 are missing at t = 0 -> set problem/env_r_top >= {rtop0:.4f} '
      f'(envelope top psi = -15 c_ph^2 at phi 0)')

# ================================================================== 5
print(SEP); print('5. penetration depth at the impact point (central trajectory)')
def rhs3(t_, s_):
    return rhs(np.asarray(s_).reshape(6, 1)).ravel()
def ev_ph(t_, s_):
    x_, y_ = (s_[0] - xa)*asep, s_[1]*asep
    return roche(np.hypot(x_, y_), np.arctan2(y_, x_)) - phis
ev_ph.terminal = True; ev_ph.direction = -1  # Phi falls from Phi_L1 to Phi_s at the photosphere
s0 = [xl1 - 1e-4, 0, 0, -eps, 0, 0]
sol = solve_ivp(rhs3, [0, 6], s0, events=ev_ph, rtol=1e-11, atol=1e-13, dense_output=True)
if sol.t_events[0].size == 0:
    raise SystemExit('central trajectory never reaches the photosphere')
si = sol.y_events[0][0]
xi_, yi_ = (si[0] - xa)*asep, si[1]*asep
r_imp, phi_imp = np.hypot(xi_, yi_), np.arctan2(yi_, xi_)
v_imp = si[3:5]*aom
gradx = (roche3(xi_ + 1e-5, yi_, 0) - roche3(xi_ - 1e-5, yi_, 0))/2e-5
grady = (roche3(xi_, yi_ + 1e-5, 0) - roche3(xi_, yi_ - 1e-5, 0))/2e-5
nrm = np.array([gradx, grady])/np.hypot(gradx, grady)     # outward normal
vn = -v_imp@nrm
vimp = np.hypot(*v_imp)
vr_imp = (v_imp[0]*xi_ + v_imp[1]*yi_)/r_imp
print(f'impact on the photosphere: phi {np.degrees(phi_imp):.3f} deg, r {r_imp:.5f}; |v| '
      f'{vimp:.2f} km/s (rotating frame), v_n {vn:.2f} (normal to the equipotential), v_r '
      f'{vr_imp:.2f}; angle from the inward normal {np.degrees(np.arccos(vn/vimp)):.1f} deg')
vout = np.hypot(vr_c, vp_c)
P_env = lambda psi: env_rho_ph*c_ph**2*(1 + psi/((NPOLY + 1)*c_ph**2))**(NPOLY + 1)
res5 = {}
for lab, rho_i in [('2-D continuity rho = |v_out|/|v_imp|', vout/vimp),
                   ('conservative rho = rho_peak(r_out) = 1', 1.0)]:
    ram = rho_i*vn**2
    psi_p = brentq(lambda ps: P_env(ps) - ram, 0, 1e8)
    r_p = r_eq(phi_imp, psi_p, lo=R_IN - 0.5)
    rph_i = r_eq(phi_imp)
    d_imp = rph_i - r_p
    d_90 = RACC_KEY - r_eq(np.pi/2, psi_p, lo=R_IN - 0.5)
    res5[lab] = (rho_i, ram, psi_p, d_imp, d_90)
    print(f'  [{lab}] rho_imp {rho_i:.4f} (code), ram {ram:.4e}, ram/P_ph '
          f'{ram/P_env(0):.3e}; psi_pen {psi_p:.4e} = {psi_p/c_ph**2:.1f} c_ph^2; d_pen at '
          f'phi_imp {d_imp:.4f} Rsun; same psi at phi 90: depth {d_90:.4f}')
d_c = res5['2-D continuity rho = |v_out|/|v_imp|'][3]
d_k = res5['conservative rho = rho_peak(r_out) = 1'][3]
wbd = 1.5*d_c
print(f'=> env_wb_depth = 1.5 d_pen (continuity) = {wbd:.4f} Rsun -> hydro/wb_rmax = R_acc - '
      f'that = {RACC_KEY - wbd:.4f}; conservative (rho 1): 1.5 x {d_k:.4f} = {1.5*d_k:.4f} '
      f'-> {RACC_KEY - 1.5*d_k:.4f}')

# ================================================================== 6
print(SEP); print(f'6. radial grid: plateau radial stretch (mesh/use_grid_stretch_r_poly), nx1 = {NX1}')
# StretchRPoly (src/coordinates/grid_stretch.hpp): r = r0 + (r1 - r0) u(xi), xi = i/nx1 for
# face i, u = xi + sum c_k xi^k (1 - xi) + a [G(xi) - xi G(1)], G = int_0^xi g,
# g = (tanh((xi - xa)/w) - tanh((xi - xb)/w))/2.  Here c_k = 0 (PURE plateau): du/dxi is
# flat inside [xa, xb] (dr_in) and flat outside (dr_out) with tanh edges of half-width w.
# (The 10-07 grid fitted c1..c8 to a log target; with ~70% of the cells in the fine zone
# that fit leaves 15% wiggles inside the zone, so it is not used here.)
Hp_min = Hp.min()
drf = 0.25*Hp_min                       # dr <= H_p/4 with the smallest photospheric H_p
z_lo, z_hi = 8.80, rph_max + 0.15
WC, SHRINK, MARGIN = 12.0, 0.97, 0.004  # blend half-width (cells), dr_in = 0.97 drf, edge margin

def lcs(x):
    ax = np.abs(x); return ax + np.log1p(np.exp(-2*ax)) - np.log(2.0)
def Gp(xa_, xb_, w_, xi):
    return 0.5*w_*(lcs((xi - xa_)/w_) - lcs(-xa_/w_) - lcs((xi - xb_)/w_) + lcs(-xb_/w_))
def faces(nx1, p, R0=R_IN, R1=r_out):
    xi = np.linspace(0, 1, nx1 + 1)
    return R0 + (R1 - R0)*(xi + p[0]*(Gp(p[1], p[2], p[3], xi) - xi*Gp(p[1], p[2], p[3], 1.0)))
def build(nx1, w_, drin, xa_, xb_):
    """amplitude a so that du/dxi = 1 + a (1 - G(1)) = drin nx1/(r1 - r0) inside"""
    G1 = Gp(xa_, xb_, w_, 1.0)
    return [(drin*nx1/(r_out - R_IN) - 1)/(1 - G1), xa_, xb_, w_]
def zone(nx1, p, drmax):
    rf = faces(nx1, p); dr = np.diff(rf)
    i0 = i1 = int(np.argmin(dr))
    while i0 > 0 and dr[i0-1] <= drmax: i0 -= 1
    while i1 < nx1 - 1 and dr[i1+1] <= drmax: i1 += 1
    return rf[i0], rf[i1+1], rf, dr
def solve_grid(nx1):
    w_ = WC/nx1; drin = drf*SHRINK
    xa_, xb_ = 0.3, 0.9
    for _ in range(8):                  # zone edges (dr <= drf) at z_lo - margin, z_hi + margin
        xa_ = brentq(lambda x: zone(nx1, build(nx1, w_, drin, x, xb_), drf)[0] - (z_lo - MARGIN),
                     0.02, xb_ - 0.05)
        xb_ = brentq(lambda x: zone(nx1, build(nx1, w_, drin, xa_, x), drf)[1] - (z_hi + MARGIN),
                     xa_ + 0.05, 0.999)
    rf = faces(nx1, build(nx1, w_, drin, xa_, xb_))
    iR = int(np.argmin(np.abs(rf - RACC_KEY)))
    # shift the plateau rigidly (xa, xb + s) so face iR is EXACTLY r_acc: the pgen needs a face
    # within 1e-6 r_acc of problem/r_acc (inner = envelope fatal otherwise)
    sft = brentq(lambda s_: faces(nx1, build(nx1, w_, drin, xa_ + s_, xb_ + s_))[iR] - RACC_KEY,
                 -2.0/nx1, 2.0/nx1, xtol=1e-17)
    p = build(nx1, w_, drin, xa_ + sft, xb_ + sft)
    el, eh, rf, dr = zone(nx1, p, drf)
    ratio = np.maximum(dr[1:]/dr[:-1], dr[:-1]/dr[1:]).max()
    return p, rf, dr, el, eh, iR, ratio

print(f'H_p,min at the photosphere = {Hp_min:.5f} Rsun (phi ~75-90) -> dr target H_p/4 = '
      f'{drf:.4e} Rsun; fine zone wanted {z_lo:.2f}..{z_hi:.4f} (= max r_ph + 0.15)')
print(' nx1  dr_in      zone(dr<=H_p/4)    n_fine  max ratio  max dr/r   dr(r_in)   dr(R_acc)  '
      'dr(9.5)    dr(r_out)')
grids = {}
for nx in sorted({500, 640, 768, NX1}):
    p, rf, dr, el, eh, iR, ratio = solve_grid(nx)
    grids[nx] = (p, rf, dr, el, eh, iR, ratio)
    di = lambda x: dr[min(np.searchsorted(rf, x, side='right') - 1, nx - 1)]
    print(f' {nx:4d} {drf*SHRINK:.3e}  {el:.4f}..{eh:.4f}  {np.sum(dr <= drf):5d}   {ratio:.4f}'
          f'    {np.max(dr/rf[:-1]):.2e}  {dr[0]:.3e}  {di(RACC_KEY):.3e}  {di(9.5):.3e}  {dr[-1]:.3e}')
p, rf, dr, el, eh, iR, ratio = grids[NX1]
RF = rf
print(f'chosen nx1 = {NX1}: fold-free {np.all(dr > 0)}; face {iR} = {rf[iR]:.15g} '
      f'(r_acc key {RACC_KEY}); dr/H_p,min at R_acc {dr[iR]/Hp_min:.3f}')
for x in (R_IN, 8.5, 8.8, RACC_KEY, 9.5, 9.65, 10.0, 11.0, r_out):
    i = min(np.searchsorted(rf, x, side='right') - 1, NX1 - 1)
    print(f'  r {x:8.4f}: dr {dr[i]:.4e}, dr/r {dr[i]/rf[i]:.3e}, r dphi(2048) '
          f'{rf[i]*2*np.pi/2048:.4e}')
print('<mesh> keys (c1..c4 must be present = 0; c5..c8 may be omitted; <meshblock> nx1 = nx1):')
print(f'nx1       = {NX1}\nx1min     = {R_IN}\nx1max     = {r_out:.4f}')
print('use_grid_stretch_r_poly = true')
for k in range(4):
    print(f'f_stretch_r_c{k+1} = 0.0')
print(f'f_stretch_r_p_amp = {p[0]:.17g}\nf_stretch_r_p_xa = {p[1]:.17g}\n'
      f'f_stretch_r_p_xb = {p[2]:.17g}\nf_stretch_r_p_w = {p[3]:.17g}')
np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'faces_nx{NX1}.txt'),
           rf, fmt='%.15g')
# resolution of the envelope where WB is switched off (plain gravity below the fine zone)
for lab, dd in [('continuity', wbd), ('conservative', 1.5*d_k)]:
    rw = RACC_KEY - dd
    psi_w = phis - roche(rw, np.pi/2)
    gw = (roche(rw + 1e-5, np.pi/2) - roche(rw - 1e-5, np.pi/2))/2e-5
    Hw = (c_ph**2 + psi_w/(NPOLY + 1))/gw
    for nx in sorted(grids):
        rr_, dd_ = grids[nx][1], grids[nx][2]
        i = np.searchsorted(rr_, rw, side='right') - 1
        print(f'  WB switch-off at r = {rw:.4f} ({lab}): envelope H_p = {Hw:.4f} Rsun; nx1 {nx}: '
              f'dr {dd_[i]:.4e} -> H_p/dr {Hw/dd_[i]:.1f}')

# ================================================================== 7
print(SEP); print('7. ambient sanity (code units: density in rho_peak(r_out))')
P_ph_code = env_rho_ph*c_ph**2
ram_out = 1.0*(vr_c**2 + vp_c**2)
inflow = np.sqrt(2*np.pi)*sig_phi*r_out*abs(vr_c)*r_out*DTH_WEDGE   # per code time, wedge
rt90 = r_eq(np.pi/2, -15*c_ph**2, lo=RACC_KEY - 0.5)
print(f'P_ph = env_rho_ph c_ph^2 = {P_ph_code:.3f}; stream ram at r_out = 1 x |v|^2 = '
      f'{ram_out:.3e}; stream inflow (wedge) {inflow:.3f}/code time = {inflow*porb:.3f}/orbit')
print('reference P_amb/P_ph: RY Per h5 (PASS) 1.9e-5; RY Per a3 (env_amb_rho 1e-4, FAILED: outflow '
      'at phi ~ 0) 1.9e-3; Plaskett env12 (0.3 x 400 = 120) 7.5e-4')
for ar, fl in ((AMB_RHO, DFLOOR), (2e-7, 2e-8), (1e-7, 1e-8), (5e-9, 5e-10)):
    P_amb = ar*CS_AMB**2
    def ramb(r, ph):
        return max(ar*np.exp(-(roche(r, ph) - roche(RACC_KEY, ph))/CS_AMB**2), fl)
    amb_out = ramb(r_out, phi_c)
    rr = np.linspace(rt90, r_out, 1001); phg = np.linspace(0, 2*np.pi, 91)[:-1]
    Mamb = np.mean([np.trapezoid([ramb(r, ph)*r*r for r in rr], rr) for ph in phg])* \
        2*np.pi*DTH_WEDGE
    print(f' env_amb_rho {ar:.0e} (dfloor = rho_amb {fl:.0e}): P_amb(R_acc)/P_ph '
          f'{P_amb/P_ph_code:.2e} (cold layer ends {np.log(P_ph_code/P_amb):.1f} H_p up); '
          f'rho_amb(r_out) {amb_out:.2e} = {amb_out/fl:.1f} floor; P_amb(r_out)/ram '
          f'{amb_out*CS_AMB**2/ram_out:.1e}; ambient mass (r > {rt90:.3f}) {Mamb:.2e} = '
          f'{Mamb/(inflow*porb):.1e} of the inflow per orbit')

# ================================================================== 8
print(SEP); print('8. measurement radius')
Hp_max = Hp.max()
for lab, hh in [('H_p,max (phi 0 side)', Hp_max), ('H_p at phi 90', Hp[phd == 90][0])]:
    target = rph_max + 10*hh
    if 'RF' in globals():
        k = np.searchsorted(RF, target)
        print(f'  max r_ph + 10 {lab} = {rph_max:.4f} + 10 x {hh:.5f} = {target:.4f} -> '
              f'face {k}: r_meas = {RF[k]:.6f}; r_out - r_meas = {r_out - RF[k]:.4f} Rsun')

# ================================================================== 9
print(SEP); print('9. run length')
print(f'P_orb = 2 pi/Omega = {porb:.6f} code (= {porb*KTIME/86400:.4f} d)')
print(f'time/tlim = 10 P_orb = {10*porb:.5f}; output dt = P/20 = {porb/20:.6f}')
print(f'old input tlim 8.525 = {8.525/porb:.2f} orbits of THIS binary; bin dt 0.04263 = P/'
      f'{porb/0.04263:.2f}  (both RY Per numbers, P_orb 0.8519)')
