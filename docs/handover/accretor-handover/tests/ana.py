#!/usr/bin/env python3
"""Analysis of ry_per_accretor bin dumps.
usage: ana.py <bin file> [cs_kms] [mode]
  mode maxv   : max |v| and max |v|/cs over the domain (test a1)
  mode stream : stream ridge phi(r) vs the ballistic orbit (roche_stream.py integrator),
                min r of dense gas, impact phi and angle at the first active cell
Code units: Rsun, km/s."""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc  # noqa: E402

fn = sys.argv[1]
cs = float(sys.argv[2]) if len(sys.argv) > 2 else 15.5
mode = sys.argv[3] if len(sys.argv) > 3 else 'maxv'
fd = bc.read_binary(fn)
nmb = fd['n_mbs']
names = list(fd['mb_data'].keys())
# assemble r-phi arrays (one block in r; blocks split in phi), theta-averaged
blocks = []
for b in range(nmb):
    x1v = fd['mb_x1v'][b]
    x3v = fd['mb_x3v'][b]
    d = fd['mb_data']['dens'][b]
    vr = fd['mb_data']['velx'][b]
    vt = fd['mb_data']['vely'][b]
    vp = fd['mb_data']['velz'][b]
    blocks.append((x3v[0], x1v, x3v, d, vr, vt, vp))
blocks.sort(key=lambda t: t[0])
r = blocks[0][1]
ph = np.concatenate([b[2] for b in blocks])
D = np.concatenate([b[3] for b in blocks], axis=0)    # (k, j, i)
VR = np.concatenate([b[4] for b in blocks], axis=0)
VT = np.concatenate([b[5] for b in blocks], axis=0)
VP = np.concatenate([b[6] for b in blocks], axis=0)
print(f"file {fn} time {fd.get('time', float('nan'))} grid nr {len(r)} nphi {len(ph)} "
      f"ntheta {D.shape[1]}; vars {names}")
nan = int(np.sum(~np.isfinite(D)) + np.sum(~np.isfinite(VR)) + np.sum(~np.isfinite(VP)))
print(f"non-finite values: {nan}")
V = np.sqrt(VR**2 + VT**2 + VP**2)
if mode == 'maxv':
    k, j, i = np.unravel_index(np.argmax(V), V.shape)
    print(f"max|v| = {V.max():.4e} km/s = {V.max()/cs:.4e} cs at r {r[i]:.3f} phi "
          f"{np.degrees(ph[k]):.1f} deg (j {j}); rms|v| = {np.sqrt(np.mean(V**2)):.3e}; "
          f"max|v_theta| {np.abs(VT).max():.3e}")
    for lab, sl in (("inner 5 cells", slice(0, 5)), ("interior", slice(5, -5)),
                    ("outer 5 cells", slice(-5, None))):
        print(f"  {lab}: max|v| {V[:, :, sl].max():.4e} km/s = {V[:, :, sl].max()/cs:.3e} cs")
    print(f"  rho min/max {D.min():.4e} {D.max():.4e}")
    sys.exit(0)

# ---- stream mode
from scipy.integrate import solve_ivp  # noqa: E402
Ma, Md, a_rs = 6.24, 1.69, 30.3
mu = Md/(Ma + Md)
xa, xd = -mu, 1 - mu
G, Msun, Rsun = 6.674e-8, 1.989e33, 6.957e10
om = np.sqrt(G*(Ma + Md)*Msun/(a_rs*Rsun)**3)
aom = a_rs*Rsun*om/1e5
cs_d = np.sqrt(1.381e-16*6250/(1.27*1.6726e-24))/1e5
eps = cs_d/aom
gx = lambda x: (-(1 - mu)*(x - xa)/abs(x - xa)**3 - mu*(x - xd)/abs(x - xd)**3 + x)
from scipy.optimize import brentq  # noqa: E402
xL1 = brentq(gx, xa + 1e-3, xd - 1e-3)


def rhs(t, s):
    x, y, vx, vy = s
    r1 = np.hypot(x - xa, y)
    r2 = np.hypot(x - xd, y)
    return [vx, vy, -(1 - mu)*(x - xa)/r1**3 - mu*(x - xd)/r2**3 + x + 2*vy,
            -(1 - mu)*y/r1**3 - mu*y/r2**3 + y - 2*vx]


sol = solve_ivp(rhs, [0, 20], [xL1 - 1e-4, 0, -eps, 0], rtol=1e-11, atol=1e-13,
                max_step=1e-4)
bx, by = sol.y[0] - xa, sol.y[1]
br = np.hypot(bx, by)*a_rs
bph = np.arctan2(by, bx)
imin = np.argmax(np.diff(np.sign(np.diff(br))) > 0) + 1
print(f"ballistic: r_min {br[imin]:.4f} Rsun at phi {np.degrees(bph[imin]):.2f} deg")
rin = r[0]
# ridge phi(r) on the inbound branch, tracked inward from r_out
Dm = D.mean(axis=1)
print(" r(Rsun)  phi_ridge  phi_ball  dphi(deg)  rho_ridge")
nphi = len(ph)
kprev = None
for i in range(len(r) - 1, -1, -1):
    # ballistic phi at this r (inbound branch only)
    sel = np.where(br[:imin + 1] <= r[i])[0]
    if len(sel) == 0:
        continue
    pb = bph[sel[0]]
    if kprev is None:
        kc = int(np.argmin(np.abs(np.angle(np.exp(1j*(ph - pb))))))
    else:
        kc = kprev
    win = [(kc + dk) % nphi for dk in range(-nphi//60, nphi//60 + 1)]
    kk = win[int(np.argmax(Dm[win, i]))]
    kprev = kk
    if i % max(1, len(r)//16) == 0 or i == 0:
        dphi = np.degrees(np.angle(np.exp(1j*(ph[kk] - pb))))
        print(f" {r[i]:7.3f}  {np.degrees(ph[kk]):8.2f}  {np.degrees(pb):8.2f}  "
              f"{dphi:8.2f}  {Dm[kk, i]:.3e}")
# min radius of dense gas
for thr in (0.3, 0.1):
    ii = np.where((Dm > thr).any(axis=0))[0]
    if len(ii):
        print(f"min r with theta-mean rho > {thr}: {r[ii[0]]:.4f} Rsun (r_in {rin:.3f})")
# impact at the first active cell: mass flux peak
fl = -(Dm[:, 0]*VR.mean(axis=1)[:, 0])
k = int(np.argmax(fl))
vr0, vp0 = VR.mean(axis=1)[k, 0], VP.mean(axis=1)[k, 0]
sel = np.where(br[:imin + 1] <= r[0])[0]
if len(sel):
    print(f"ballistic at r = {r[0]:.4f}: phi {np.degrees(bph[sel[0]]):.2f} deg")
print(f"run impact (max inward mass flux at r {r[0]:.4f}): phi {np.degrees(ph[k]):.2f} deg, "
      f"v_r {vr0:.1f} v_phi {vp0:.1f} km/s, angle from inward normal "
      f"{np.degrees(np.arctan2(abs(vp0), -vr0)):.1f} deg; total inward flux share within "
      f"+-10 deg: {fl[np.abs(np.angle(np.exp(1j*(ph - ph[k])))) < np.radians(10)].sum()/max(fl[fl > 0].sum(), 1e-300):.3f}")
