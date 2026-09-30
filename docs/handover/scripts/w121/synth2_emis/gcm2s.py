"""gcm2s.py: explain the emis11-vs-rt_Fb gap by re-solving the SAME 11-band columns with the
GCM's own spherical two-stream closure (two_stream_rt.hpp, ck_spherical: plane-parallel
cells, face mixing c = beta (d_above - u_below), beta = (A_above - A_below)/(A_above + A_below),
lagged down ray from a plane-parallel probe up-sweep, the up-sweep reusing the down-sweep's c;
2-point Gauss angles = ck_nquad 2).  Cells here use the linear-in-tau source of rtcore.

  python gcm2s.py w1x|w10x     -> prints L and per-band ratios vs the olr dump
"""
import sys
import numpy as np
import rtcore as rc
import ck11
from emis11 import CK

arm = sys.argv[1]
od = {'w1x': 'synth_rot300', 'w10x': 'synth_rot300_10x'}[arm]
st = np.load(rc.SD + 'state/%s.npz' % arm)
T, p, rho, rf = st['T'], st['p'], st['rho'], st['rf']
dOm = st['dOm']
ncol, nlay = T.shape
ck = ck11.CK(**CK[arm])
kr = ck.kr(T.ravel(), p.ravel()/1e6, rho.ravel()).reshape(ncol, nlay, 11, 8)
B = rc.BandPlanck(ck.wl)(T)[..., ::-1]
dz = np.diff(rf)
dtau = kr*dz[None, :, None, None]                       # (ncol, nlay, 11, 8)
Sf = np.moveaxis(rc.face_source(np.moveaxis(B[..., None], 1, 0), np.moveaxis(dtau, 1, 0)), 0, 1)
mu, w = rc.gauss01(2)
Af = rf**2                                              # face areas (up to a constant)
Ac = (rf[1:]**3 - rf[:-1]**3)/(3*dz)                    # cell frame areas V/dx
# beta at face f (0..nlay): below = cell f-1 (or the face itself at f=0), above = cell f
Ab = np.concatenate([[Af[0]], Ac])
Aa = np.concatenate([Ac, [Af[-1]]])
beta = (Aa - Ab)/(Aa + Ab)


def step(Iin, t, s_in, s_out):
    e = np.exp(-t)
    a = np.where(t > 1e-4, -np.expm1(-t)/np.maximum(t, 1e-30), 1 - 0.5*t)
    return Iin*e + s_out*(1 - a) + s_in*(a - e)


res = {}
for sph in (False, True):
    F = np.zeros((ncol, 11))
    for m, (mu_, w_) in enumerate(zip(mu, w)):
        t = dtau/mu_
        # probe up-sweep (plane-parallel)
        up = np.zeros((nlay + 1,) + dtau.shape[:1] + dtau.shape[2:])
        up[0] = Sf[:, 0]
        for l in range(nlay):
            up[l + 1] = step(up[l], t[:, l], Sf[:, l], Sf[:, l + 1])
        if sph:
            c = np.zeros_like(up)
            d = np.zeros_like(up[0])                   # d_above at the top face = 0
            for f in range(nlay, 0, -1):
                c[f] = beta[f]*(d - up[f])
                d = d + c[f]                           # d_below at face f
                d = step(d, t[:, f - 1], Sf[:, f], Sf[:, f - 1])   # through cell f-1
            u = Sf[:, 0].copy()
            for l in range(nlay):
                u = step(u, t[:, l], Sf[:, l], Sf[:, l + 1])
                u = u + c[l + 1]
        else:
            u = up[-1]
        F += 2*np.pi*w_*mu_*np.einsum('qbg,g->qb', u, ck.gw)
    res[sph] = F
d = np.loadtxt('/viper/ptmp2/jinma/w121prod_0929/%s/diag/olr_00600.txt' % od)
idx = {(a, b, c_): i for i, (a, b, c_) in enumerate(zip(st['gid'], st['k'], st['j']))}
ii = np.array([idx[(int(a), int(b), int(c_))] for a, b, c_ in d[:, :3]])
Fo = np.zeros((ncol, 11))
Fo[ii] = d[:, 5:]
rt = rf[-1]
Lo = (Fo*dOm[:, None]).sum(0)*rt**2
for sph, F in res.items():
    L = (F*dOm[:, None]).sum(0)*rt**2
    rq = F.sum(1)/Fo.sum(1)
    print('%s %s: L = %.4e (olr %.4e) ratio %.4f; per band %s; per column ratio median %.4f '
          '5-95%% %.4f..%.4f' % (arm, 'GCM spherical 2-stream' if sph else 'plane-parallel 2-stream',
                                 L.sum(), Lo.sum(), L.sum()/Lo.sum(), np.round(L/Lo, 3),
                                 np.median(rq), *np.percentile(rq, [5, 95])))
x = Af[-1]/Af
print('beta_top %.4f; A_top/A(r) at R_p: %.3f -> 2x/(1+x) = %.3f' % (
    beta[-1], (rt/rc.RP)**2, 2*(rt/rc.RP)**2/(1 + (rt/rc.RP)**2)))
