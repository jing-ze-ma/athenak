#!/usr/bin/env python3
"""diag_wall.py state.npz [extra toyhd args]: one step from a saved state, decompose the wall-cell
radiation flux and the energy/momentum billing per term."""
import sys, numpy as np
import toyhd
from toy import CRAT, PRAT, DAY, T_UNIT
st = sys.argv[1]; extra = sys.argv[2:]
a = toyhd.argparse.Namespace()
p = None
# build args via the toyhd parser
import argparse
sys.argv = ['toyhd'] + extra
parser_args = None
def get_args(argv):
    import toyhd as T
    pr = argparse.ArgumentParser()
    # reuse main()'s parser by calling it with --days 0 would run; replicate minimal by parsing through main's code
    return None
# simplest: monkeypatch RHD.step loop by calling main with days tiny and debug on
toyhd_main = toyhd.main
captured = {}
orig_step = toyhd.RHD.step
def step(self):
    orig_step(self)
    captured['s'] = self
toyhd.RHD.step = step
toyhd.main(['--load', st, '--days', '1e-9', '--every', '1', '--debug', '--out', 'scratch/diag'] + extra)
s = captured['s']; d = s.dbg
w, mu = s.w, s.mu
cdt = d['cdt']; tc = d['tc']; I = d['I']; chi = d['chi']
E = I @ w; F = I @ (w * mu)
Lref = s.Lref
def L(x, i): return x * s.rv[i]**2 / Lref
print('wall cells after one dt step (stage 2):')
for i in range(5):
    b = d['v'][i] / CRAT
    # thick-limit decomposition of the cell flux: F = (4/3) beta E  + [transport-driven]
    Fadv = 4.0 / 3.0 * b * E[i]
    # transport-driven anisotropy: per angle  I - S/tc^4 ~ -(trans - advdiv + I - I_ini)/(cdt chi)
    rhs = -(d['trans'][i] - d['advdiv'][i] + (I[i] - d['I_ini'][i]))
    Ftr = (rhs / (cdt * chi[i])) @ (w * mu)
    rhs_t = -(d['trans'][i]); rhs_a = d['advdiv'][i]; rhs_dt = -(I[i] - d['I_ini'][i])
    print(f" i{i} v {d['v'][i]*36.618:+.3f} km/s  F/L {L(F[i],i):+8.1f}  (4/3)vE/L {L(Fadv,i):+8.1f}  "
          f"transport-part {L((rhs_t/(cdt*chi[i]))@(w*mu),i):+8.1f}  adv-part {L((rhs_a/(cdt*chi[i]))@(w*mu),i):+8.1f}  "
          f"dI/dt-part {L((rhs_dt/(cdt*chi[i]))@(w*mu),i):+8.1f}   E {E[i]:.4e} Tgas {d['T'][i]*1e5:.4e}")
# the I_ini anisotropy: the previous step's flux
Fini = d['I_ini'] @ (w * mu)
print('F of I_ini (start of step) /L:', [f"{L(Fini[i],i):+.1f}" for i in range(5)])
print('dF_src (billed) * prat/crat per cell (momentum to gas per step):', (-PRAT * d['dF'][:5] / CRAT))
print('dE_src*prat (energy to gas per step):', (-PRAT * d['dE'][:5]))
print('gas eint per cell:', (d['rho'] * d['Tg0'] / (5/3-1))[:5])
# wall-face decomposition (test a)
Lf = lambda x: x * s.rf[0]**2 / Lref / CRAT
Ac, Bc = d['Ac'], d['Bc']
Fdo = (Ac[0] * d['Ig'] + Bc[0] * I[0]) @ w        # implicit DO part of the face-0 energy flux
Fad = s.Fadv[0]                                   # explicit advective part adv_face * I_upwind
print(f"face 0: total Lw {Lf(Fdo+Fad):+.1f} L = DO part {Lf(Fdo):+.1f} + explicit advective part {Lf(Fad):+.1f};"
      f" v_face(adv) = {d['adv'][0]*36.618:+.3f} km/s; (4/3)E0 v_face A / L = {Lf(4/3*E[0]*d['adv'][0]):+.1f};"
      f" half of cell-0 Doppler flux (2/3) v*_0 E0 = {Lf(CRAT*2/3*d['v'][0]/CRAT*E[0]):+.1f} L; half ghost F_g = {Lf(CRAT*0.5*s.Fg[0]):+.1f}")
