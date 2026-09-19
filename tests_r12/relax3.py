#!/usr/bin/env python3
"""tests_r12/relax3.py (2026-09-19): make the DEEP RADIATIVE ZONE of the IC carry L under
the code's own operator.  e_ledger_shell.txt (face luminosities integrated by the sweep)
shows faces 1-14 (0.51-0.62 R) carry 1.0047 L and faces 15-16 1.012/1.025 L while the MLT
closure above is pinned to F_req = L: ~0.5 % L is dumped at the base of the MLT zone.
In the diffusion regime F ~ dT^4/dr, so below the anchor (first face with l <= 1, the MLT
base) the IC's dT^4/dr is divided by l(r) = L_face/L; T is unchanged at and above the
anchor.  rho is re-solved DOWNWARD from the anchor as a difference against the old IC
(same trick as relax2.py).  Usage: relax3.py <ledger_arm> <out_tag> [--start ic]"""
import sys
import numpy as np
import relax_ic as R
import relax2


def main():
    arm, tag = sys.argv[1], sys.argv[2]
    ic = relax2.arg('--start', R.IC0)
    fn = arm + '/e_ledger_shell.txt'
    t = float(open(fn).readline().split('=')[1].split()[0])
    lf = np.loadtxt(fn)[:, 2]/(R.LSTAR*t)
    rc = relax2.rd(arm + '/rt_profile.bin')[0][1]
    x1min, x1max = 1.18585e11, 2.4057e11
    rf = np.concatenate(([x1min], 0.5*(rc[1:] + rc[:-1]), [x1max]))   # ~faces
    lf[0] = lf[1]                         # the wall face is the BC's own artefact
    ia = int(np.argmax(lf[1:] <= 1.0)) + 1
    print('anchor face %d at %.4f R ; l below it: min %.4f max %.4f'
          % (ia, rf[ia]/R.RSTAR, lf[:ia].min(), lf[:ia].max()))
    d = np.loadtxt(ic, comments='#')
    r, rho, eint = d[:, 0], d[:, 1], d[:, 2]
    T = R.temp_of(rho, eint)
    # a CONSTANT deep factor (mean of faces 1-14), ramped to 1 over 0.60-0.64 R: the
    # 1.01-1.025 bump on the two faces under the MLT base is left to adjust by itself
    ldeep = float(lf[1:15].mean())
    ell = 1.0 + (ldeep - 1.0)*(1.0 - R.smooth((r - 0.60*R.RSTAR)/(0.04*R.RSTAR)))
    print('deep factor %.5f' % ldeep)
    k0 = int(np.searchsorted(r, rf[ia]))
    T4 = T**4
    T4n = T4.copy()
    for i in range(k0 - 1, -1, -1):
        T4n[i] = T4n[i+1] + (T4[i] - T4[i+1])/(0.5*(ell[i] + ell[i+1]))
    Tn = T4n**0.25

    def rhs(rh, Tt, dpr, rr):
        w = R.taper_w(rh, Tt)
        return (-rh*R.GM/rr**2 - w*dpr
                + (1.0 - w)*rh*R.kappa(Tt, rh)*R.LSTAR/(4*np.pi*rr**2)/R.CLIGHT)
    dpo = np.gradient(R.A_RAD*T**4/3.0, r)
    dpn = np.gradient(R.A_RAD*Tn**4/3.0, r)
    pgo = R.pgas(rho, T)
    ro = rhs(rho, T, dpo, r)
    rn = rho.copy()
    # UPWARD from the first point: with Prad/Pgas ~ 50 the gas scale height is beta H_p,
    # the downward march is exponentially unstable and the upward one forgets its start
    # within a few beta H_p (rho g ~ -dPrad/dr locally)
    one = np.ones(1)
    rn[0] = rho[0]/ell[0]
    pg = float(R.pgas(rn[0]*one, Tn[0]*one)[0])
    kend = int(np.searchsorted(r, 0.72*R.RSTAR))
    for i in range(0, kend):
        h = r[i+1] - r[i]
        k1 = rhs(rn[i]*one, Tn[i]*one, dpn[i]*one, r[i]*one)[0]
        pe = max(pg + h*k1, 1e-3*pg)
        re = float(R.rho_of(pe, Tn[i+1], rn[i]))
        k2 = rhs(re*one, Tn[i+1]*one, dpn[i+1]*one, r[i+1]*one)[0]
        pg = max(pg + 0.5*h*(k1 + k2), 1e-3*pg)
        rn[i+1] = float(R.rho_of(pg, Tn[i+1], re))
    # the same march on the OLD T measures the offline force law's own error: remove it
    ro_ = rho.copy()
    pg = float(R.pgas(ro_[0]*one, T[0]*one)[0])
    for i in range(0, kend):
        h = r[i+1] - r[i]
        k1 = rhs(ro_[i]*one, T[i]*one, dpo[i]*one, r[i]*one)[0]
        pe = max(pg + h*k1, 1e-3*pg)
        re = float(R.rho_of(pe, T[i+1], ro_[i]))
        k2 = rhs(re*one, T[i+1]*one, dpo[i+1]*one, r[i+1]*one)[0]
        pg = max(pg + 0.5*h*(k1 + k2), 1e-3*pg)
        ro_[i+1] = float(R.rho_of(pg, T[i+1], re))
    print('  offline re-march of the OLD ic vs the ic itself: max |drho/rho| = %.3e'
          % np.abs(ro_[:kend]/rho[:kend] - 1).max())
    rn[:kend] = rho[:kend]*rn[:kend]/ro_[:kend]
    rn[kend:] = rho[kend:]
    print('  rho_new/rho at the hand-back 0.72 R: %.5f' % (rn[kend-1]/rho[kend-1]))
    for s in (0.45, 0.50, 0.55, 0.60, 0.63, 0.65, 0.70):
        k = int(np.argmin(np.abs(r/R.RSTAR - s)))
        print('  %.3f R: T_new/T %.5f  rho_new/rho %.5f' % (s, Tn[k]/T[k], rn[k]/rho[k]))
    en = eint + (R.egas(rn, Tn) + R.taper_w(rn, Tn)*R.A_RAD*Tn**4
                 - R.egas(rho, T) - R.taper_w(rho, T)*R.A_RAD*T**4)
    out = R.HERE + '/ic_rd%s.txt' % tag
    with open(out, 'w') as fh:
        fh.write('# tests_r12/relax3.py %s from %s: r[cm] rho[g/cm^3] eint[erg/cm^3]\n'
                 % (tag, arm))
        for a, b_, c_ in zip(r, rn, en):
            fh.write('%.10e %.10e %.10e\n' % (a, b_, c_))
    print('  wrote', out)


if __name__ == '__main__':
    main()
