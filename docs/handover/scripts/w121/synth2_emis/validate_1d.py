"""validate_1d.py: 1-D checks of the new RT (brief item 7b).

  python validate_1d.py w1x
(a) isothermal column (T = 1500, 2500 K, the run's rho and grid, plane-parallel source B):
    I(mu) / B_lambda(T) for all bins and mu -> must be 1.
(b) one dayside (substellar) and one nightside (antistellar) column of w1x: pRT's own
    Radtrans.calculate_flux (on-the-fly opacities + its RT, hydrostatic dtau = k dp / g) vs
    this code (premixed table 1x_eq, same dtau rule, same angle quadrature, plane-parallel).
    Reported: bolometric (0.3-28 um) and band-integrated flux ratios.
"""
import os
import sys
import numpy as np
os.environ['pRT_emcee_mode'] = 'True'
import rtcore as rc   # noqa: E402
import prt_table as pt   # noqa: E402
from emis_prt import interp   # noqa: E402

arm = sys.argv[1] if len(sys.argv) > 1 else 'w1x'
st = np.load(rc.SD + 'state/%s.npz' % arm)
t = np.load(rc.SD + 'tables/prt_1x_eq.npz')
lam = t['lam']
le = t['lam_edges']
dl = np.abs(np.diff(le))*1e-4
gw = t['gw']
out = []


def kcol(T, pbar):
    iT, fT = interp(np.log10(T), t['lT'])
    iP, fP = interp(np.log10(pbar), t['lP'])
    a, b = fT[:, None, None], fP[:, None, None]
    lk = t['lk']
    return 10.0**((1 - a)*((1 - b)*lk[iT, iP] + b*lk[iT, iP + 1])
                  + a*((1 - b)*lk[iT + 1, iP] + b*lk[iT + 1, iP + 1]))


# (a) isothermal
q = 0
dz = np.diff(st['rf'])
for Tiso in (1500.0, 2500.0):
    T = np.full(st['T'].shape[1], Tiso)
    k = kcol(T, st['p'][q]/1e6)
    dtau = k*(st['rho'][q]*dz)[:, None, None]
    B = rc.planck_lam(lam[None, :]*1e-4, T[:, None])[:, :, None]
    mus = np.array([0.05, 0.3, 0.7, 1.0])
    J = rc.formal(dtau, rc.face_source(B, dtau), mus)
    r = np.einsum('fgm,g->fm', J, gw)/B[0, :, 0][:, None]
    out.append('(a) isothermal T=%.0f K: I/B over all %d bins x mu %s: min %.6f max %.6f' % (
        Tiso, len(lam), mus, r.min(), r.max()))

# (b) pRT's own 1-D emission
from petitRADTRANS.radtrans import Radtrans   # noqa: E402
pt.cp.set_input_data_path('/viper/ptmp2/jinma/prt_data')
X = st['X']
g = 1167.8749*(1.126575e10/rc.RP)**2        # the run's point-mass g at R_p
bands = [(0.6, 0.85), (0.85, 2.85), (1.12, 1.64), (2.70, 3.72), (3.82, 5.15), (5.0, 28.0)]
for nm, q in (('day (substellar)', np.argmax(-X[:, 0])), ('night (antistellar)', np.argmax(X[:, 0]))):
    T = st['T'][q][::-1]
    pb = st['p'][q][::-1]/1e6                       # top -> bottom, increasing p
    sel = (pb > 1e-7) & (np.diff(np.concatenate([[0.0], pb])) > 0)  # sponge-top cells out
    T, pb = T[sel], pb[sel]
    rt = Radtrans(pressures=pb, wavelength_boundaries=pt.WL,
                  line_species=[s for s, _, _ in pt.LS],
                  gas_continuum_contributors=['H2--H2', 'H2--He', 'H-'],
                  rayleigh_species=['H2', 'He'], line_opacity_mode='c-k')
    mf, mu, _ = pt.chem('1x', T, pb)
    wl, F, _ = rt.calculate_flux(temperatures=T, mass_fractions=mf, mean_molar_masses=mu,
                                 reference_gravity=g)
    Fp = np.asarray(F)                              # erg/s/cm2/cm
    ang = rt.emission_angle_grid
    # this code, same layers: cell-centred T at pressure points, dtau between points by k dp/g
    k = kcol(T, pb)                                 # (n, nf, ng) top -> bottom
    kb = k[::-1]
    Tb = T[::-1]
    pbb = pb[::-1]*1e6
    # layer i spans the half-way points in p between neighbours
    pf = np.concatenate([[pbb[0]], np.sqrt(pbb[1:]*pbb[:-1]), [pbb[-1]]])
    dtau = kb*(-np.diff(pf)/g)[:, None, None]
    B = rc.planck_lam(lam[None, :]*1e-4, Tb[:, None])[:, :, None]
    J = rc.formal(dtau, rc.face_source(B, dtau), ang['cos_angles'])
    Jg = np.einsum('fgm,g->fm', J, gw)
    Fm = 2*np.pi*(Jg*(ang['cos_angles']*ang['weights'])).sum(1)
    Fm8 = None
    out.append('(b) %s column %d: T range %.0f-%.0f K, %d layers, angle grid %s' % (
        nm, q, T.min(), T.max(), len(T), np.round(ang['cos_angles'], 4)))
    out.append('    bolometric 0.3-28 um: pRT %.5e  this code %.5e  ratio %.4f' % (
        (Fp*dl).sum(), (Fm*dl).sum(), (Fm*dl).sum()/(Fp*dl).sum()))
    for l0, l1 in bands:
        m = (lam >= l0) & (lam < l1)
        out.append('    %.2f-%.2f um: ratio this/pRT %.4f   T_b-equivalent flux %.4e' % (
            l0, l1, (Fm[m]*dl[m]).sum()/(Fp[m]*dl[m]).sum(), (Fp[m]*dl[m]).sum()))
    rr = Fm/np.maximum(Fp, 1e-300)
    out.append('    per-bin ratio: median %.4f, 5-95 %% %.4f..%.4f' % (
        np.median(rr), *np.percentile(rr, [5, 95])))
open(rc.SD + 'out/validate_1d.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
