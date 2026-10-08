"""Available gas (line + molecular + continuum) radiation force on the observed Betelgeuse
profile.  usage:  python gas.py TABLE SET PROF TCASE MASS
  TABLE: A   (DACE low-P table A, Premixed_1x_g8_11_hiT2_dace_lowP_SiOfc + lowP CE table)
         old (original Premixed_1x_g8_11_hiT2, clamped at its 1e-8 bar edge)
  SET A|B, PROF dent|harper, TCASE obs | <T [K]> (constant T for r >= 1.2 R*: cool bracket)
(1) static ck column (rt.py logic, LTE S = B(T), T held fixed = observed / bracket):
    grey-Lucy hydrostatic photosphere (rt.hse, s = 1) below, observed rho beyond the radius
    where it exceeds the hydrostatic one; formal solution per (band, g); Gamma_F = k_F/k_E.
(2) thin flux-mean Gamma and Sobolev Gamma (sobolev.py) on 1.2-30 R* with |dv/dr| from the
    steady wind (Mdot 1, 2, 4e-6) and the turbulent bracket v_t/H_rho (v_t 5, 10 km/s);
    flux A (unattenuated) and B (window-attenuated from 1.2 R*).
"""
import os
import sys
tab = sys.argv[1]
if tab == 'A':
    os.environ['CK_KTABLE'] = 'Premixed_1x_g8_11_hiT2_dace_lowP_SiOfc.txt'
    os.environ['CK_CETABLE'] = 'FastChem_ck_1x_int_hiT2_lowP.txt'
else:
    os.environ['CK_KTABLE'] = 'Premixed_1x_g8_11_hiT2.txt'
    os.environ['CK_CETABLE'] = 'FastChem_ck_1x_int_hiT2.txt'
os.environ['CK_DATA'] = '/orion/ptmp/jinma/rsg_wind_1008/dace/data/'
import numpy as np  # noqa: E402
import rsglib as L  # noqa: E402
import rt  # noqa: E402
import sobolev as S  # noqa: E402
import prof as P  # noqa: E402
import wind as W  # noqa: E402

OUT = '/orion/ptmp/jinma/rsg_wind_1008/betelgeuse/out/'
MODE = 'clamp'
N_LI = 14
NHSE = int(os.environ.get("NHSE", 100000))
RE_ALL = False


def column(st, prof, Tc):
    a = rt.hse(st, 1., MODE)
    R = st['R']
    xo = np.exp(np.linspace(np.log(1.0), np.log(30.), 260))
    ro = P.rho(xo, prof, st['d'])
    To = P.Tobs(xo, prof) if Tc is None else np.where(xo >= 1.2, Tc, P.Tobs(xo, prof))
    rh = np.exp(np.interp(xo*R, a['r'], np.log(a['rho']), right=-300.))
    # join: first radius beyond R where the observed density exceeds the hydrostatic one
    j = np.where((xo >= 1.0) & (ro > rh))[0][0]
    xj = xo[j]
    keep_h = (a['r'] < xj*R) & (a['tau'] < 20.)   # drop the H- runaway base point
    rH, TH, rhoH = a['r'][keep_h], a['T'][keep_h], a['rho'][keep_h]
    if len(rH) > NHSE:
        sel = np.unique(np.concatenate([np.linspace(0, len(rH) - 1, NHSE).astype(int)]))
        rH, TH, rhoH = rH[sel], TH[sel], rhoH[sel]
    r = np.concatenate([rH, xo[j:]*R])
    T = np.concatenate([TH, To[j:]])
    rho = np.concatenate([rhoH, ro[j:]])
    N = len(r)
    nh = len(rH)
    # lambda iterations to LTE radiative equilibrium (rt.run logic) in the hydrostatic
    # photosphere only (observed / bracket T held fixed above), or everywhere for 're'
    free = np.arange(N) < nh if not RE_ALL else np.ones(N, bool)
    for it in range(N_LI + 1):
        K, kR, kP = rt.opac(T, rho, MODE)
        chi = (K*rho[:, None, None]).reshape(N, -1)
        Sb = np.repeat(rt.Bband(T), L.NG, axis=1)
        J, H = rt.formal(r, chi, Sb)
        dtr = np.zeros(N)
        dtr[:-1] = 0.5*(kR[1:]*rho[1:] + kR[:-1]*rho[:-1])*np.diff(r)
        tauR = np.cumsum(dtr[::-1])[::-1]
        if it == N_LI:
            break
        w = np.tile(L.GW, L.NB)
        kq = K.reshape(N, -1)
        Tn = T.copy()
        for i in np.where((tauR < rt.TAU_RE) & free)[0]:
            rhs = (w*kq[i]*J[i]).sum()
            kb = (w*kq[i]).reshape(L.NB, L.NG).sum(1)
            lo, hi = 300., 12000.
            for _ in range(40):
                tm = 0.5*(lo + hi)
                f = (kb*rt.Bband(np.array([tm]))[0]).sum() - rhs
                lo, hi = (lo, tm) if f > 0 else (tm, hi)
            Tn[i] = 0.5*(lo + hi)
        T = np.sqrt(Tn*T)
    w = np.tile(L.GW, L.NB)
    kq = K.reshape(N, -1)
    Fq = 4*np.pi*H
    F = (w*Fq).sum(1)
    kF = (w*kq*Fq).sum(1)/F
    Fb = (w*Fq).reshape(N, L.NB, L.NG).sum(2)
    kFb = (w*kq*Fq).reshape(N, L.NB, L.NG).sum(2)/F[:, None]/st['kE']
    return dict(cGFb=kFb, cr=r/R, crho=rho, cT=T, ctauR=tauR, ckR=kR, cGF=kF/st['kE'],
                cFratio=F/(st['L']/(4*np.pi*r**2)), cGthin=np.array(
                    [L.kappa_thin_flux(K[i], st['Teff']) for i in range(N)])/st['kE'],
                cFb=Fb/F[:, None], xjoin=xj)


def sobolev_part(st, prof, Tc):
    x = W.xgrid(1.2, 30., 240)
    out = {}
    for md in (1e-6, 2e-6, 4e-6):
        w = W.steady(st, prof, md, x=x, Tfix=Tc)
        if md == 1e-6:
            rho, T = w['rho'], w['T']
            K, kR = S.opacities(T, rho, MODE)
            out['sx'], out['srho'], out['sT'] = x, rho, T
            out['Gthin'] = S.gamma_thin(st, K)
            out['Gthin_b'] = (L.fband(st['Teff'])[0]*(L.GW*K).sum(2))/st['kE']
            Hr = w['Hrho']
        # the wind density differs from rho_obs only beyond x_vinf; use each wind's own rho
        Kw = K if np.isnan(w['x_vinf']) else S.opacities(w['T'], w['rho'], MODE)[0]
        for fl in ('A', 'B'):
            out[f'Gsob_w{md:.0e}_{fl}'] = S.gamma_sob(st, w['r'], w['rho'], w['T'],
                                                      w['dvdr'], MODE, flux=fl, K=Kw,
                                                      r_ph=w['r'][0])[0]
    r = x*st['R']
    for vt in (5e5, 10e5):
        dv = vt/np.abs(Hr)
        for fl in ('A', 'B'):
            out[f'Gsob_t{vt/1e5:.0f}_{fl}'] = S.gamma_sob(st, r, rho, T, dv, MODE, flux=fl,
                                                          K=K, r_ph=r[0])[0]
    return out


if __name__ == '__main__':
    _, tab, sset, prof, tcase, M = sys.argv
    RE_ALL = tcase == 're'
    Tc = None if tcase in ('obs', 're') else float(tcase)
    st = P.stellar(sset, float(M))
    st['name'] = 'betelgeuse'
    res = dict(sobolev_part(st, prof, Tc))
    res.update(column(st, prof, Tc))
    os.makedirs(OUT, exist_ok=True)
    fn = OUT + f'gas_{tab}_{sset}{M}_{prof}_{tcase}.npz'
    np.savez(fn, **res)
    xs = [1.2, 1.5, 2, 3, 5, 10, 20]
    ii = [np.argmin(abs(res['sx'] - v)) for v in xs]
    jj = [np.argmin(abs(res['cr'] - v)) for v in xs]
    print(fn, 'join', round(res['xjoin'], 3), 'F/F* top', round(res['cFratio'][-1], 3))
    print(' static GF ', np.round(res['cGF'][jj], 3))
    print(' thin      ', np.round(res['Gthin'][ii], 3))
    print(' sob 2e-6 A', np.round(res['Gsob_w2e-06_A'][ii], 3))
    print(' sob t10 A ', np.round(res['Gsob_t10_A'][ii], 3))
