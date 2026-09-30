"""products.py: instrument phase curves, fits, T_b, T_day/T_night, epsilon, and the old-vs-new
ablation chain, from the emis11 / emis_prt outputs.

  python products.py            -> out/products.txt, out/curves.npz, out/*.png

Arms (each a disk flux x d^2 P(phase, lambda-bin) plus a recipe):
  OLD   olr rt_Fb (GCM two-stream), Lambertian, 11 bands, piecewise-Planck fill, BB star
        (= synth_rot300/synth.py phase on olr_00600)
  A1    this RT, GCM 11-band opacity, Lambertian from the new column flux, old recipe
  A2    this RT, GCM 11-band opacity, ANGLE-DEPENDENT intensity, old recipe
  A3b   this RT, pRT R1000 opacity + FastChem, angle-dependent, BINNED to the 11 bands, old recipe
  A3    same at full R1000 resolution, top-hat energy windows (old windows), BB star
  A3t   + real throughputs, photon weighting, BB star
  NEW   + PHOENIX star (T_eff 6628, log g 4.24, [Fe/H] +0.13)
  NEWc  (1x) NEW with chemistry at min(T, 2000 K) (no thermal dissociation / ionisation)
  NEWl  NEW but Lambertian (I = F/pi from the pRT column flux) -> angle effect at full res
"""
import numpy as np
from astropy.io import fits
import rtcore as rc

SD = rc.SD
PH = '/viper/ptmp2/jinma/prt_data/phoenix/'
OLR = {'w1x': '/viper/ptmp2/jinma/w121prod_0929/synth_rot300/diag/olr_00600.txt',
       'w10x': '/viper/ptmp2/jinma/w121prod_0929/synth_rot300_10x/diag/olr_00600.txt'}
# old top-hat windows (synth.py INSTR); SOSS o2 and G141 added
WIN = {'NRS1': (2.70, 3.72), 'NRS2': (3.82, 5.15), 'SOSS1': (0.85, 2.85), 'SOSS2': (0.6, 0.85),
       'WFC3': (1.12, 1.64), 'IRAC1': (3.13, 3.96), 'IRAC2': (3.92, 5.06)}
T0 = rc.TS/np.sqrt(rc.AOR)


# ---------------------------------------------------------------- star
def phoenix_bins(le_um):
    """PHOENIX-ACES-AGSS-COND-2011 HiRes surface flux F_lambda [erg/s/cm2/cm], trilinear in
    (T_eff 6600/6700, log g 4.0/4.5, [Fe/H] 0/+0.5) at (6628, 4.24, 0.13), averaged over bins
    with edges le_um (ascending); beyond 5.5 um a BB at T_eff scaled to 5.0-5.49 um."""
    w = fits.getdata(PH + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')*1e-4      # A -> um
    ft, fg, fz = (rc.TS - 6600)/100, (rc.LOGG_S - 4.0)/0.5, rc.FEH_S/0.5
    F = 0.0
    for T, a in ((6600, 1 - ft), (6700, ft)):
        for g, b in (('4.00', 1 - fg), ('4.50', fg)):
            for z, c in (('-0.0', 1 - fz), ('+0.5', fz)):
                F = F + a*b*c*fits.getdata(
                    PH + 'lte%05d-%s%s.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' % (T, g, z)
                ).astype(float)
    bol = np.trapezoid(F, w*1e-4)
    tail = np.pi*rc.planck_lam(np.geomspace(5.5, 1e4, 4000)*1e-4, rc.TS)
    m = (w > 5.0) & (w < 5.49)
    scl = F[m].mean()/(np.pi*rc.planck_lam(w[m]*1e-4, rc.TS)).mean()
    bol += scl*np.trapezoid(tail, np.geomspace(5.5, 1e4, 4000)*1e-4)
    cF = np.concatenate([[0], np.cumsum(0.5*(F[1:] + F[:-1])*np.diff(w))])
    out = np.zeros(len(le_um) - 1)
    for i in range(len(out)):
        l0, l1 = le_um[i], le_um[i + 1]
        if l1 <= 5.49:
            out[i] = (np.interp(l1, w, cF) - np.interp(l0, w, cF))/(l1 - l0)
        else:
            out[i] = scl*np.pi*rc.planck_lam(0.5*(l0 + l1)*1e-4, rc.TS)
    return out, bol/(rc.sig*rc.TS**4), scl


# ---------------------------------------------------------------- throughputs
def throughputs(lam):
    th = {}
    d = np.loadtxt(SD + 'thru/JWST_NIRSpec.G395H_F290LP.dat')
    g395 = np.interp(lam, d[:, 0]*1e-4, d[:, 1], left=0, right=0)
    th['NRS1'] = g395*((lam >= 2.70) & (lam <= 3.72))
    th['NRS2'] = g395*((lam >= 3.82) & (lam <= 5.15))
    h = fits.open(SD + 'thru/jwst_niriss_spectrace_0023.fits')
    for o, (l0, l1) in ((1, (0.85, 2.85)), (2, (0.6, 0.85))):
        t = np.interp(lam, h[o].data['WAVELENGTH'], h[o].data['THROUGHPUT'], left=0, right=0)
        th['SOSS%d' % o] = t*((lam >= l0) & (lam <= l1))
    d = np.loadtxt(SD + 'thru/HST_WFC3_IR.G141.dat')
    th['WFC3'] = np.interp(lam, d[:, 0]*1e-4, d[:, 1], left=0, right=0)*(
        (lam >= 1.12) & (lam <= 1.64))
    for i in (1, 2):
        d = np.loadtxt(SD + 'thru/Spitzer_IRAC.I%d.dat' % i)
        th['IRAC%d' % i] = np.interp(lam, d[:, 0]*1e-4, d[:, 1], left=0, right=0)
    return th


# ---------------------------------------------------------------- old recipe (synth.py)
def planck_int(l0, l1, T, n=400):
    lam = np.linspace(l0, l1, n)*1e-4
    T = np.atleast_1d(np.asarray(T, float))
    return np.trapezoid(rc.planck_lam(lam[None, :], T[:, None]), lam, axis=1)


def tb_invert(target, l0, l1):
    target = np.atleast_1d(target)
    lo, hi = np.full(target.shape, 100.0), np.full(target.shape, 2e4)
    for _ in range(60):
        mid = np.sqrt(lo*hi)
        f = planck_int(l0, l1, mid) - target
        lo = np.where(f < 0, mid, lo)
        hi = np.where(f >= 0, mid, hi)
    return np.sqrt(lo*hi)


def old_recipe(P, edges):
    """P (nph, nb) band disk flux x d^2 in TABLE order (descending edges)."""
    blo, bhi = edges[1:], edges[:-1]
    nb = P.shape[1]
    Tbb = np.zeros_like(P)
    for b in range(nb):
        Tbb[:, b] = tb_invert(P[:, b]/(np.pi*rc.RP**2)/np.pi, blo[b], bhi[b])
    cur = {'bol': P.sum(1)}
    for nm, (l0, l1) in WIN.items():
        Pi = np.zeros(P.shape[0])
        for b in range(nb):
            a, z = max(l0, blo[b]), min(l1, bhi[b])
            if z > a:
                Pi += P[:, b]*planck_int(a, z, Tbb[:, b])/planck_int(blo[b], bhi[b], Tbb[:, b])
        cur[nm] = Pi/(np.pi*planck_int(l0, l1, rc.TS)[0]*rc.RS**2)
    return cur


# ---------------------------------------------------------------- full-resolution recipes
def full_recipe(P, lam, le, weights, Fstar):
    """P (nph, nf) per unit wavelength [x d^2]; weights{nm: w(bin)}; Fstar (nf) surface F_lambda."""
    dl = np.abs(np.diff(le))*1e-4
    cur = {'bol': (P*dl).sum(1)}
    for nm, w in weights.items():
        cur[nm] = (P*w*dl).sum(1)/(rc.RS**2*(Fstar*w*dl).sum())
    return cur


def bol_tail(P, lam, le):
    """add lambda > 28 um: pi B(T_b) beyond, T_b from the 23-28 um disk flux."""
    dl = np.abs(np.diff(le))*1e-4
    m = lam > 23.0
    ll = np.geomspace(28.0, 2000.0, 3000)*1e-4
    out = np.zeros(P.shape[0])
    for i in range(P.shape[0]):
        target = (P[i, m]*dl[m]).sum()/(np.pi*rc.RP**2)/np.pi
        T = tb_invert(np.array([target]), lam[m][0], lam[m][-1])[0]
        out[i] = np.pi*rc.RP**2*np.pi*np.trapezoid(rc.planck_lam(ll, T), ll)
    return out


# ---------------------------------------------------------------- statistics
def stats(phi, y):
    A = np.stack([np.ones_like(phi), np.cos(2*np.pi*phi), np.sin(2*np.pi*phi),
                  np.cos(4*np.pi*phi), np.sin(4*np.pi*phi)], 1)
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    pf = np.linspace(0, 1, 7201)
    Af = np.stack([np.ones_like(pf), np.cos(2*np.pi*pf), np.sin(2*np.pi*pf),
                   np.cos(4*np.pi*pf), np.sin(4*np.pi*pf)], 1)
    yf = Af @ c
    pk = pf[np.argmax(yf)]
    # first harmonic alone: c1 cos + s1 sin, max at 2 pi phi = atan2(s1, c1)
    ph1 = (np.arctan2(c[2], c[1])/(2*np.pi)) % 1.0
    i5 = np.argmin(abs(phi - 0.5))
    return dict(day=y[i5], night=y[0], ampl=yf.max() - yf.min(),
                off=(0.5 - pk)*360.0, off1=(0.5 - ph1)*360.0, A1=2*np.hypot(c[1], c[2]),
                c=c)


def tb_band(y, w, lam, le, Fstar):
    """brightness T: pi <B>_w (T_b) = y <F*>_w (R*/R_p)^2 (the observers' definition)."""
    dl = np.abs(np.diff(le))*1e-4
    Fs = (Fstar*w*dl).sum()/(w*dl).sum()
    target = np.atleast_1d(y)*Fs*(rc.RS/rc.RP)**2
    lo, hi = np.full(target.shape, 100.0), np.full(target.shape, 2e4)
    m = w > 0
    for _ in range(60):
        mid = np.sqrt(lo*hi)
        Bm = (np.pi*rc.planck_lam(lam[m][None, :]*1e-4, mid[:, None])*(w*dl)[m]).sum(1)/(
            (w*dl)[m].sum())
        lo = np.where(Bm < target, mid, lo)
        hi = np.where(Bm >= target, mid, hi)
    return np.sqrt(lo*hi)


def teff(Pbol, i):
    return (Pbol[i]/(rc.RP**2*rc.sig))**0.25


def eps_ab(Td, Tn):
    """Cowan & Agol (2011): T_d^4 = T0^4 (1-A)(2/3 - 5 eps/12), T_n^4 = T0^4 (1-A) eps/4."""
    eps = 8.0/(3.0*(Td/Tn)**4 + 5.0)
    AB = 1.0 - 1.5*(Td**4 + 5.0/3.0*Tn**4)/T0**4
    return eps, AB


def read_olr(fn, st):
    d = np.loadtxt(fn)
    idx = {(a, b, c): i for i, (a, b, c) in enumerate(zip(st['gid'], st['k'], st['j']))}
    ii = np.array([idx[(int(a), int(b), int(c))] for a, b, c in d[:, :3]])
    F = np.zeros((len(st['gid']), d.shape[1] - 5))
    F[ii] = d[:, 5:]
    return F


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lines = []
    curves = {}
    phi, obs = rc.phases(72)
    for arm, met in (('w1x', '1x'), ('w10x', '10x')):
        st = np.load(SD + 'state/%s.npz' % arm)
        e11 = np.load(SD + 'out/emis11_%s.npz' % arm)
        rtop = float(e11['rtop'])
        mup = np.maximum(st['X'] @ obs.T, 0)
        Fo = read_olr(OLR[arm], st)
        Pold = np.einsum('qb,qp->pb', Fo*rtop**2/np.pi, mup*st['dOm'][:, None])
        edges = e11['edges']
        arms = {'OLD': old_recipe(Pold, edges), 'A1': old_recipe(e11['Plam'], edges),
                'A2': old_recipe(e11['Pang'], edges)}
        tabs = [met + '_eq'] + (['1x_nodiss'] if arm == 'w1x' else [])
        for tab in tabs:
            fn = SD + 'out/prt_%s_%s.npz' % (arm, tab)
            try:
                pr = np.load(fn)
            except FileNotFoundError:
                lines.append('# missing %s' % fn)
                continue
            lam, le = pr['lam'], pr['lam_edges']
            if le[1] < le[0]:
                lam, le = lam[::-1], le[::-1]
                Pang = pr['Pang'][:, ::-1]
                Fcol = pr['Fcol'][:, ::-1].astype(float)
            else:
                Pang, Fcol = pr['Pang'], pr['Fcol'].astype(float)
            dl = np.diff(le)*1e-4
            Fst, bolchk, scl = phoenix_bins(le)
            Fbb = np.pi*rc.planck_lam(lam*1e-4, rc.TS)
            th = throughputs(lam)
            tophat = {nm: ((lam >= l0) & (lam <= l1)).astype(float) for nm, (l0, l1) in WIN.items()}
            photon = {nm: t*lam for nm, t in th.items()}
            Plam = np.einsum('qf,qp->pf', Fcol*rtop**2/np.pi, mup*st['dOm'][:, None])
            tag = '' if tab.endswith('eq') else 'c'
            if tag == '':
                # 11-band binning of the pRT spectrum (table order = descending edges)
                P11 = np.zeros((len(phi), 11))
                for b in range(11):
                    lo_, hi_ = edges[b + 1], edges[b]
                    m = (lam >= lo_) & (lam < hi_)
                    P11[:, b] = (Pang[:, m]*dl[m]).sum(1)
                # bands outside 0.3-28 um: take the 11-band angle-dependent value (A2)
                for b in (0, 10):
                    P11[:, b] = e11['Pang'][:, b]
                arms['A3b'] = old_recipe(P11, edges)
                arms['A3'] = full_recipe(Pang, lam, le, tophat, Fbb)
                arms['A3t'] = full_recipe(Pang, lam, le, photon, Fbb)
                arms['NEWl'] = full_recipe(Plam, lam, le, photon, Fst)
            arms['NEW' + tag] = full_recipe(Pang, lam, le, photon, Fst)
            tail = bol_tail(Pang, lam, le)
            arms['NEW' + tag]['bol'] = arms['NEW' + tag]['bol'] + tail
            if tag == '':
                arms['NEWl']['bol'] = arms['NEWl']['bol'] + bol_tail(Plam, lam, le)
                arms['A3']['bol'] = arms['NEW']['bol']
                arms['A3t']['bol'] = arms['NEW']['bol']
                lines.append('# %s PHOENIX check: int F dlambda / sigma T_eff^4 = %.4f (tail '
                             'scale %.3f); bol tail > 28 um = %.2f %% (day) %.2f %% (night)' % (
                                 arm, bolchk, scl, 100*tail[36]/arms['NEW']['bol'][36],
                                 100*tail[0]/arms['NEW']['bol'][0]))
                spec = dict(lam=lam, day=Pang[36]/(rc.RS**2*Fst), night=Pang[0]/(rc.RS**2*Fst))
                curves[arm + '_spec'] = spec
                WSTAR = dict(lam=lam, le=le, Fst=Fst, Fbb=Fbb, th=th)
        curves[arm] = arms
        lines.append('## %s (rot-300 bin 00150, t = %.6e s)' % (arm, st['time']))
        lines.append('%-5s %-6s %9s %9s %9s %8s %8s %7s %7s' % (
            'arm', 'band', 'day ppm', 'night ppm', 'ampl ppm', 'off2 E', 'off1 E', 'Tb_day',
            'Tb_ngt'))
        for an, cur in arms.items():
            for nm in ('NRS1', 'NRS2', 'SOSS1', 'SOSS2', 'WFC3', 'IRAC1', 'IRAC2'):
                y = cur[nm]
                s = stats(phi, y)
                if an in ('A3t', 'NEW', 'NEWc', 'NEWl'):
                    Fs_ = WSTAR['Fst'] if an.startswith('NEW') else WSTAR['Fbb']
                    Td_, Tn_ = tb_band(np.array([s['day'], s['night']]), photon[nm], WSTAR['lam'],
                                       WSTAR['le'], Fs_)
                elif an == 'A3':
                    Td_, Tn_ = tb_band(np.array([s['day'], s['night']]), tophat[nm],
                                       WSTAR['lam'], WSTAR['le'], WSTAR['Fbb'])
                else:
                    l0, l1 = WIN[nm]
                    sb = planck_int(l0, l1, rc.TS)[0]
                    Td_, Tn_ = tb_invert(np.array([s['day'], s['night']])*sb/(rc.RP/rc.RS)**2,
                                         l0, l1)
                lines.append('%-5s %-6s %9.1f %9.1f %9.1f %8.2f %8.2f %7.0f %7.0f' % (
                    an, nm, 1e6*s['day'], 1e6*s['night'], 1e6*s['ampl'], s['off'], s['off1'],
                    Td_, Tn_))
            Pb = cur['bol']
            Td, Tn = teff(Pb, 36), teff(Pb, 0)
            eps, AB = eps_ab(Td, Tn)
            sb_ = stats(phi, Pb/Pb.max())
            lines.append('%-5s bol    T_day %.0f  T_night %.0f  Td/Tn %.3f  eps %.3f  A_B %.3f  '
                         'L_disk-avg %.4e erg/s  off2 %.2f E' % (
                             an, Td, Tn, Td/Tn, eps, AB, 4*np.pi*Pb.mean(), sb_['off']))
    open(SD + 'out/products.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    np.savez(SD + 'out/curves.npz', phi=phi, **{
        '%s__%s__%s' % (a, an, nm): v for a in ('w1x', 'w10x') if a in curves
        for an, cur in curves[a].items() for nm, v in cur.items()})
    # plots
    obsd = {'NRS1': (3924, 136), 'NRS2': (4924, 630), 'SOSS1': (1156, None),
            'SOSS2': (363, None), 'IRAC1': (4077, 553), 'IRAC2': (5121, 1045)}
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for i, nm in enumerate(('NRS1', 'NRS2', 'SOSS1', 'WFC3', 'IRAC1', 'IRAC2')):
        a = ax.flat[i]
        for arm, ls in (('w1x', '-'), ('w10x', '--')):
            if arm not in curves:
                continue
            for an, col in (('OLD', '0.5'), ('NEW', 'C3')):
                if an in curves[arm]:
                    a.plot(phi, 1e6*curves[arm][an][nm], ls, color=col, label='%s %s' % (arm, an))
        if nm in obsd:
            a.plot([0.5], [obsd[nm][0]], 'k*', ms=10)
            if obsd[nm][1]:
                a.plot([0, 1], [obsd[nm][1]]*2, 'kx', ms=8)
        a.set_title(nm)
        a.set_xlabel('phase')
        a.set_ylabel('Fp/F* [ppm]')
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(SD + 'out/phase_curves_new.png', dpi=110)
    fig, a = plt.subplots(figsize=(10, 4.5))
    for arm, ls in (('w1x', '-'), ('w10x', '--')):
        k = arm + '_spec'
        if k in curves:
            a.semilogx(curves[k]['lam'], 1e6*curves[k]['day'], ls, color='C3', lw=0.6,
                       label=arm + ' day (phase 0.5)')
            a.semilogx(curves[k]['lam'], 1e6*curves[k]['night'], ls, color='C0', lw=0.6,
                       label=arm + ' night (phase 0)')
    a.set_xlim(0.5, 12)
    a.set_xlabel('lambda [um]')
    a.set_ylabel('Fp/F* [ppm] (PHOENIX star)')
    a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(SD + 'out/spectra_day_night.png', dpi=110)
