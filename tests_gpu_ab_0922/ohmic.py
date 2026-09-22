"""Which limit holds dt: the MHD fast-mode CFL, or the Ohmic cap?

Rebuilds BOTH per-cell limits of the last dump of each arm, exactly as the code
computes them:
  CFL    dt_i = cfl * dx_i / (|v_i| + cf_i)          (mhd_newdt.cpp; cf_2,3 /sin_cell)
  Ohmic  dt_i = cfl * dx_i^2 / (6 eta)   (Resistivity::NewTimeStepGeneralResist)
with eta = ResistivityEOS(x_e(rho,T), T, max_eta) and x_e from a python port of
eos_composition.hpp (xe_model.py).  Read-only.
"""
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert                                    # noqa: E402
import dhjcs                                          # noqa: E402
import common as C                                    # noqa: E402
import xe_model as X                                  # noqa: E402
import cfl as CFL                                     # noqa: E402

MAX_ETA = 1.0e13
RP = 9.44e9
ARMS = ['prodbin', 'sph', 'sphbeam', 'offfix']
XE = X.XeTable()


def analyse(name):
    path = sorted(glob.glob(dict(C.ARMS1)[name] + '/bin/*.bin'))[-1]
    raw = bin_convert.read_binary(path)
    mb = raw['mb_data']
    d = {k: np.asarray(mb[k], dtype=np.float64) for k in raw['var_names']}
    g = C.cs_metric(raw)
    geo = dhjcs.cs_geometry(raw)
    T, p, g1 = CFL.EOS.invert(d['dens'], d['eint'])
    dtc, dt1, dt2, dt3, cf1, cf2, cf3 = CFL.cell_dt(d, g, T, p, g1)
    dt_cfl = C.CFL * dtc                                    # per cell, seconds

    xe = XE(d['dens'], T)
    eta = X.eta_eos(xe, T, MAX_ETA)
    dx1 = np.broadcast_to(g['dx1'][None, None, None, :], eta.shape)
    dxmin2 = np.minimum(np.minimum(dx1**2, g['dx2']**2), g['dx3']**2)
    dt_ohm = C.CFL * dxmin2 / (6.0 * eta)

    dt_all = np.minimum(dt_cfl, dt_ohm)
    capped = eta >= 0.999999 * MAX_ETA
    lat = np.broadcast_to(np.degrees(geo['lat'])[..., None], eta.shape)
    lon = np.broadcast_to(np.degrees(geo['lon'])[..., None], eta.shape)
    rr = np.broadcast_to(g['rc'][None, None, None, :] / RP, eta.shape)
    bsq = d['bcc1']**2 + d['bcc2']**2 + d['bcc3']**2
    beta = np.where(bsq > 0, 2 * p / bsq, np.inf)
    # run's own dt at the end
    cl, _ = C.cyclelines(dict(C.ARMS1)[name])
    return dict(name=name, raw=raw, d=d, g=g, T=T, p=p, xe=xe, eta=eta, beta=beta,
                dt_cfl=dt_cfl, dt_ohm=dt_ohm, dt_all=dt_all, capped=capped,
                lat=lat, lon=lon, rr=rr, dt_run=cl[-1, 3],
                dt_run_med=np.median(cl[-200:, 3]))


def zone(lon):
    a = np.abs(lon)
    return np.where(a < 75, 'day', np.where(a > 105, 'night', 'term'))


def top_table(r, fh, n=20):
    flat = r['dt_all'].ravel()
    order = np.argsort(flat)[:n]
    idx = np.unravel_index(order, r['dt_all'].shape)
    fh.write('\n**%s** -- last dump t = %.5e, cycle %d.  Reconstructed limits: '
             'CFL-only %.2f s, Ohmic-only %.2f s, combined %.2f s; run dt %.2f s '
             '(median of last 200 cycles %.2f s).\n\n'
             % (r['name'], r['raw']['time'], r['raw']['cycle'], r['dt_cfl'].min(),
                r['dt_ohm'].min(), r['dt_all'].min(), r['dt_run'], r['dt_run_med']))
    fh.write('| # | dt [s] | binds | r/Rp | lat | lon_ss | zone | p [bar] | T [K] '
             '| rho | x_e | eta | eta/max | beta |\n')
    fh.write('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n')
    for q in range(n):
        m, k, j, i = [int(x[q]) for x in idx]
        b = 'CFL' if r['dt_cfl'][m, k, j, i] < r['dt_ohm'][m, k, j, i] else 'Ohmic'
        fh.write('| %d | %.3f | %s | %.3f | %+.1f | %+.1f | %s | %.3e | %.0f | %.2e '
                 '| %.2e | %.2e | %.3f | %.2e |\n'
                 % (q + 1, r['dt_all'][m, k, j, i], b, r['rr'][m, k, j, i],
                    r['lat'][m, k, j, i], r['lon'][m, k, j, i],
                    str(zone(r['lon'][m, k, j, i])), r['p'][m, k, j, i] / 1e6,
                    r['T'][m, k, j, i], r['d']['dens'][m, k, j, i], r['xe'][m, k, j, i],
                    r['eta'][m, k, j, i], r['eta'][m, k, j, i] / MAX_ETA,
                    r['beta'][m, k, j, i]))


def census(r, fh):
    c = r['capped']
    z = zone(r['lon'])
    tot = c.size
    fh.write('\n%s: %d of %d cells (%.2f %%) at eta = max_eta. '
             % (r['name'], c.sum(), tot, 100.0 * c.mean()))
    if c.any():
        for lab in ['day', 'term', 'night']:
            sel = c & (z == lab)
            fh.write('%s %.1f %% ' % (lab, 100.0 * sel.sum() / max(c.sum(), 1)))
        fh.write('| capped T: median %.0f K, 5-95%% %.0f-%.0f K; '
                 'r/Rp of capped: %.3f-%.3f\n'
                 % (np.median(r['T'][c]), np.percentile(r['T'][c], 5),
                    np.percentile(r['T'][c], 95), r['rr'][c].min(), r['rr'][c].max()))
    else:
        fh.write('\n')
    # the shell that actually binds
    ib = np.unravel_index(np.argmin(r['dt_all']), r['dt_all'].shape)[3]
    shell = (slice(None), slice(None), slice(None), ib)
    zs = z[shell]
    cs = c[shell]
    fh.write('  binding shell i = %d (r/Rp %.3f): %.1f %% of its cells capped '
             '(day %.1f %%, term %.1f %%, night %.1f %%); T on that shell: '
             'day %.0f K, night %.0f K\n'
             % (ib, r['g']['rc'][ib] / RP, 100.0 * cs.mean(),
                100.0 * (cs & (zs == 'day')).sum() / max((zs == 'day').sum(), 1),
                100.0 * (cs & (zs == 'term')).sum() / max((zs == 'term').sum(), 1),
                100.0 * (cs & (zs == 'night')).sum() / max((zs == 'night').sum(), 1),
                np.mean(r['T'][shell][zs == 'day']),
                np.mean(r['T'][shell][zs == 'night'])))
    return ib


def tcap_curve(rhos):
    """T at which eta(rho,T) first reaches max_eta, by bisection on log T."""
    lo = np.full_like(rhos, np.log10(200.0))
    hi = np.full_like(rhos, np.log10(1.0e5))
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        e = X.eta_eos(XE(rhos, 10**mid), 10**mid, MAX_ETA)
        capped = e >= 0.999999 * MAX_ETA
        lo = np.where(capped, mid, lo)
        hi = np.where(capped, hi, mid)
    return 10**(0.5 * (lo + hi))


if __name__ == '__main__':
    res = {a: analyse(a) for a in ARMS}
    fh = open(C.OUT + 'ohmic_tables.md', 'w')
    fh.write('# Ohmic cap vs MHD CFL, per cell (ohmic.py)\n')
    fh.write('\nmax_eta = %.1e, cfl = %.2f, eta = 230 sqrt(T)/max(x_e, 230 sqrt(T)'
             '/max_eta) + 5.2e11*20/T^1.5, capped at max_eta.\n' % (MAX_ETA, C.CFL))

    fh.write('\n## Summary\n\n')
    fh.write('| arm | run dt [s] | reconstructed dt [s] | CFL-only [s] | Ohmic-only [s]'
             ' | binds | capped cells | dt if cap removed |\n'
             '|---|---|---|---|---|---|---|---|\n')
    for a in ARMS:
        r = res[a]
        b = 'CFL' if r['dt_cfl'].min() < r['dt_ohm'].min() else 'Ohmic'
        fh.write('| %s | %.2f | %.2f | %.2f | %.2f | %s | %.2f %% | %.2f |\n'
                 % (a, r['dt_run_med'], r['dt_all'].min(), r['dt_cfl'].min(),
                    r['dt_ohm'].min(), b, 100.0 * r['capped'].mean(),
                    r['dt_cfl'].min()))

    fh.write('\n## 20 most limiting cells per arm\n')
    for a in ARMS:
        top_table(res[a], fh)

    fh.write('\n## Day/night/terminator split of the Ohmic-capped cells\n\n')
    fh.write('day |lon_ss| < 75 deg, terminator 75-105 deg, night > 105 deg.\n\n')
    ibs = {}
    for a in ARMS:
        ibs[a] = census(res[a], fh)

    fh.write('\n## Where the cap switches on\n\n')
    rhos = 10.0**np.arange(-12.0, -2.9, 0.5)
    tc = tcap_curve(rhos)
    fh.write('| rho [g/cm^3] | ' + ' | '.join('%.0e' % x for x in rhos) + ' |\n')
    fh.write('|---|' + '---|' * len(rhos) + '\n')
    fh.write('| T_cap [K] | ' + ' | '.join('%.0f' % x for x in tc) + ' |\n')
    fh.write('\nHow far the capped cells sit below their own T_cap:\n\n')
    fh.write('| arm | median T of capped | median T_cap there | median deficit [K] '
             '| deficit at the binding shell |\n|---|---|---|---|---|\n')
    for a in ARMS:
        r = res[a]
        c = r['capped']
        if not c.any():
            fh.write('| %s | - | - | - | - |\n' % a)
            continue
        tcc = tcap_curve(r['d']['dens'][c])
        ib = ibs[a]
        sh = (slice(None), slice(None), slice(None), ib)
        cs = r['capped'][sh]
        dsh = (tcap_curve(r['d']['dens'][sh][cs]) - r['T'][sh][cs]
               if cs.any() else np.array([np.nan]))
        fh.write('| %s | %.0f | %.0f | %.0f | %.0f |\n'
                 % (a, np.median(r['T'][c]), np.median(tcc),
                    np.median(tcc - r['T'][c]), np.median(dsh)))
    # ---------------------------------------------------- night-side T on the shell
    fh.write('\n## Temperature distribution on the binding shell (night side only)\n\n')
    fh.write('| arm | i | r/Rp | T_cap there [K] | night T min | 1 %% | 5 %% | median '
             '| n cells below T_cap |\n|---|---|---|---|---|---|---|---|---|\n')
    for a in ARMS:
        r = res[a]
        ib = ibs[a]
        sh = (slice(None), slice(None), slice(None), ib)
        night = zone(r['lon'][sh]) == 'night'
        tn = r['T'][sh][night]
        tcs = tcap_curve(r['d']['dens'][sh][night])
        fh.write('| %s | %d | %.3f | %.0f | %.0f | %.0f | %.0f | %.0f | %d |\n'
                 % (a, ib, r['g']['rc'][ib] / RP, np.median(tcs), tn.min(),
                    np.percentile(tn, 1), np.percentile(tn, 5), np.median(tn),
                    (tn < tcs).sum()))

    fh.write('\n## Capped cells and night-side T, shell by shell (i = 38..49)\n\n')
    fh.write('| i | r/Rp | dx_min [cm] | ' + ' | '.join(
        '%s ncap / nightT 1%%' % a for a in ARMS)
        + ' |\n|---|---|---|' + '---|' * len(ARMS) + '\n')
    for i in range(38, 50):
        cells = []
        for a in ARMS:
            r = res[a]
            sh = (slice(None), slice(None), slice(None), i)
            n = zone(r['lon'][sh]) == 'night'
            cells.append('%d / %.0f' % (r['capped'][sh].sum(),
                                        np.percentile(r['T'][sh][n], 1)))
        r0 = res['prodbin']
        fh.write('| %d | %.3f | %.3e | ' % (
            i, r0['g']['rc'][i] / RP,
            min(r0['g']['dx1'][i], r0['g']['dx2'][..., i].min(),
                r0['g']['dx3'][..., i].min())) + ' | '.join(cells) + ' |\n')

    # ---------------------------------------------------- what would change the limit
    fh.write('\n## Would raising max_eta, or warming the twilight, change the limit?\n')
    fh.write('\n### dt vs max_eta (the Ohmic dt scales as 1/eta, so RAISING it is'
             ' worse)\n\n')
    etas = [1.0e12, 3.0e12, 5.0e12, 1.0e13, 3.0e13, 1.0e14]
    fh.write('| arm | CFL-only | ' + ' | '.join('max_eta %.0e' % e for e in etas)
             + ' |\n|---|---|' + '---|' * len(etas) + '\n')
    for a in ARMS:
        r = res[a]
        row = []
        for e in etas:
            et = X.eta_eos(r['xe'], r['T'], e)
            dx1 = np.broadcast_to(r['g']['dx1'][None, None, None, :], et.shape)
            dxm = np.minimum(np.minimum(dx1**2, r['g']['dx2']**2), r['g']['dx3']**2)
            row.append(min(r['dt_cfl'].min(), (C.CFL * dxm / (6.0 * et)).min()))
        fh.write('| %s | %.2f | ' % (a, r['dt_cfl'].min())
                 + ' | '.join('%.2f' % v for v in row) + ' |\n')

    fh.write('\n### dt vs a uniform temperature boost dT (what more twilight heating '
             'would buy), max_eta fixed at 1e13\n\n')
    boosts = [0, 25, 50, 100, 200, 400]
    fh.write('| arm | ' + ' | '.join('dT = %d K' % b for b in boosts) + ' |\n|---|'
             + '---|' * len(boosts) + '\n')
    for a in ARMS:
        r = res[a]
        row = []
        for b in boosts:
            Tb = r['T'] + b
            et = X.eta_eos(XE(r['d']['dens'], Tb), Tb, MAX_ETA)
            dx1 = np.broadcast_to(r['g']['dx1'][None, None, None, :], et.shape)
            dxm = np.minimum(np.minimum(dx1**2, r['g']['dx2']**2), r['g']['dx3']**2)
            row.append(min(r['dt_cfl'].min(), (C.CFL * dxm / (6.0 * et)).min()))
        fh.write('| %s | ' % a + ' | '.join('%.2f' % v for v in row) + ' |\n')
    fh.close()

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    for n, a in enumerate(ARMS):
        r = res[a]
        ib = ibs[a]
        A = ax.ravel()[n]
        sh = (slice(None), slice(None), slice(None), ib)
        lo = r['lon'][sh].ravel()
        la = r['lat'][sh].ravel()
        val = np.log10(r['eta'][sh].ravel())
        s = A.scatter(lo, la, c=val, s=7, cmap='inferno', vmin=8, vmax=13)
        A.set_title('%s  i=%d  r/Rp %.3f  dt %.1f s'
                    % (a, ib, r['g']['rc'][ib] / RP, r['dt_all'].min()), fontsize=9)
        A.set_xlabel('lon from substellar [deg]', fontsize=8)
        A.set_ylabel('lat [deg]', fontsize=8)
        A.axvline(-90, color='c', lw=0.6)
        A.axvline(90, color='c', lw=0.6)
        plt.colorbar(s, ax=A, label='log10 eta')
    fig.suptitle('log10 eta on the binding radial shell (max_eta = 1e13 is the ceiling)')
    fig.tight_layout()
    fig.savefig(C.OUT + 'fig_ohmic.png', dpi=110)
    print(open(C.OUT + 'ohmic_tables.md').read())
