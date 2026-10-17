#!/usr/bin/env python3
"""gate_aprime.py DECK_RUN PAIR_RUN [options]: the ADOPTED atmosphere gate (a') (USER 2026-10-17).

Replaces the absolute atmosphere thresholds of gate (a) (gate_atm.py: 0.99-1.1 R_ph |<v_r>| < 1 km/s,
|drho| < 0.1) for stars whose photosphere/thin atmosphere carries a converged, bounded, self-excited
oscillation (BSG new setup, bsgres-1017).  DECK_RUN = the gate run on the deck grid, PAIR_RUN = the same
gate on its 2x-atmosphere grid (seed 0, 16x16 lateral, >= 8 d, 0.05-d hydro_w + m1 dumps, hst).

Items (all must PASS; items 3-5 are evaluated on BOTH runs):
 1 RESOLUTION   >= 4.4 cells per atmospheric (density) scale height over 0.95-1.1 R_ph on the deck grid:
                min over the band of H_rho/dr, H_rho = -dr/dln<rho> of the t = 0 dump (the code's startup
                atmosphere), cells at the floor plateau (rho < 1.5 x min rho) excluded.
 2 CONVERGENCE  deck vs pair, saturated state t >= TSAT (4 d), band 0.99-1.1 R_ph, both runs >= 8 d:
                band mean rms <v_r> within 15 %; time-rms <v_r>(r) profile L2 difference <= 10 %
                (pair profile interpolated to the deck radii, ||p_d - p_p||/||p_p||); dominant period
                (FFT, median of the 5 top-rms radii) within 10 %.
 3 BOUNDED      no growth (max band |<v_r>| over the last third <= 1.5 x over the middle third) and
                saturated band mean rms <= c_s,iso(R_ph) (from the deck's he_sponge_rate rule comment
                'c_s ...' or --cs-iso).
 4 STAR         interior 0.5-0.9 R_ph max |<v_r>| < 1 km/s and 0.9-0.99 R_ph < 3 km/s (whole run);
                (0.5-0.9: whole run; 0.9-0.99: saturated state t >= TSAT, as the bsgres-1017
                proposal measured it; the whole-run 0.9-0.99 max is printed for information);
                R_ph(tau_R = 2/3) within +-2 % of its t = 0 value at every dump; L_top/L_in within +-5 %
                for t >= 2 d; max |M_tot/M_tot(0) - 1| < 1e-6; max |Mdot_top| < 1e-8 Msun/yr.
 5 ROBUSTNESS   summary NON-CONVERGED = 0, no FATAL / DIVERGED / NaN in run.log, rc 0 (slurm log
                'rc=N', else run.log 'Terminating on time limit'); dt_min >= 0.25 x dt_median (hst dt, without
                the first row and the final step truncated to tlim);
                dt_retry retries reported.
R_ph = the code's startup sph1d photosphere (run.log 'R_ph A -> B', B) unless --rph.
Shell means: mass-weighted <v_r>, volume-weighted <rho>, <E> (as gate_atm.py / bsgres res_ana.py).
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                                'vis', 'python'))
import bin_convert as bc  # noqa: E402

RSUN = 6.957e10
DAY = 86400.0
C = 2.99792458e10
A_RAD = 7.565733e-15
MSUN_YR = 1.989e33/3.15576e7
LIM = dict(cells_per_h=4.4, rms_rel=0.15, prof_l2=0.10, per_rel=0.10, t_min=8.0, grow=1.5,
           v_int=1.0, v_sub=3.0, rph_rel=0.02, lrat=0.05, lrat_t=2.0, dm=1e-6, mdot=1e-8,
           dt_ratio=0.25)


def deck_of(run, override=None):
    if override:
        return override
    txt = open(os.path.join(run, 'cmd.txt')).read()
    return re.search(r'-i\s+(\S+)', txt).group(1)


def deck_key(deck, key):
    for ln in open(deck):
        m = re.match(r'^\s*%s\s*=\s*([^#\s]+)(.*)$' % key, ln)
        if m:
            return m.group(1), m.group(2)
    return None, ''


def read_table(fn):
    hdr = [ln for ln in open(fn) if ln.startswith('#')]
    nT, nD, lT0, dlT, lD0, dlD = [float(x) for x in hdr[-2].lstrip('#').split()]
    nT, nD = int(nT), int(nD)
    kt = np.loadtxt(fn).reshape(nT, nD)

    def kap(T, rho):
        x = np.clip((np.log10(T) - lT0)/dlT, 0, nT - 1.001)
        y = np.clip((np.log10(rho) - lD0)/dlD, 0, nD - 1.001)
        i, j = x.astype(int), y.astype(int)
        fx, fy = x - i, y - j
        return 10**(kt[i, j]*(1-fx)*(1-fy) + kt[i+1, j]*fx*(1-fy) + kt[i, j+1]*(1-fx)*fy +
                    kt[i+1, j+1]*fx*fy)
    return kap


def period(t, s):
    m = len(s)
    if m < 16:
        return np.nan
    s = s - np.polyval(np.polyfit(t, s, 1), t)
    nfft = 8*m
    fr = np.fft.rfftfreq(nfft, d=np.median(np.diff(t)))
    p = np.abs(np.fft.rfft(s*np.hanning(m), nfft))**2
    ok = fr > 1.0/(t[-1] - t[0])
    k = np.where(ok)[0][np.argmax(p[ok])]
    return float(1.0/fr[k])


def run_log_checks(run, slurm_log):
    lg = open(os.path.join(run, 'run.log'), errors='replace').read()
    nc = re.findall(r'NON-CONVERGED=([0-9.eE+-]+)', lg)
    rc = None
    if slurm_log and os.path.exists(slurm_log):
        m = re.findall(r'\brc=(\d+)', open(slurm_log, errors='replace').read())
        rc = int(m[-1]) if m else None
    rc_src = 'slurm log' if rc is not None else 'run.log'
    if rc is None:
        rc = 0 if 'Terminating on time limit' in lg else -1
    rt = re.findall(r'dt_retry: \d+ checkpoints.*?, (\d+) retries', lg)
    m = re.search(r'R_ph ([0-9.eE+-]+) -> ([0-9.eE+-]+)', lg)
    return dict(nc=float(nc[-1]) if nc else None, fatal=len(re.findall(r'FATAL', lg)),
                diverged=len(re.findall(r'DIVERGED', lg)),
                nan=len(re.findall(r'(?i)(?<![a-z])(nan|-nan|inf)(?![a-z])', lg)), rc=rc,
                rc_src=rc_src, retries=int(rt[-1]) if rt else None,
                rph_start=(float(m.group(1)), float(m.group(2))) if m else None)


def load(run, kap, rph_cli, tsat):
    fh = sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin')))
    fm = sorted(glob.glob(os.path.join(run, 'bin', '*.m1.*.bin')))
    n = min(len(fh), len(fm))
    T, V, D, E = [], [], [], []
    for a, b in zip(fh[:n], fm[:n]):
        h = bc.read_binary_as_athdf(a)
        mm = bc.read_binary_as_athdf(b)
        x2f = np.asarray(h['x2f'], float)
        dom = (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :, None] * \
            np.diff(np.asarray(h['x3f'], float))[:, None, None]
        rho = np.asarray(h['dens'], float)
        w = np.broadcast_to(dom, rho.shape)
        ms = (rho*w).sum(axis=(0, 1))
        T.append(float(h['Time'])/DAY)
        D.append(ms/w.sum(axis=(0, 1)))
        V.append((rho*np.asarray(h['velx'], float)*w).sum(axis=(0, 1))/ms/1e5)
        E.append((np.asarray(mm['m1_e'], float)*w).sum(axis=(0, 1))/w.sum(axis=(0, 1)))
    T, V, D, E = np.array(T), np.array(V), np.array(D), np.array(E)
    r = np.asarray(h['x1v'], float)
    rf = np.asarray(h['x1f'], float)
    return T, r, rf, V, D, E


def analyse(run, deck, slurm_log, rph_cli, cs_cli, tsat):
    lg = run_log_checks(run, slurm_log)
    kap = read_table(deck_key(deck, 'he_opac_table')[0])
    rph = rph_cli or (lg['rph_start'][1] if lg['rph_start'] else None)
    if rph is None:
        raise SystemExit('%s: no startup R_ph in run.log; pass --rph' % run)
    cs = cs_cli
    if cs is None:
        _, com = deck_key(deck, 'he_sponge_rate')
        m = re.search(r'c_s ([0-9.eE+-]+)', com)
        cs = float(m.group(1))/1e5 if m else None
    T, r, rf, V, D, E = load(run, kap, rph, tsat)
    x = r/rph
    dr = np.diff(rf)
    b = (x >= 0.99) & (x < 1.1)
    vmax = np.abs(V[:, b]).max(axis=1)
    vrms = np.sqrt((V[:, b]**2).mean(axis=1))
    n = len(vmax)
    grow = bool(vmax[2*n//3:].max() > LIM['grow']*max(vmax[n//3:2*n//3].max(), 1e-3)) \
        if n >= 6 else None
    sat = T >= tsat
    prof = np.sqrt((V[sat]**2).mean(axis=0)) if sat.any() else np.full(len(r), np.nan)
    ib = np.where(b)[0]
    top5 = ib[np.argsort(prof[ib])[-5:]]
    per = float(np.median([period(T[sat], V[sat, i]) for i in top5])) if sat.sum() >= 16 \
        else np.nan
    # resolution: cells per H_rho at t = 0 (startup atmosphere), floor plateau excluded
    rho0 = D[0]
    fl = rho0.min()
    with np.errstate(divide='ignore', invalid='ignore'):
        H = -(r[1:] - r[:-1])/np.log(rho0[1:]/rho0[:-1])
    xm = 0.5*(x[1:] + x[:-1])
    drm = 0.5*(dr[1:] + dr[:-1])
    ok = (xm >= 0.95) & (xm < 1.1) & (rho0[1:] > 1.5*fl) & (rho0[:-1] > 1.5*fl) & (H > 0) & \
        np.isfinite(H)
    nh = H/drm
    i_nh = np.where(ok)[0][np.argmin(nh[ok])] if ok.any() else None
    # photosphere radius tau_R = 2/3 from the top, every dump
    rphT = []
    for j in range(len(T)):
        Tr = (np.maximum(E[j], 1e-30)/A_RAD)**0.25
        dtau = kap(Tr, np.maximum(D[j], 1e-30))*D[j]*dr
        tau = np.cumsum(dtau[::-1])[::-1]
        k = np.where(tau <= 2/3)[0]
        i0 = max(k[0] if len(k) else len(r) - 1, 1)
        ta, tb = tau[i0 - 1], tau[i0]
        rphT.append(rf[i0 - 1] + (rf[i0] - rf[i0 - 1])*(ta - 2/3)/max(ta - tb, 1e-30))
    rphT = np.array(rphT)
    # hst
    fu = glob.glob(os.path.join(run, '*.user.hst'))[0]
    cols = {}
    for ln in open(fu):
        if ln.startswith('#') and '[1]=' in ln:
            cols = {nm: int(i) - 1 for i, nm in re.findall(r'\[(\d+)\]=(\S+)', ln)}
            break
    u = np.loadtxt(fu)
    th = u[:, cols['time']]/DAY
    lr = u[:, cols['L_top']]/u[:, cols['L_in']]
    late = th >= LIM['lrat_t']
    dm = np.abs(u[:, cols['M_tot']]/u[0, cols['M_tot']] - 1).max()
    mdot = np.abs(u[:, cols['Mdot_top']]).max()/MSUN_YR
    # dt: rows 1 .. n-2 (row 0 = the startup dt, the last row = the step truncated to tlim)
    dtc = u[1:-1, cols['dt']] if len(u) > 3 else u[:, cols['dt']]
    bands = {}
    for lo, hi in [(0.5, 0.9), (0.9, 0.99), (0.99, 1.1)]:
        k = (x >= lo) & (x < hi)
        bands['%.2f-%.2f' % (lo, hi)] = float(np.abs(V[:, k]).max()) if k.any() else np.nan
        bands['%.2f-%.2f sat' % (lo, hi)] = float(np.abs(V[sat][:, k]).max()) \
            if k.any() and sat.any() else np.nan
    return dict(run=run, deck=deck, nx1=len(r), ndump=len(T), t_end=float(T[-1]), rph=rph,
                cs_iso=cs, x=x, prof=prof, sat_vrms_mean=float(vrms[sat].mean()) if sat.any()
                else np.nan, sat_vmax=float(vmax[sat].max()) if sat.any() else np.nan,
                vmax=float(vmax.max()), growth=grow, period=per,
                cells_per_h=float(nh[i_nh]) if i_nh is not None else np.nan,
                cells_per_h_x=float(xm[i_nh]) if i_nh is not None else np.nan,
                h_rho_rsun=float(H[i_nh]/RSUN) if i_nh is not None else np.nan,
                dr_rsun=float(drm[i_nh]/RSUN) if i_nh is not None else np.nan,
                bands=bands, rphT0=float(rphT[0]), rphT_dev=float(np.abs(rphT/rphT[0] - 1).max()),
                rphT_min=float(rphT.min()/RSUN), rphT_max=float(rphT.max()/RSUN),
                lrat_min=float(lr[late].min()) if late.any() else np.nan,
                lrat_max=float(lr[late].max()) if late.any() else np.nan, dM=float(dm),
                mdot_msun_yr=float(mdot), dt_min=float(dtc.min()), dt_med=float(np.median(dtc)),
                log=lg)


def verdict(d, p):
    it = []

    def add(item, name, ok, val):
        it.append((item, name, bool(ok), val))
    # 1 resolution (deck grid; the pair's value reported)
    add(1, 'cells per H_rho (deck grid, min over 0.95-1.1 R_ph) >= %.1f' % LIM['cells_per_h'],
        d['cells_per_h'] >= LIM['cells_per_h'],
        '%.2f at %.4f R_ph (H_rho %.4f Rsun, dr %.4f); pair %.2f' %
        (d['cells_per_h'], d['cells_per_h_x'], d['h_rho_rsun'], d['dr_rsun'], p['cells_per_h']))
    # 2 convergence
    add(2, 'both runs >= %.0f d' % LIM['t_min'], min(d['t_end'], p['t_end']) >= LIM['t_min'] - 1e-3,
        '%.2f / %.2f d' % (d['t_end'], p['t_end']))
    rr = abs(d['sat_vrms_mean']/p['sat_vrms_mean'] - 1)
    add(2, 'band mean rms within %.0f %%' % (100*LIM['rms_rel']), rr <= LIM['rms_rel'],
        '%.3f vs %.3f km/s (%.1f %%)' % (d['sat_vrms_mean'], p['sat_vrms_mean'], 100*rr))
    b = (d['x'] >= 0.99) & (d['x'] < 1.1)
    pp = np.interp(d['x'][b], p['x'], p['prof'])
    l2 = float(np.sqrt(((d['prof'][b] - pp)**2).sum()/(pp**2).sum()))
    add(2, 'rms profile L2 diff <= %.0f %%' % (100*LIM['prof_l2']), l2 <= LIM['prof_l2'],
        '%.1f %%' % (100*l2))
    pr = abs(d['period']/p['period'] - 1)
    add(2, 'period within %.0f %%' % (100*LIM['per_rel']), pr <= LIM['per_rel'],
        '%.3f vs %.3f d (%.1f %%)' % (d['period'], p['period'], 100*pr))
    for tag, q in (('deck', d), ('pair', p)):
        add(3, '%s: no growth' % tag, q['growth'] is False, 'growth %s, max %.2f km/s, sat max %.2f'
            % (q['growth'], q['vmax'], q['sat_vmax']))
        add(3, '%s: sat band rms <= c_s,iso(R_ph)' % tag,
            q['cs_iso'] is not None and q['sat_vrms_mean'] <= q['cs_iso'],
            '%.2f vs %s km/s' % (q['sat_vrms_mean'], '%.2f' % q['cs_iso'] if q['cs_iso'] else '?'))
        add(4, '%s: 0.5-0.9 R_ph max|<v_r>| < %g km/s' % (tag, LIM['v_int']),
            q['bands']['0.50-0.90'] < LIM['v_int'], '%.3f' % q['bands']['0.50-0.90'])
        add(4, '%s: 0.9-0.99 R_ph sat max|<v_r>| < %g km/s' % (tag, LIM['v_sub']),
            q['bands']['0.90-0.99 sat'] < LIM['v_sub'], '%.3f (whole run %.3f)' %
            (q['bands']['0.90-0.99 sat'], q['bands']['0.90-0.99']))
        add(4, '%s: R_ph(tau 2/3) within +-%.0f %%' % (tag, 100*LIM['rph_rel']),
            q['rphT_dev'] <= LIM['rph_rel'], 't0 %.3f Rsun, range %.3f-%.3f, max dev %.2f %%' %
            (q['rphT0']/RSUN, q['rphT_min'], q['rphT_max'], 100*q['rphT_dev']))
        add(4, '%s: L_top/L_in within +-%.0f %% (t >= %g d)' % (tag, 100*LIM['lrat'], LIM['lrat_t']),
            abs(q['lrat_min'] - 1) <= LIM['lrat'] and abs(q['lrat_max'] - 1) <= LIM['lrat'],
            '%.4f-%.4f' % (q['lrat_min'], q['lrat_max']))
        add(4, '%s: |dM/M| < %g' % (tag, LIM['dm']), q['dM'] < LIM['dm'], '%.2e' % q['dM'])
        add(4, '%s: Mdot_top < %g Msun/yr' % (tag, LIM['mdot']), q['mdot_msun_yr'] < LIM['mdot'],
            '%.2e' % q['mdot_msun_yr'])
        lg = q['log']
        add(5, '%s: NON-CONVERGED = 0' % tag, lg['nc'] == 0, '%s' % lg['nc'])
        add(5, '%s: no FATAL/DIVERGED/NaN' % tag, lg['fatal'] + lg['diverged'] + lg['nan'] == 0,
            '%d/%d/%d' % (lg['fatal'], lg['diverged'], lg['nan']))
        add(5, '%s: rc 0' % tag, lg['rc'] == 0, '%d (%s)' % (lg['rc'], lg['rc_src']))
        add(5, '%s: dt_min >= %.2f dt_median' % (tag, LIM['dt_ratio']),
            q['dt_min'] >= LIM['dt_ratio']*q['dt_med'], '%.4g / %.4g = %.2f; dt_retry retries %s' %
            (q['dt_min'], q['dt_med'], q['dt_min']/q['dt_med'], lg['retries']))
    return it, l2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('deck_run')
    ap.add_argument('pair_run')
    ap.add_argument('--deck-in', help='athinput of the deck run (default: -i of DECK_RUN/cmd.txt)')
    ap.add_argument('--pair-in', help='athinput of the pair run (default: -i of PAIR_RUN/cmd.txt)')
    ap.add_argument('--deck-log', help='slurm log of the deck run (rc=N)')
    ap.add_argument('--pair-log', help='slurm log of the pair run (rc=N)')
    ap.add_argument('--rph', type=float, help='R_ph [cm] (default: startup sph1d R_ph of each run)')
    ap.add_argument('--cs-iso', type=float, help='c_s,iso(R_ph) [km/s] (default: deck comment)')
    ap.add_argument('--tsat', type=float, default=4.0, help='saturated state from t [d] (4)')
    ap.add_argument('--json', help='write the results here')
    a = ap.parse_args()
    d = analyse(a.deck_run, deck_of(a.deck_run, a.deck_in), a.deck_log, a.rph, a.cs_iso, a.tsat)
    p = analyse(a.pair_run, deck_of(a.pair_run, a.pair_in), a.pair_log, a.rph, a.cs_iso, a.tsat)
    it, l2 = verdict(d, p)
    print("# gate (a') (ADOPTED 2026-10-17): deck %s (nx1 %d) vs pair %s (nx1 %d); R_ph %.6e / %.6e"
          % (a.deck_run, d['nx1'], a.pair_run, p['nx1'], d['rph'], p['rph']))
    allok = True
    for item, name, ok, val in it:
        allok &= ok
        print('  [%d] %-58s %-4s %s' % (item, name, 'PASS' if ok else 'FAIL', val))
    fails = sorted(set(i for i, _, ok, _ in it if not ok))
    print("# gate (a') %s%s" % ('PASS' if allok else 'FAIL',
                                '' if allok else ' (items %s)' % ','.join(map(str, fails))))
    if a.json:
        def strip(q):
            return {k: v for k, v in q.items() if k not in ('x', 'prof')}
        json.dump(dict(deck=strip(d), pair=strip(p), profile_l2=l2, limits=LIM,
                       items=[dict(item=i, name=nm, ok=ok, value=v) for i, nm, ok, v in it],
                       verdict='PASS' if allok else 'FAIL', failed_items=fails),
                  open(a.json, 'w'), indent=1, default=float)
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
