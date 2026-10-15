#!/usr/bin/env python3
"""make_star_ic.py SPEC [--out DIR] [--compare REF] [--no-png]

THE standard initial condition for he_star_m1 stellar-envelope runs (BSG, He giants, AG Car /
LBVs, future stars): spec (YAML/JSON) -> tables (cached) -> IC column + provenance header ->
deck block -> diagnostics (text + PNG).  The method is the AG Car one (make_ic_agcar_ge.py):
see README.md and icstd/column.py.

Outputs in DIR (default spec 'out_dir', else ./<name>):
  ic_<name>.txt          the column (r rho eint F_r E T F_MLT/F), provenance header
  deck_<name>.athinput   the deck block (star/IC/background keys)
  diag_<name>.txt/.png   diagnostics
  prov_<name>.json       resolved inputs, md5s, git hash
--compare REF: max relative difference per column vs a reference column (rows matched by
  r; the regression check).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from icstd import analyze as AN  # noqa: E402
from icstd import column as COL  # noqa: E402
from icstd import deck as DK  # noqa: E402
from icstd import mesa as MS  # noqa: E402
from icstd import mesacol as MC  # noqa: E402
from icstd import spec as SP  # noqa: E402
from icstd import tables as TB  # noqa: E402
from icstd.micro import RSUN, LSUN, Star, SIG, q_hopf  # noqa: E402

COLS = ['r', 'rho', 'eint', 'F_r', 'E', 'T', 'F_MLT/F']


def git_info():
    try:
        h = subprocess.run(['git', '-C', HERE, 'rev-parse', 'HEAD'], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True).stdout.strip()
        d = subprocess.run(['git', '-C', HERE, 'status', '--porcelain', '--', '.'],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout
        return h + ('+dirty' if d.strip() else '')
    except Exception:                                    # noqa: BLE001
        return 'unknown'


def scripts_md5():
    h = hashlib.md5()
    fs = [os.path.join(HERE, 'make_star_ic.py')] + sorted(
        os.path.join(HERE, 'icstd', f) for f in os.listdir(os.path.join(HERE, 'icstd'))
        if f.endswith('.py'))
    for f in fs:
        h.update(open(f, 'rb').read())
    return h.hexdigest()


def compare(fn_new, fn_ref):
    a = np.loadtxt(fn_new, comments='#')
    b = np.loadtxt(fn_ref, comments='#')
    out = ['regression: %s vs %s' % (fn_new, fn_ref),
           '  rows new %d ref %d' % (a.shape[0], b.shape[0])]
    if a.shape == b.shape:
        dr = np.max(np.abs(a[:, 0]/b[:, 0] - 1))
        out.append('  same row count; max |r_new/r_ref - 1| = %.3e' % dr)
        res = {}
        for j, nm in enumerate(COLS):
            x, y = a[:, j], b[:, j]
            nz = np.abs(y) > 0
            rel = np.max(np.abs(x[nz]/y[nz] - 1)) if nz.any() else 0.0
            ab = np.max(np.abs(x[~nz] - y[~nz])) if (~nz).any() else 0.0
            res[nm] = rel
            out.append('  %-8s max rel diff %.3e, max abs diff %.3e%s'
                       % (nm, rel, np.max(np.abs(x - y)),
                          '' if ab == 0 else ' ; abs diff where ref = 0: %.3e' % ab))
        la = [ln for ln in open(fn_new) if not ln.startswith('#')]
        lb = [ln for ln in open(fn_ref) if not ln.startswith('#')]
        same = sum(1 for p, q in zip(la, lb) if p == q)
        out.append('  text-identical data rows: %d of %d' % (same, len(lb)))
        out.append('  MAX REL DIFF over all columns: %.3e' % max(res.values()))
    else:
        rr = b[:, 0]
        m = (rr >= a[0, 0]) & (rr <= a[-1, 0])
        for j, nm in enumerate(COLS[1:], 1):
            yi = np.exp(np.interp(np.log(rr[m]), np.log(a[:, 0]),
                                  np.log(np.maximum(np.abs(a[:, j]), 1e-300))))
            y = np.abs(b[m, j])
            nz = y > 0
            out.append('  %-8s (interpolated) max rel diff %.3e'
                       % (nm, np.max(np.abs(yi[nz]/y[nz] - 1))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('spec')
    ap.add_argument('--out')
    ap.add_argument('--compare')
    ap.add_argument('--no-png', action='store_true')
    ap.add_argument('--no-scan', action='store_true', help="mode 'global': skip the alpha / "
                    "MLT++ branch scan")
    a = ap.parse_args()
    t0 = time.time()
    s = SP.load(a.spec)
    name = s['name']
    out = a.out or s.get('out_dir') or os.path.join(os.getcwd(), name)
    os.makedirs(out, exist_ok=True)
    rep = []

    def log(x):
        print(x, flush=True)
        rep.append(x)
    log('# make_star_ic.py %s (git %s, scripts md5 %s)' % (a.spec, git_info(), scripts_md5()))
    cache = s['_cfg']['cache']
    comp = s['composition']
    if comp.get('Y') is None:
        comp['Y'] = 1.0 - comp['X'] - comp['Z']
    # ---------------------------------------------------------------- star
    mesa = None
    if s.get('mesa_profile'):
        mp = s['mesa_profile']
        mesa = MS.read(mp['path'], mp.get('format', 'mesa'), mp.get('r_max_rsun'))
        log('MESA profile %s (md5 %s): %s' % (mp['path'], TB.md5(mp['path']), mesa['note']))
    mode = s.get('mode') or ('mesa' if mesa is not None else 'global')
    if mode == 'mesa' and mesa is None:
        raise SystemExit("mode 'mesa' needs mesa_profile")
    log('MODE %s (%s)' % (mode, 'MESA T(P) structure, our EOS/kappa/gravity' if mode == 'mesa'
                          else 'surface-inward MLT/MLT++ column (AG Car)'))
    sv, ssrc = SP.resolve_star(s, mesa)
    L, GM, R = sv['L'], sv['GM'], sv['R']
    Tph = (3.0*(L/(4*np.pi*R**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
    log('star %s: L %.8e erg/s (log L/Lsun %.4f; %s), GM %.8e (%s), R_ph %.8e cm = %.4f Rsun '
        '(%s), Teff %.2f K, T_ph(Hopf, tau 2/3) %.2f K' % (name, L, np.log10(L/LSUN), ssrc['L'],
                                                         GM, ssrc['GM'], R, R/RSUN, ssrc['R'],
                                                         sv['Teff'], Tph))
    rt = s['r_top']
    x1max = float(rt['x1max_cm']) if rt.get('x1max_cm') else rt['rtop_fac']*R
    if rt.get('col_top_cm'):
        col_top = float(rt['col_top_cm'])
    elif rt.get('x1max_cm'):
        col_top = rt['col_top_fac']*x1max
    else:                               # (make_ic_agcar_ge.py op order: 1.15*rout_fac*R)
        col_top = rt['col_top_fac']*rt['rtop_fac']*R
    # tfloor rule first (analytic): it sets the EOS range
    trad = AN.t_rad_at(type('o', (), dict(RSTAR=R, LUM=L))(), x1max, Tph)
    bg = s['background']
    tfloor = float(bg['tfloor']) if bg.get('tfloor') else 100.0*np.floor(
        bg['tfloor_rel']*trad/100.0)
    e = s['eos']
    if os.environ.get('ICSTD_EOS_RANGE'):           # the automatic re-run with a wider table
        e.update(json.loads(os.environ['ICSTD_EOS_RANGE']))
        log('EOS range extended (automatic re-run): %s' % os.environ['ICSTD_EOS_RANGE'])
    elif mode == 'mesa' and not e.get('dump'):       # cover MESA's base with margin up front
        rb_ = SP.rsun_or_cm(s['r_in'], 'bottom') or SP.rsun_or_cm(s['r_in'], 'r')
        if rb_ is None:
            rb_ = mesa['rc'][0]
        kk = int(np.argmin(abs(mesa['rc'] - rb_)))
        lt_b = np.log10(mesa['T'][:kk + 1].max())
        ld_b = np.log10(mesa['rho'][:kk + 1].max())
        e['logt_max'] = max(e['logt_max'], np.ceil(10*(lt_b + 0.1))/10)
        e['logd_max'] = max(e['logd_max'], np.ceil(ld_b + 0.5))
    if e.get('logt_min') is None:
        e['logt_min'] = float(np.floor(10*np.log10(tfloor) + 1e-9)/10.0)
    # ---------------------------------------------------------------- tables
    log('tables:')
    ot = TB.opacity_tables(s, cache, log)
    ed = TB.eos_dump(s, cache, log)
    log('  Rosseland %s md5 %s' % (ot['ross'], TB.md5(ot['ross'])))
    if ot.get('planck'):
        log('  Planck    %s md5 %s' % (ot['planck'], TB.md5(ot['planck'])))
    log('  EOS dump  %s md5 %s (%s)' % (ed['dump'], TB.md5(ed['dump']), ed['note']))
    ml = s['mlt']
    ross_ic = s['opacity'].get('ic_rosseland') or ot['ross']
    if ross_ic != ot['ross']:
        log('  IC built with kappa_R %s md5 %s (opacity.ic_rosseland; the deck names the '
            'tables above)' % (ross_ic, TB.md5(ross_ic)))
    legacy = ml['mltpp'] == 'agcar_legacy'
    if ml['mltpp'] not in ('agcar_legacy', 'mesa_gated'):
        raise SystemExit("mlt.mltpp must be 'mesa_gated' or 'agcar_legacy'")
    st = Star(GM, L, R, ed['dump'], ross_ic, s['opacity']['ld_ext'], ml['alpha'],
              ml['pp_alpha'] if legacy else 0.0, ml['f1'], ml['f2'],
              ml['pp_form'] if legacy else 'mesa', verbose=False,
              eos_clamp=(s['grid']['scheme'] != 'mltpp_inward'),
              mlt_option=ml['mlt_option'], alt_hp=ml['alt_scale_height'])
    st.PP_MAXLOGT = ml['max_logT']
    log('MLT: alpha %g, %s, H_p %s; MLT++ %s' % (
        ml['alpha'], ml['mlt_option'], 'min(P/(rho g), sqrt(P/G)/rho) (MESA alt_scale_height)'
        if ml['alt_scale_height'] else 'P/(rho g)',
        ('agcar_legacy: fixed pp_alpha %g (%s form) where nabla - nabla_ad > %g'
         % (ml['pp_alpha'], ml['pp_form'], ml['f1'])) if legacy else
        ('MESA-gated (use_mltpp %s; f1 %g f2 %g; lambda1/beta1 %g/%g lambda2/beta2 %g/%g '
         'dlambda/dbeta %g/%g; xa1_center %g)' % (ml['use_mltpp'], ml['f1'], ml['f2'],
                                                  ml['lambda1'], ml['beta1'], ml['lambda2'],
                                                  ml['beta2'], ml['dlambda'], ml['dbeta'],
                                                  ml['xa1_center']))))
    if st.eos.comp:
        log('  EOS dump composition (X Y Z H2 ion metal_ion mh): %s; spec X %g Y %g Z %g'
            % (' '.join(st.eos.comp), comp['X'], comp['Y'], comp['Z']))
    gr = s['grid']
    c = dict(rtol=gr['rtol'], atm_rad=gr['atm_rad'], rho_lo=gr['rho_lo'], rho_hi=gr['rho_hi'],
             tauj=gr['tauj'], dlnr=gr['dlnr'], drn_atm=gr['drn_atm_rph']*R,
             fl_step=gr['fl_step'], r_col_top=col_top, scheme=gr['scheme'],
             drn_atm_rph=gr['drn_atm_rph'], rho_floor_rel=bg['rho_floor_rel'])
    # ---------------------------------------------------------------- column floor
    if bg.get('rho_floor'):
        c['rho_floor'] = float(bg['rho_floor'])
        log('column floor %.3e (spec)' % c['rho_floor'])
    elif mode == 'mesa':
        c['rho_floor'] = None
        log('column floor: rule decade(%g rho_ph) applied after the build (our rho_ph)'
            % bg['rho_floor_rel'])
    else:
        rph0, _ = AN.iso_rho_ph(st, c)
        c['rho_floor'] = AN.decade(bg['rho_floor_rel']*rph0)
        log('column floor %.1e = decade(%g x rho_ph %.4e) (rule)' % (c['rho_floor'],
                                                                   bg['rho_floor_rel'], rph0))
    # ---------------------------------------------------------------- r_in rule
    ri = s['r_in']
    if mode == 'global' and not legacy and ml['use_mltpp']:
        log('MLT++ gate (MESA set_gradT_excess_alpha fixed point, coarse columns to T %.1e K '
            'or the column bottom):' % ri['explore_tmax'])
        rstop = SP.rsun_or_cm(ri, 'bottom') or 0.01*R
        ctrl = {k: ml[k] for k in ('lambda1', 'beta1', 'lambda2', 'beta2', 'dlambda', 'dbeta')}
        a_ = COL.solve_gated(st, dict(c, rho_floor=c.get('rho_floor') or 1e-16), rstop,
                             ri['explore_tmax'], ml['xa1_center'], ctrl, ml['pp_tol'], log)
        log('  -> gradT_excess_alpha* = %.5f' % a_)
    ms = MC.mesa_struct(st, mesa) if mesa is not None else None
    if ri['rule'] == 'explicit':
        x1min = SP.rsun_or_cm(ri, 'r')
        rin_rule = 'explicit'
    elif ri['rule'] == 'hp_below_cz' and mode == 'mesa':
        x1min, fzm = MC.mesa_rin(ms, ri['n_hp'], ri['cz_fmlt'])
        if x1min is None:
            raise SystemExit('r_in rule: the MESA profile is too shallow for %g H_p' % ri['n_hp'])
        rin_rule = ('%g H_p below MESA\'s Fe-CZ base %.3f Rsun (Fe peak %.3f Rsun, %s in MESA)'
                    % (ri['n_hp'], fzm['r_base']/RSUN, fzm['r_fe']/RSUN,
                       'convective' if fzm['fe_convective'] else 'RADIATIVE'))
        log('  r_in rule (MESA structure): ' + rin_rule + ' -> x1min %.4f Rsun' % (x1min/RSUN))
    elif ri['rule'] == 'hp_below_cz':
        log('r_in rule: exploring inward to T %.2e K' % ri['explore_tmax'])
        ex = COL.explore(st, dict(c, legacy_deep=False), 0.01*R, ri['explore_tmax'])
        prx = AN.profile(st, ex)
        fz = AN.fe_cz(ex, prx, ri['cz_fmlt'])
        if fz is None:
            raise SystemExit('r_in rule: no Fe opacity peak (1e5 < T < 1e6 K) in the column')
        x1min = AN.r_at_hp_below(ex['r'], prx['Hp'], fz['r_base'], ri['n_hp'])
        if x1min is None:
            raise SystemExit('r_in rule: column too shallow for %g H_p (raise explore_tmax)'
                             % ri['n_hp'])
        rin_rule = '%g H_p below the Fe-CZ base %.3f Rsun' % (ri['n_hp'], fz['r_base']/RSUN)
        log('  Fe peak %.3f Rsun (convective %s), CZ base %.3f Rsun -> x1min %.4f Rsun'
            % (fz['r_fe']/RSUN, fz['fe_convective'], fz['r_base']/RSUN, x1min/RSUN))
    else:
        raise SystemExit('r_in.rule must be explicit or hp_below_cz')
    rbot = SP.rsun_or_cm(ri, 'bottom') or ri['rlo_fac']*x1min
    c.update(r_wall=x1min, r_bot=rbot)
    if ri.get('legacy_deep'):
        ld = ri['legacy_deep']
        rw1 = ld['r_wall1_rsun']*RSUN
        c.update(legacy_deep=True, r_wall1=rw1, r_bot1=ld['rlo_fac1']*rw1)
    log('x1min %.10e (%.4f Rsun, %s); column bottom %.10e (%.4f Rsun); x1max %.10e '
        '(%.4f R_ph); column top %.6e' % (x1min, x1min/RSUN, rin_rule, rbot, rbot/RSUN, x1max,
                                         x1max/R, col_top))
    # ---------------------------------------------------------------- build
    t1 = time.time()
    if mode == 'mesa':
        c['R_target'] = R
        col = MC.column_mesa(st, c, mesa, ms, log)
        c['rho_floor'] = col['rho_floor']
        st.GM, st.RSTAR = col['GM_wall'], col['rph']
        GM = col['GM_wall']
        ssrc['GM'] = 'G m(x1min) of the column (he_gm_column = true)'
        log('column floor %.1e = decade(%g x rho_ph %.4e) (rule)' % (col['rho_floor'],
                                                                   bg['rho_floor_rel'],
                                                                   col['rho_R']))
    else:
        col = COL.build(st, c)
    col['L'] = L
    check_tables(log, s, st, col, ed, a)
    log('column (%s): %d rows, %.0f s; rho_ph %.6e, tau-%g join %.6e (%.5f R_ph), floor from '
        '%.6e (%.4f R_ph)' % (col['scheme'], col['r'].size, time.time() - t1, col['rho_R'],
                              c['tauj'], col['rj'], col['rj']/R, col['r_floor'],
                              col['r_floor']/R))
    # ---------------------------------------------------------------- analysis
    pr = AN.profile(st, col)
    if mode == 'mesa':             # the column's own gradient (MESA T(P)), not local MLT
        ii = col['interior']
        n_ = ii['r'].size
        pr['nab'][:n_], pr['nrad'][:n_] = ii['nab'], ii['nrad']
    fz = AN.fe_cz(col, pr, ri['cz_fmlt'])
    bgv, why = AN.background(st, s, col['tph'], col['rho_R'], x1max, c['rho_floor'])
    F_rin = AN.flux_at(col, x1min)
    r = col['r']
    # ---------------------------------------------------------------- column file
    fn = os.path.join(out, 'ic_%s.txt' % name)
    spec_line = json.dumps({k: v for k, v in s.items() if not k.startswith('_')},
                           sort_keys=True, default=str)
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC (tools/ic_star make_star_ic.py, mode %s, scheme '
                 '%s; he_ic_cols = 5): %s; L %.8e GM %.8e R_ph %.8e Teff %.2f; alpha_MLT %g, '
                 'MLT++ pp_alpha %g f1 %g f2 %g; x1min %.10e; column floor %.1e\n'
                 % (mode, col['scheme'], name, L, GM, col['rph'],
                    (L/(4*np.pi*col['rph']**2*SIG))**0.25, ml['alpha'], ml['pp_alpha'],
                    ml['f1'], ml['f2'], x1min, c['rho_floor']))
        fh.write('# provenance: git %s; scripts md5 %s; EOS dump %s md5 %s; kappa_R %s md5 %s; '
                 'spec %s md5 %s\n' % (git_info(), scripts_md5(), ed['dump'],
                                       TB.md5(ed['dump']), ross_ic, TB.md5(ross_ic),
                                       s['_file'], TB.md5(s['_file'])))
        fh.write('# spec: %s\n' % spec_line)
        fh.write('# r[cm] rho eint F_r(radiative part) E  | T[K] F_MLT/F (cols 6,7)\n')
        for row in zip(r, col['rho'], col['eint'], col['Fr'], col['E'], col['T'],
                       col['fmlt']):
            fh.write(COL.fmt_row(row))
    icmd5 = TB.md5(fn)
    log('wrote %s (md5 %s)' % (fn, icmd5))
    # ---------------------------------------------------------------- deck block
    tt = fz['t_turn'] if fz else float('nan')
    rp = s['ramp']
    dv = dict(name=name, spec=s['_file'], x1min=x1min, x1max=x1max, rtop_fac=x1max/R,
              rin_rule=rin_rule, X=comp['X'], Y=comp['Y'], Z=comp['Z'],
              a_metal=comp.get('a_metal', 16.0), eos_logd_max=e['logd_max'],
              eos_logt_max=e['logt_max'], eos_dlog=e['dlog'], eos_logt_min=e['logt_min'],
              Tmax=float(col['T'].max()), F_rin=F_rin, F_L_rin=L/(4*np.pi*x1min**2), GM=GM,
              gm_src=ssrc['GM'], ic_file=fn, ic_md5=icmd5, ross=ot['ross'],
              ross_md5=TB.md5(ot['ross']), planck=ot.get('planck') or 'NONE',
              planck_md5=TB.md5(ot['planck']) if ot.get('planck') else '-',
              R=col['rph'], rj=col['rj'], gm_column=(mode == 'mesa'), mode=mode, t_turn=tt, start_turnovers=rp['start_turnovers'],
              time_turnovers=rp['time_turnovers'], ramp_start=rp['start_turnovers']*tt,
              ramp_time=rp['time_turnovers']*tt,
              seed_rlo=1.007*(fz['r_base'] if fz else x1min), seed_rhi=0.9786*R,
              why_tfloor=why['tfloor_kelvin'], why_pfloor=why['pfloor'],
              why_vceil=why['vceil'], why_sponge=why['he_sponge_rate'],
              wb_rmax_rel=bg['wb_rmax_rel'], dfloor_rel=bg['dfloor_rel'],
              dmax_rel=bg['dmax_rel'], sponge_r0_rel=bg['sponge_r0_rel'],
              balance_rmax_rel=bg['balance_rmax_rel'], **{k: v for k, v in bgv.items()
                                                          if k != 'eos_logt_min'})
    dfn = os.path.join(out, 'deck_%s.athinput' % name)
    open(dfn, 'w').write(DK.block(dv))
    log('wrote %s' % dfn)
    # ---------------------------------------------------------------- diagnostics
    diag(log, s, st, col, pr, fz, bgv, ot, ed, mesa, x1min, x1max, rbot, F_rin, Tph, tt, dv)
    if mode == 'mesa':
        diag_mesa(log, s, st, col, ms, mesa, x1min, ri['cz_fmlt'])
    else:
        diag_global(log, s, st, col, pr, c, x1min, a.no_scan)
    if a.compare:
        for x in compare(fn, a.compare):
            log(x)
    if not a.no_png:
        from icstd import diag as DG
        info = dict(name=name, x1min=x1min, x1max=x1max, dmax=bgv['he_sponge_dmax'],
                    rho_floor=c['rho_floor'], r_fe=fz['r_fe'] if fz else None,
                    r_base=fz['r_base'] if fz else None, r_top=fz['r_top'] if fz else None)
        pfn = os.path.join(out, 'diag_%s.png' % name)
        DG.png(pfn, col, pr, info, mesa)
        log('wrote %s' % pfn)
    log('total %.0f s' % (time.time() - t0))
    open(os.path.join(out, 'diag_%s.txt' % name), 'w').write('\n'.join(rep) + '\n')
    json.dump(dict(spec=s, star=sv, star_src=ssrc, x1min=x1min, x1max=x1max, r_bot=rbot,
                   F_rin=F_rin, ic=fn, ic_md5=icmd5, tables=dict(opacity=ot, eos=ed),
                   background=bgv, t_turn=tt, git=git_info(), scripts_md5=scripts_md5()),
              open(os.path.join(out, 'prov_%s.json' % name), 'w'), indent=1, default=str)


def diag(log, s, st, col, pr, fz, bgv, ot, ed, mesa, x1min, x1max, rbot, F_rin, Tph, tt, dv):
    r, T, rho, R = col['r'], col['T'], col['rho'], col['rph']
    L = st.LUM
    log('---- coverage')
    lt, lr = np.log10(T), np.log10(rho)
    eos = st.eos
    log('  column log T %.4f..%.4f, log rho %.4f..%.4f' % (lt.min(), lt.max(), lr.min(),
                                                           lr.max()))
    log('  EOS dump grid log T %.3f..%.3f, log rho %.3f..%.3f: column %s; builder clamped '
        'calls %d of %d' % (eos.yr[0], eos.yr[1], eos.xr[0], eos.xr[1],
                            'INSIDE' if (lt.min() >= eos.yr[0] and lt.max() <= eos.yr[1] and
                                         lr.min() >= eos.xr[0] and lr.max() <= eos.xr[1])
                            else 'OUTSIDE (clamped!)', eos.nclamp, eos.ncall))
    e = s['eos']
    ok = (lt.max() <= e['logt_max'] and lr.min() >= e['logd_min'] and
          lr.max() <= e['logd_max'] and lt.min() >= e['logt_min'])
    log('  deck EOS range (eos_logt %.2f..%.2f, eos_logd %.1f..%.1f): column %s (he_star_m1 '
        'fatals outside)' % (e['logt_min'], e['logt_max'], e['logd_min'], e['logd_max'],
                             'inside' if ok else 'OUTSIDE'))
    t0, t1, d0, d1 = st.ross_lims
    log('  Rosseland grid (with he_opac_logd_min %g) log T %.3f..%.3f, log rho %.3f..%.3f: '
        'column %s' % (s['opacity']['ld_ext'], t0, t1, d0, d1,
                       'inside' if (lt.min() >= t0 and lt.max() <= t1 and lr.min() >= d0 and
                                    lr.max() <= d1) else 'OUTSIDE (edge-clamped)'))
    if ot.get('fields') and os.path.exists(ot['fields']):
        f = np.load(ot['fields'])
        fl = np.interp(lt, f['tops_lt'], f['tops_floor'])
        lR = lr - 3*lt + 18.0
        reg = np.where(lt >= 4.2, np.where(lr >= fl, 'TOPS data', 'TOPS extension'),
                       np.where(lt < 4.0, np.where((lR >= f['ferg_lr'][0]) &
                                                   (lR <= f['ferg_lr'][1]) &
                                                   (lt >= f['ferg_lt'][0]),
                                                   'Ferguson data', 'Ferguson ext/edge'),
                                'blend 4.0-4.2'))
        for nm in ('TOPS data', 'TOPS extension', 'blend 4.0-4.2', 'Ferguson data',
                   'Ferguson ext/edge'):
            mm = reg == nm
            if mm.any():
                log('  opacity source %-18s %6d rows, r/R_ph %.4f..%.4f' % (nm, mm.sum(),
                                                                         r[mm].min()/R,
                                                                         r[mm].max()/R))
        below = pr['below']
        log('  below R_ph: rows on model extension (not data): %d'
            % int(np.sum(below & np.isin(reg, ['TOPS extension', 'Ferguson ext/edge']))))
    else:
        log('  opacity: explicit table, source map n/a (%s)' % ot['note'])
    log('---- structure (rows below R_ph)')
    log('  r/R_ph   r[Rsun]     T[K]       rho       kappa  nabla   n_ad   n_rad  FMLT/F '
        'Gam_L  Gam_rad  Hp[Rsun]  tau')
    for xx in list(np.linspace(r[0]/R, 1.0, 14)):
        i = int(np.argmin(abs(r/R - xx)))
        log('  %6.4f %9.3f %10.4e %10.4e %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %8.3f %9.3e'
            % (r[i]/R, r[i]/RSUN, T[i], rho[i], pr['kap'][i], pr['nab'][i], pr['nad'][i],
               pr['nrad'][i], col['fmlt'][i], pr['gam_L'][i], pr['gam_rad'][i],
               pr['Hp'][i]/RSUN, pr['tau'][i]))
    log('---- convection and the wall')
    if fz:
        log('  Fe opacity peak: r %.4f Rsun (%.4f R_ph) T %.4e rho %.4e kappa %.4f, %s'
            % (fz['r_fe']/RSUN, fz['r_fe']/R, fz['T_fe'], fz['rho_fe'], fz['kap_fe'],
               'convective' if fz['fe_convective'] else 'RADIATIVE (no Fe CZ)'))
        log('  Fe CZ (F_MLT/F > %g): %.4f-%.4f Rsun (%.4f-%.4f R_ph)'
            % (s['r_in']['cz_fmlt'], fz['r_base']/RSUN, fz['r_top']/RSUN, fz['r_base']/R,
               fz['r_top']/R))
        for z in fz['zones']:
            log('    convective layer %.4f-%.4f Rsun, max F_MLT/F %.3f' % (z[0]/RSUN,
                                                                         z[1]/RSUN, z[2]))
        nb = AN.hp_between(r, pr['Hp'], x1min, fz['r_base'])
        nb0 = AN.hp_between(r, pr['Hp'], r[0], fz['r_base'])
        log('  H_p buffer below the CZ base: from x1min %.3f H_p, from the column bottom %.3f'
            % (nb, nb0))
        # Schwarzschild boundary: the lowest row above x1min with nabla_rad > nabla_ad that is
        # contiguous with the Fe CZ (the 1 %-F_MLT base sits slightly above it)
        k0 = int(np.searchsorted(r, x1min))
        kb = fz['i_base']
        while kb - 1 > k0 and pr['nrad'][kb-1] > pr['nad'][kb-1]:
            kb -= 1
        rS = r[kb]
        mb = (r >= x1min) & (r < rS)
        if mb.any():
            q = np.nanmax(pr['nrad'][mb]/pr['nad'][mb])
            log('  Schwarzschild CZ base %.4f Rsun; x1min..there: %.3f H_p, max nabla_rad/nabla_ad '
                '%.4f -> %s; nabla_rad/nabla_ad at x1min %.4f'
                % (rS/RSUN, AN.hp_between(r, pr['Hp'], x1min, rS), q,
                   'RADIATIVE' if q < 1 else 'NOT radiative (FAIL)',
                   pr['nrad'][k0]/pr['nad'][k0]))
            if nb < 4.0:
                log('  WARNING: buffer %.2f H_p < 4 (standard: >= 4-5 H_p)' % nb)
        log('  t_Fe = H_p/v_MLT (textbook alpha %g) at r %.3f Rsun: %.4e s -> mlt_ramp_start '
            '%.4e, mlt_ramp_time %.4e' % (st.ALPHA, r[fz['i_turn']]/RSUN, tt, dv['ramp_start'],
                                          dv['ramp_time']))
    i0 = int(np.argmin(abs(r - x1min)))
    log('  at x1min: T %.4e rho %.4e nabla_rad/nabla_ad %.4f F_MLT/F %.2e tau %.4e'
        % (T[i0], rho[i0], pr['nrad'][i0]/pr['nad'][i0], col['fmlt'][i0], pr['tau'][i0]))
    log('  ghost margin: the column starts %.4f Rsun below x1min; he_star_m1 needs it to cover '
        'x1min - 3 dx1(wall) - 2 dx1 -> the grid must have dx1(wall) < ~%.4f Rsun'
        % ((x1min - r[0])/RSUN, (x1min - r[0])/5.5/RSUN))
    log('  implicit_flux_x1min = F_r(x1min) %.12e; L/(4 pi x1min^2) %.12e (ratio %.8f)'
        % (F_rin, L/(4*np.pi*x1min**2), F_rin/(L/(4*np.pi*x1min**2))))
    dV = 4*np.pi*r**2*np.gradient(r)
    me = np.sum((rho*dV)[(r >= x1min) & (r <= R)])
    log('  envelope mass x1min..R_ph %.4e g = %.4f Msun (point-mass GM: M_env/M %.3e)'
        % (me, me/1.98847e33, me*6.67408e-8/st.GM))
    log('---- luminosity consistency')
    m = pr['below'] & (r >= x1min)
    log('  4 pi r^2 (F_r + F_MLT)/L - 1: max |.| %.2e (x1min..R_ph)'
        % np.max(np.abs(pr['Lsum'][m]/L - 1)))
    mr = m & (col['fmlt'] < 1e-6) & (r < col['rj']*0.999) & (r > r[2])
    if mr.any():
        q = pr['Fd'][mr]/col['Fr'][mr]
        log('  radiative rows: F_diff(column T gradient, kappa)/F_r = %.5f..%.5f (consistency '
            'of F_r with the column)' % (q.min(), q.max()))
    log('  tau at x1min %.4e; tau(R_ph) of the column %.4f (2/3 = the shot); R_ph %.6e'
        % (pr['tau'][i0], pr['tau'][int(np.argmin(abs(r - R)))], R))
    log('---- background / floor layer (rules: README)')
    log('  rho_ph %.4e, T_ph %.2f K, g_ph %.4e, c_s,iso %.4e -> g/c_s %.4e /s; v_esc %.4e'
        % (bgv['rho_ph'], Tph, bgv['g_ph'], bgv['cs_ph'], bgv['nu_ph'], bgv['v_esc']))
    log('  column floor %.1e from r %.6e (%.4f R_ph) to the column top: floor layer %.4f-%.4f '
        'R_ph in the domain' % (bgv['rho_floor'], col['r_floor'], col['r_floor']/R,
                                col['r_floor']/R, x1max/R))
    rd = r[(rho < bgv['he_sponge_dmax']) & (r > R)]
    log('  sponge: dmax %.1e %s column floor %.1e -> %s; gas with rho < dmax from r %.4f R_ph;'
        ' r0 %.4f R_ph; rate %.1e /s = %.1f g/c_s (1/rate = %.2e s; dynamical time at R_ph '
        '%.3e s)' % (bgv['he_sponge_dmax'], '>' if bgv['he_sponge_dmax'] > bgv['rho_floor'] else
                     '<=', bgv['rho_floor'],
                     'ACTIVE on the floor gas' if bgv['he_sponge_dmax'] > bgv['rho_floor'] else
                     'INACTIVE (floor above dmax!)', rd.min()/R if rd.size else np.nan,
                     bgv['he_sponge_r0']/R, bgv['he_sponge_rate'],
                     bgv['he_sponge_rate']/bgv['nu_ph'], 1/bgv['he_sponge_rate'],
                     np.sqrt(R**3/st.GM)))
    log('  dfloor %.1e, tfloor %.0f K (T_rad(x1max) %.1f K, Teff %.1f), eos_logt_min %.1f, '
        'pfloor %.1e, vceil %.1e, wb_rmax %.3f R_ph, balance_rmax %.4f R_ph (tau-3 join %.4f)'
        % (bgv['dfloor'], bgv['tfloor_kelvin'], bgv['T_rad_x1max'],
           (L/(4*np.pi*R**2*SIG))**0.25, s['eos']['logt_min'], bgv['pfloor'], bgv['vceil'],
           bgv['wb_rmax']/R, bgv['he_ic_balance_rmax']/R, col['rj']/R))
    if bgv['he_ic_balance_rmax'] >= col['rj']:
        log('  WARNING: he_ic_balance_rmax above the tau-3 join')
    if fz and bgv['he_ic_balance_rmax'] > fz['r_top']:
        log('  NOTE: he_ic_balance_rmax above the Fe-CZ top (E1: inside the MLT++ adiabat)')
    log('---- regime flags')
    teff = (L/(4*np.pi*R**2*SIG))**0.25
    flags = []
    if teff < 5000:
        flags.append('Teff %.0f K < 5000 K: molecules/H2 matter; the EOS dumps have H2 OFF '
                     '(eos_h2 = false)' % teff)
    if bgv['tfloor_kelvin'] < 3500 or T.min() < 4000:
        flags.append('T down to %.0f K / tfloor %.0f K: low-T opacity is Ferguson (grains below'
                     ' ~2000 K); check coverage above' % (T.min(), bgv['tfloor_kelvin']))
    hph = bgv['cs_ph']**2/bgv['g_ph']
    if teff < 5000 or hph/R > 0.03:
        flags.append('extended/cool atmosphere (H_ph/R_ph %.3e): the grey RE/Hopf atmosphere '
                     'and isothermal cap are crude (molecular, non-grey)' % (hph/R))
    if np.nanmax(pr['gam_L'][pr['below']]) > 1:
        flags.append('super-Eddington Gamma_L max %.2f below R_ph (MLT++ carries the excess)'
                     % np.nanmax(pr['gam_L'][pr['below']]))
    if x1max/R > 5:
        flags.append('x1max %.1f R_ph' % (x1max/R))
    for fl_ in flags or ['none (hot-star regime: the floor/sponge rules are calibrated here)']:
        log('  FLAG: ' + fl_)
    if mesa is not None:
        log('---- MESA comparison (%s)' % mesa['path'])
        log('  MESA L_ph %.6e (column L/L_MESA %.6f), M %.6e g (GM/(G M) %.6f), R_ph %s, '
            'Teff %s' % (mesa['L_ph'], L/mesa['L_ph'], mesa['M_tot'],
                         st.GM/(6.67408e-8*mesa['M_tot']),
                         '%.6e (column R/R_MESA %.6f)' % (mesa['R_ph'], R/mesa['R_ph'])
                         if mesa['R_ph'] else 'n/a', mesa['Teff']))
        Lm = MS.at(mesa, 'L', x1min)
        log('  L_MESA(x1min)/L_ph = %.5f (%s, tol %g); m_MESA(x1min) %.6e g'
            % (Lm/mesa['L_ph'], 'FLAT' if abs(Lm/mesa['L_ph'] - 1) <
               s['r_in']['l_flat_tol'] else 'NOT FLAT', s['r_in']['l_flat_tol'],
               MS.at(mesa, 'm', x1min)))
        log('  r/R_ph   T_col/T_MESA  rho_col/rho_MESA')
        for xx in np.linspace(max(r[0], mesa['rc'][0])/R, min(0.995, mesa['rc'][-1]/R), 12):
            rr = xx*R
            i = int(np.argmin(abs(r - rr)))
            log('  %7.4f  %8.4f  %8.4f' % (r[i]/R, T[i]/MS.at(mesa, 'T', r[i], True),
                                          rho[i]/MS.at(mesa, 'rho', r[i], True)))
        for rr_ in s.get('report_radii_rsun', []):
            i = int(np.argmin(abs(r - rr_*RSUN)))
            log('  at %g Rsun: T %.4e (MESA %.4e, ratio %.4f), rho %.4e (MESA %.4e, ratio %.4f)'
                % (rr_, T[i], MS.at(mesa, 'T', r[i], True), T[i]/MS.at(mesa, 'T', r[i], True),
                   rho[i], MS.at(mesa, 'rho', r[i], True),
                   rho[i]/MS.at(mesa, 'rho', r[i], True)))
        mm = (mesa['r'] >= x1min)
        menv = mesa['M_tot'] - MS.at(mesa, 'm', x1min)
        log('  envelope mass above x1min: MESA %.4e g (column above, x1min..R_ph, see above)'
            % menv)
        del mm



def diag_mesa(log, s, st, col, ms, mesa, x1min, czf):
    """mode 'mesa' checks: m_env, R_ph, T_ph vs MESA, CZs coincide, F_rad > F layers"""
    from icstd.micro import G_CODE
    log('---- mode mesa checks')
    R, L = col['rph'], col['L']
    Rm = mesa['R_ph'] or s['star'].get('R_ph_rsun', 0)*RSUN or ms['r'][-1]
    log('  R_ph ours %.6e vs MESA %.6e: ratio %.5f' % (R, Rm, R/Rm))
    Tm = mesa['Teff'] or (L/(4*np.pi*Rm**2*SIG))**0.25
    log('  Teff ours %.2f vs MESA %.2f: ratio %.5f; T_ph(tau 2/3) %.2f'
        % ((L/(4*np.pi*R**2*SIG))**0.25, Tm, (L/(4*np.pi*R**2*SIG))**0.25/Tm, col['tph']))
    mw = col['m_wall']
    menv = col['m_top'] - mw
    mm = MS.at(mesa, 'm', x1min)
    menvM = mesa['M_tot'] - mm
    log('  m(x1min) ours %.8e vs MESA %.8e (%.2e rel); m_env(x1min..tau_j) ours %.5e g = %.5f '
        'Msun vs MESA (to its top) %.5e g = %.5f Msun: ratio %.4f'
        % (mw, mm, mw/mm - 1, menv, menv/1.98847e33, menvM, menvM/1.98847e33, menv/menvM))
    log('  deck: he_gm = G m(x1min) = %.10e, he_gm_column = true' % (G_CODE*mw))
    ii = col['interior']
    ex = col['excess']
    if ex['n']:
        log('  F_rad(our kappa) > F on %d interior rows, %.3f-%.3f Rsun, max excess %.3f: there '
            'F_r = F, F_MLT = 0 (no negative scaffold); radiation will carry %.1f %% more than '
            'L there in the run (our kappa below MESA\'s)' % (ex['n'], ex['rlo']/RSUN,
                                                               ex['rhi']/RSUN, ex['max'],
                                                               100*ex['max']))
    else:
        log('  F_rad(our kappa) <= F everywhere (no clipping)')
    if col.get('ext'):
        log('  OUTWARD MLT++ extension above the MESA profile top: %.4f-%.4f Rsun'
            % (col['ext']['r_e']/RSUN, col['ext']['rj']/RSUN))
    zo = []
    fm = np.clip(1 - ii['frad'], 0, 1)
    i = 0
    while i < fm.size:
        if fm[i] > czf:
            k = i
            while k + 1 < fm.size and fm[k+1] > czf:
                k += 1
            zo.append((ii['r'][i], ii['r'][k], fm[i:k+1].max()))
            i = k + 1
        else:
            i += 1
    zm = [z for z in MC.mesa_zones(ms, czf) if z[1] >= col['r'][0]]   # the column's range
    log('  convective zones (F_MLT/F > %g), ours (our kappa on MESA T(P)):' % czf)
    for z in zo:
        log('    %.4f-%.4f Rsun max %.3f' % (z[0]/RSUN, z[1]/RSUN, z[2]))
    log('  MESA\'s own (MESA kappa):')
    for z in zm:
        log('    %.4f-%.4f Rsun max %.3f' % (z[0]/RSUN, z[1]/RSUN, z[2]))
    if not zm:
        log('    none (radiative everywhere)')
    r = col['r']
    log('  r[Rsun]   T/T_MESA  rho/rho_MESA  m/m_MESA')
    for xx in np.linspace(r[0], min(ii['r'][-1], ms['r'][-1]), 10):
        k = int(np.argmin(abs(ii['r'] - xx)))
        rr = ii['r'][k]
        log('  %8.3f  %8.5f  %8.5f  %10.7f' % (rr/RSUN, col['T'][k]/MS.at(mesa, 'T', rr, True),
                                              col['rho'][k]/MS.at(mesa, 'rho', rr, True),
                                              ii['m'][k]/MS.at(mesa, 'm', rr)))


def classify(col, pr, x1min, M):
    r, rho = col['r'], col['rho']
    m = pr['below'] & (r >= x1min)
    dV = 4*np.pi*r**2*np.gradient(r)
    menv = float(np.sum((rho*dV)[m]))
    inv = bool(np.any(np.diff(np.log(rho[m])) > 1e-6))
    i0 = int(np.argmin(abs(r - x1min)))
    fb = float(col['fmlt'][i0])
    if menv > 0.3*M or fb > 0.9:
        br = 'DENSE ADIABAT (convective to the base / m_env ~ M)'
    elif inv:
        br = 'INFLATED (density inversion below R_ph)'
    else:
        br = 'normal (radiative base, monotone rho)'
    return menv, inv, fb, br


def diag_global(log, s, st, col, pr, c, x1min, no_scan):
    """mode 'global' physicality checks (flag, never accept silently) + branch scan"""
    from icstd.micro import G_CODE
    log('---- mode global physicality checks')
    M = st.GM/G_CODE
    menv, inv, fb, br = classify(col, pr, x1min, M)
    if inv:
        mm_ = pr['below'] & (col['r'] >= x1min)
        rr_ = col['r'][mm_][1:][np.diff(np.log(col['rho'][mm_])) > 1e-6]
        log('  density inversion rows: %.4f-%.4f R_ph (%d rows)' % (rr_.min()/col['rph'],
                                                                 rr_.max()/col['rph'],
                                                                 rr_.size))
    log('  m_env(x1min..R_ph) %.4e g = %.4e M (%s)' % (menv, menv/M,
                                                    'OK << M' if menv < 0.05*M else
                                                    'FLAG: not << M (point-mass gravity wrong)'))
    r, rho = col['r'], col['rho']
    nb = max(10, int(0.05*r.size))
    d2 = np.abs(np.diff(np.log(rho[:nb]), 2))
    d2r = np.abs(np.diff(np.log(rho[nb:nb + 2000]), 2))
    sm = d2.max()/max(np.median(d2r), 1e-30)
    log('  smooth base: max |d2 ln rho| in the bottom %d rows = %.2e (%.1f x the median above) '
        '-> %s' % (nb, d2.max(), sm, 'smooth' if sm < 100 else 'FLAG: kink/jump at the base'))
    log('  branch: %s; F_MLT/F(x1min) %.3f; density inversion %s' % (br, fb, inv))
    if no_scan:
        return
    log('  branch scan (coarse: dlnr 1e-3, rtol 1e-9; same EOS/kappa/L/M/R_ph, bottom %.3f Rsun)'
        % (c['r_bot']/RSUN))
    log('    alpha  MLT++   T(x1min)    rho(x1min)  F_MLT/F(x1min)  m_env/M     branch')
    a0, p0 = st.ALPHA, st.PP_ALPHA
    for al in (1.0, 1.5, 2.0):
        for pp in (st.PP_ALPHA if st.PP_ALPHA > 0 else 1.0, 0.0):
            st.ALPHA, st.PP_ALPHA = al, pp
            cc = dict(c, scheme='agcar_ge', legacy_deep=False, dlnr=1e-3, rtol=1e-9)
            try:
                cs = COL.column_agcar_ge(st, cc)
                ps = AN.profile(st, cs)
                me, iv, f0, b = classify(cs, ps, x1min, M)
                k = int(np.argmin(abs(cs['r'] - x1min)))
                log('    %4.1f  %-5s  %10.4e  %10.4e  %8.3f       %10.3e  %s'
                    % (al, 'on' if pp > 0 else 'off', cs['T'][k], cs['rho'][k], f0, me/M, b))
            except Exception as e:                       # noqa: BLE001
                log('    %4.1f  %-5s  FAILED (%s)' % (al, 'on' if pp > 0 else 'off', str(e)[:80]))
    st.ALPHA, st.PP_ALPHA = a0, p0


def check_tables(log, s, st, col, ed, a):
    """HARD CHECK: every column row inside the EOS dump and the opacity grid (with margin).
    EOS outside and the dump made by the recipe: re-run once with the range widened (0.1 dex
    in T, 0.5 dex in rho beyond the column).  Otherwise FATAL (spec allow_outside_tables:
    true only for reproducing an old file)."""
    lt, lr = np.log10(col['T']), np.log10(col['rho'])
    ex, ey = st.eos.xr, st.eos.yr
    t0, t1, d0, d1 = st.ross_lims
    me, mt = 0.2, 0.05
    eo = (lr < ex[0] + me) | (lr > ex[1] - me) | (lt < ey[0] + mt) | (lt > ey[1] - mt)
    oo = (lr < d0) | (lr > d1) | (lt < t0) | (lt > t1)
    r = col['r']
    msg = []
    if eo.any():
        msg.append('EOS dump %s (log rho %.2f..%.2f, log T %.2f..%.2f; margin %g/%g dex): %d '
                   'rows OUTSIDE, r %.4f..%.4f Rsun, column log rho %.3f..%.3f log T %.3f..%.3f'
                   % (st.eos.fn, ex[0], ex[1], ey[0], ey[1], me, mt, eo.sum(),
                      r[eo].min()/RSUN, r[eo].max()/RSUN, lr.min(), lr.max(), lt.min(),
                      lt.max()))
    if oo.any():
        msg.append('opacity grid (log rho %.2f..%.2f, log T %.2f..%.2f): %d rows OUTSIDE, '
                   'r %.4f..%.4f Rsun' % (d0, d1, t0, t1, oo.sum(), r[oo].min()/RSUN,
                                          r[oo].max()/RSUN))
    if not msg:
        log('table check: every column row inside the EOS dump (with margin) and the '
            'opacity grid')
        return
    for m_ in msg:
        log('TABLE CHECK: ' + m_)
    if s.get('allow_outside_tables'):
        log('TABLE CHECK: allow_outside_tables = true (reproduction of an old file only)')
        return
    if eo.any() and not oo.any() and not s['eos'].get('dump') and \
            not os.environ.get('ICSTD_EOS_RANGE'):
        e = s['eos']
        new = dict(logd_max=float(max(e['logd_max'], np.ceil(lr.max() + 0.5))),
                   logd_min=float(min(e['logd_min'], np.floor(lr.min() - 0.5))),
                   logt_max=float(max(e['logt_max'], np.ceil(10*(lt.max() + 0.1))/10)),
                   logt_min=float(min(e['logt_min'], np.floor(10*(lt.min() - 0.1))/10)))
        log('TABLE CHECK: re-running with the EOS range %s' % new)
        os.environ['ICSTD_EOS_RANGE'] = json.dumps(new)
        sys.stdout.flush()
        os.execv(sys.executable, [sys.executable] + sys.argv)
    raise SystemExit('FATAL table check: ' + ' | '.join(msg))


if __name__ == '__main__':
    main()
