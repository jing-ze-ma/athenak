#!/usr/bin/env python3
"""ACC3D page 10-20: run status, anomaly checks, conservation budgets, hst history and the map figures (READ ONLY on the run).
usage: mk_prod.py <imgdir>      (env: see common.py). Writes <imgdir>/prod.json for build_page.py.
Maps (only for new dumps; FORCE=1 redraws): log rho, log T_gas (EOS table), log P_rad/P_gas (P_rad = E/3 from the m1 dump,
P_gas from the EOS table), v_phi; at t = TS + 0.05, mid-run and the latest hydro dump that has an m1 dump at the same time."""
import os
import re
import sys
import glob
import json
import subprocess
import numpy as np
import common as C
from common import cum, rate, TS, PORB, MDOT, DIRS

IMG = sys.argv[1]
os.makedirs(IMG, exist_ok=True)
plt = C.mpl()
U = C.hstcat('user')
Hh = C.hstcat('hydro')
if U is None:
    raise SystemExit('no ryper.user.hst under %s' % C.RUN)
t = U['time']
st = dict(t=float(t[-1]), dt=float(U['dt'][-1]), run=C.RUN, dirs=C.stream_dirs_label(), chain_cum=C.CHAIN)

# ---------------- run logs ----------------
nretry, nan, fatal, nc, rt, wall_s, segs, lastline = 0, 0, 0, [], [], 0.0, 0, None
for d in DIRS:
    lf = d + '/run.log'
    if not os.path.exists(lf):
        continue
    txt = open(lf, errors='replace').read()
    el = re.findall(r'elapsed=([0-9.e+-]+) cycle=(\d+) time=([0-9.e+-]+) dt=([0-9.e+-]+)', txt)
    prev = 0.0
    for e in el:                        # elapsed restarts at every segment (link / phase): sum the segment maxima
        v = float(e[0])
        if v < prev:
            wall_s += prev; segs += 1
        prev = v
    wall_s += prev; segs += 1 if el else 0
    if el:
        lastline = el[-1]
    r = re.findall(r'^DT_RETRY cycle=\d+ t=([0-9.e+-]+) dt', txt, re.M)
    nretry += len(r); rt += [float(x) for x in r]
    nan += len(re.findall(r'\bnan\b', txt, re.I)); fatal += len(re.findall(r'FATAL|DIVERGED', txt))
    nc += [int(x) for x in re.findall(r'NON-CONVERGED=(\d+)', txt)]
st['cycle'] = int(lastline[1]) if lastline else -1
st['elapsed_h'] = float(lastline[0]) / 3600 if lastline else 0.0
st['wall_h'] = wall_s / 3600; st['segments'] = segs
st['retries'] = nretry; st['retry_t'] = sorted(set(round(x, 3) for x in rt))
st['nan_lines'] = nan; st['fatal_lines'] = fatal; st['nc_max'] = max(nc) if nc else 0
ch = os.path.join(C.RUN, 'chain.txt')
st['chain'] = open(ch).read().strip().splitlines()[-4:] if os.path.exists(ch) else []
st['stopped'] = any('STOPPED' in x for x in st['chain'])
q = []
if C.E('JOBNAME'):
    try:
        q = subprocess.run(['squeue', '-u', os.environ.get('USER', ''), '-h', '-n', C.E('JOBNAME'), '-o', '%i %T %r %M %l'],
                           capture_output=True, text=True, timeout=30).stdout.strip().splitlines()
    except Exception:  # noqa: BLE001
        q = []
st['queue'] = q
st['nrst'] = len(glob.glob(os.path.join(C.RUN, 'stream', 'rst', '*.rst')))
for k in ('Mdom', 'Jdom', 'Etot'):
    st[k + '0'] = float(U[k][0]); st[k + '1'] = float(U[k][-1])
st['Lmeas'] = float(U['Lmeas'][-1]); st['Ltop'] = float(U['Ltop'][-1])
st['hst_nan'] = int(sum(np.isnan(v).sum() for v in U.values()) + (sum(np.isnan(v).sum() for v in Hh.values()) if Hh else 0))
tm, mdin = rate(t, U['Min']); _, mdR = rate(t, U['MR']); _, mdo = rate(t, U['Mout'])
late = tm > max(t[-1] - 0.05, TS)
st['Mdot_in'] = float(np.nanmean(mdin[late])) if late.any() else float('nan')
st['Mdot_R'] = float(np.nanmean(mdR[late])) if late.any() else float('nan')
st['Mdot_out'] = float(np.nanmean(mdo[late])) if late.any() else float('nan')
s = t > TS + 0.02
st['Pic'] = float(np.nanmean(U['Pic'][s][-50:])) if s.any() else float(np.nanmean(U['Pic'][-50:]))
st['Picmax'] = float(np.nanmax(U['Picmax'])) if 'Picmax' in U else float('nan')
jumps = {}
for k in ('Mdom', 'Jdom', 'Etot', 'Erad', 'Lmeas'):
    v = U[k]
    rel = np.abs(np.diff(v)) / (np.abs(v[:-1]) + 1e-30)
    i = int(np.nanargmax(rel)) if len(rel) else 0
    jumps[k] = (float(rel[i]) if len(rel) else 0.0, float(t[i + 1]) if len(rel) else 0.0)
st['jumps'] = jumps
st['dt_min_late'] = float(np.min(U['dt'][-200:]))

# ---------------- conservation budgets (user.hst; all from t = 0) ----------------
d0 = lambda c: cum(U[c]) - cum(U[c])[0]   # noqa: E731
dM = U['Mdom'] - U['Mdom'][0]
rM = dM - (d0('Min') - d0('Mout') - d0('Mmsk') + d0('Mfl'))
dE = U['Etot'] - U['Etot'][0]
rE = dE - d0('Ebnd') - d0('Erel') + d0('Emsk')
dJ = U['Jdom'] - U['Jdom'][0]
rJ = dJ - (d0('Jin') - d0('Jout'))
st['budget'] = dict(dM=float(dM[-1]), Min=float(d0('Min')[-1]), Mout=float(d0('Mout')[-1]), Mmsk=float(d0('Mmsk')[-1]),
                    Mfl=float(d0('Mfl')[-1]), rM=float(rM[-1]), rM_max=float(np.max(np.abs(rM))),
                    dE=float(dE[-1]), Ebnd=float(d0('Ebnd')[-1]), Erel=float(d0('Erel')[-1]), Emsk=float(d0('Emsk')[-1]),
                    rE=float(rE[-1]), Etot0=float(U['Etot'][0]),
                    dJ=float(dJ[-1]), Jin=float(d0('Jin')[-1]), Jout=float(d0('Jout')[-1]), rJ=float(rJ[-1]),
                    MR=float(d0('MR')[-1]), JR=float(d0('JR')[-1]))
ss = t > TS + 1e-6
if ss.sum() > 2:
    i0 = np.argmax(ss)
    st['budget']['rE_stream'] = float(rE[-1] - rE[i0]); st['budget']['Ebnd_stream'] = float(d0('Ebnd')[-1] - d0('Ebnd')[i0])
    st['budget']['rJ_rate'] = float((rJ[-1] - rJ[i0]) / max(t[-1] - t[i0], 1e-30))
    st['budget']['Jin_rate'] = float((d0('Jin')[-1] - d0('Jin')[i0]) / max(t[-1] - t[i0], 1e-30))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
a = ax[0]; a.plot(t, rM, color=C.C1); a.set_title(r'mass: $\Delta M_{dom}$ - (M$_{in}$ - M$_{out}$ - M$_{msk}$ + M$_{fl}$)')
a = ax[1]; a.plot(t, rE, color=C.C1, label='residual'); a.plot(t, d0('Ebnd'), color=C.C2, lw=.8, label=r'E$_{bnd}$ (throughput)')
a.set_title(r'energy: $\Delta E_{tot}$ - E$_{bnd}$ - E$_{rel}$ + E$_{msk}$'); a.legend(); a.set_yscale('symlog', linthresh=1e3)
a = ax[2]; a.plot(t, rJ, color=C.C1, label=r'r$_J$ = $\Delta J_{dom}$ - (J$_{in}$ - J$_{out}$)')
a.plot(t, d0('Jin'), color=C.C2, lw=.8, label=r'J$_{in}$ (nozzle)'); a.plot(t, dJ, color=C.C3, lw=.8, label=r'$\Delta$J$_{dom}$')
a.set_title(r'L$_z$ (inertial, about the accretor axis)'); a.legend()
for a in ax.flat:
    a.axvline(TS, color='#a3adb8', lw=.6, ls=':'); a.grid(alpha=.5); a.set_xlabel('t [code]')
fig.suptitle(f'Conservation budgets from user.hst ({C.SITE}, {C.RUN_LABEL}); code units', color='#e8ecf1')
fig.tight_layout(); fig.savefig(IMG + '/fig_budget.png', dpi=80); plt.close(fig)

# ---------------- history ----------------
fig, ax = plt.subplots(3, 3, figsize=(15, 11.5))
a = ax[0, 0]; a.plot(t, U['dt'], color=C.C1, lw=1); a.set_yscale('log'); a.set_title('dt [code]')
for x in st['retry_t']:
    a.axvline(x, color=C.C5, lw=.6, alpha=.6)
a = ax[0, 1]; a.plot(t[1:], U['Pic'][1:], color=C.C1, lw=1, label='mean'); a.plot(t[1:], U['Picmax'][1:], color=C.C2, lw=.8, label='max')
a.set_yscale('log'); a.set_title('Picard iterations per step (hst interval)'); a.legend()
a = ax[0, 2]; a.plot(t, U['Lmeas'] / 1e4, color=C.C1, label=r'L$_{\rm meas}$ (r 9.57)')
a.plot(t, U['Ltop'] / 1e4, color=C.C2, label=r'L$_{\rm top}$ (r$_{\rm out}$, incl. donor light)')
a.set_title(r'luminosity [10$^4$ L$_\odot$]'); a.legend()
a = ax[1, 0]
a.plot(t, dM, color=C.C1, label=r'M$_{\rm dom}$ - M$_{\rm dom}$(0)')
a.plot(t, d0('Min'), color=C.C2, lw=1, label=r'$\int\dot M_{\rm in}$ (nozzle)')
a.plot(t, d0('Mout'), color=C.C3, lw=1, label=r'$\int\dot M_{\rm out}$ (boundaries)')
a.plot(t, d0('Mfl'), color=C.C4, lw=1, label='floor mass added')
a.set_title('domain mass change and cumulative boundary terms [code]'); a.legend(fontsize=7)
a = ax[1, 1]
a.plot(tm, mdin, color=C.C2, lw=1, label=r'$\dot M_{\rm in}$ nozzle')
a.plot(tm, mdR, color=C.C1, lw=1.2, label=r'$\dot M_{\rm R}$ (hst MR)')
a.plot(tm, mdo, color=C.C3, lw=1, label=r'$\dot M_{\rm out}$ boundaries')
a.axhline(MDOT, color='#a3adb8', lw=.6, ls='--')
a.set_title(r'mass rates [code; %.1f = 1e-4 M$_\odot$/yr]' % MDOT); a.legend(fontsize=7); a.set_xlim(TS - 0.01, None)
a = ax[1, 2]
a.plot(t, U['Jdom'] / U['Jdom'][0] - 1, color=C.C1, label=r'J$_{\rm dom}$')
a.plot(t, U['Jenv'] / U['Jenv'][0] - 1, color=C.C2, lw=1, label=r'J$_{\rm env}$')
a.set_title(r'angular momentum J$_z$: relative change since t = 0'); a.legend()
a = ax[2, 0]
a.plot(t, U['Etot'] / U['Etot'][0] - 1, color=C.C1, label=r'E$_{\rm tot}$ (user)')
if Hh:
    a.plot(Hh['time'], Hh['tot-E'] / Hh['tot-E'][0] - 1, color=C.C2, lw=1, label='tot-E (hydro.hst)')
a.set_title('total energy: relative change since t = 0'); a.legend()
a = ax[2, 1]
a.plot(t, U['KE'], color=C.C1, label='KE (user, r < r$_{meas}$)')
if Hh:
    a.plot(Hh['time'], Hh['1-KE'] + Hh['2-KE'] + Hh['3-KE'], color=C.C2, lw=1, ls='--', label='KE (hydro.hst, domain)')
a.set_yscale('log'); a.set_title('kinetic energy (rotating frame) [code]'); a.legend()
a = ax[2, 2]
a.plot(t, U['Erad'], color=C.C1, label=r'E$_{\rm rad}$'); a.set_title('radiation energy in the domain [code]'); a.legend()
for a in ax.flat:
    a.axvline(TS, color='#a3adb8', lw=.6, ls=':'); a.grid(alpha=.5)
    a.set_xlabel(f't [code]  (P$_{{orb}}$ = {PORB})')
fig.suptitle(f'{C.RUN_LABEL} ({C.SITE}, {C.LAYOUT}): hst to t = {t[-1]:.4f} = {t[-1] / PORB:.3f} P$_{{orb}}$; '
             f'dotted = stream switch-on (t {TS}); red ticks in dt = dt retries', color='#e8ecf1', fontsize=12)
fig.tight_layout(); fig.savefig(IMG + '/fig_prod_hist.png', dpi=80); plt.close(fig)
print('hist + budget done', flush=True)

# ---------------- maps ----------------
SPEC = [('rho', r'log$_{10}\rho$ [code]', 'magma', -8, 2.5),
        ('lT', r'log$_{10}T_{\rm gas}$ [K] (EOS table)', 'inferno', 3.6, 6.0),
        ('lbeta', r'log$_{10}P_{\rm rad}/P_{\rm gas}$ (P$_{\rm rad}$ = E/3)', 'RdBu_r', -2, 5),
        ('vphi', r'v$_\phi$ [km/s] (rotating frame)', 'coolwarm', -300, 300)]


def fields(h, m1f):
    th = C.thermo(h, m1f)
    return {'rho': np.log10(np.maximum(h['dens'], 1e-30)), 'lT': th['lT'], 'lbeta': th.get('lbeta'), 'vphi': h['velz']}, th


def prodslice(hf, m1f, fn, center=None):
    h = C.load(hf, ['dens', 'eint', 'velz'])
    q, th = fields(h, m1f)
    tt = h['Time']; nj = h['dens'].shape[1]; js = [nj // 2 - 1, nj // 2]
    rf, pf, tf = h['x1f'], h['x3f'], h['x2f']
    RMAX = rf[-1]
    ir = len(rf) - 1
    PP, RR = np.meshgrid(pf, rf, indexing='ij'); X, Y = RR * np.cos(PP), RR * np.sin(PP)
    TT, RR2 = np.meshgrid(tf, rf, indexing='ij')
    if center is None:   # impact zoom: densest stream gas at r 9.8-13 on the equator
        rv = h['x1v']; sel = (rv > 9.8) & (rv < 13.0)
        z = h['dens'][:, js, :].mean(axis=1)[:, sel]
        kk, _ = np.unravel_index(np.argmax(z), z.shape)
        ph = h['x3v'][kk]; r = 10.6
        center = (float(r * np.cos(ph)), float(r * np.sin(ph)), float(ph))
    cx, cy, cph = center
    kc = int(np.argmin(np.abs(((h['x3v'] - cph + np.pi) % (2 * np.pi)) - np.pi)))
    k0 = 0; k1 = int(np.argmin(np.abs(h['x3v'] - np.pi)))
    fig, ax = plt.subplots(4, 4, figsize=(18, 19.5), gridspec_kw=dict(height_ratios=[1.2, 1.0, 0.62, 1.0]))
    for i, (key, lab, cm, lo, hi) in enumerate(SPEC):
        if q[key] is None:
            for row in range(4):
                ax[row, i].set_axis_off()
            continue
        ze = q[key][:, js, :ir].mean(axis=1)
        for row, (xl, yl, ttl) in enumerate([((-RMAX, RMAX), (-RMAX, RMAX), 'equator'),
                                             ((cx - 2.5, cx + 2.5), (cy - 2.5, cy + 2.5), 'impact zoom')]):
            a = ax[row, i]
            pc = a.pcolormesh(X, Y, ze, cmap=cm, **C.nrm(key, lo, hi), shading='flat', rasterized=True)
            if row == 0:
                a.add_patch(plt.Circle((0, 0), rf[0], color='#3a424c'))
                a.add_patch(plt.Rectangle((cx - 2.5, cy - 2.5), 5, 5, fill=False, ec=C.C2, lw=1))
                a.annotate(r'donor / L1 $\rightarrow$', xy=(RMAX * .97, -RMAX * .93), ha='right', fontsize=8, color=C.C2)
            a.set_xlim(*xl); a.set_ylim(*yl); a.set_aspect('equal')
            a.set_title(f'{ttl}: ' + lab, fontsize=9); a.set_xlabel('x'); a.set_ylabel('y')
            fig.colorbar(pc, ax=a, shrink=.8)
        for row, (kks, xl, yl, ttl) in enumerate([(((k0, 1), (k1, -1)), (-RMAX, RMAX), (-9, 9), r'merid. $\phi$ 0 / 180'),
                                                  (((kc, 1),), (RIN_, RMAX), (-9, 9), 'merid. through impact')]):
            a = ax[2 + row, i]
            for kk, sg in kks:
                pc = a.pcolormesh(sg * RR2 * np.sin(TT), RR2 * np.cos(TT), q[key][kk, :, :ir], cmap=cm, **C.nrm(key, lo, hi),
                                  shading='flat', rasterized=True)
            a.set_xlim(*xl); a.set_ylim(*yl); a.set_aspect('equal')
            a.set_title(f'{ttl}: ' + lab + (f' ($\\phi$ {np.degrees(h["x3v"][kc]):.0f})' if row else ''), fontsize=9)
            a.set_xlabel('x' if row == 0 else r'R$_{cyl}$'); a.set_ylabel('z')
            fig.colorbar(pc, ax=a, shrink=.8)
    fig.suptitle(f"{C.RUN_LABEL} ({C.SITE}): t = {tt:.4f} ({tt / PORB:.3f} P_orb); hydro {os.path.basename(hf)}"
                 f" + {os.path.basename(m1f) if m1f else 'no m1 dump'}  (Cartesian, R_sun, donor +x)", color='#e8ecf1', fontsize=12)
    fig.tight_layout(); fig.savefig(fn, dpi=72, pil_kwargs={'quality': 80}); plt.close(fig)
    rng = {k: [float(np.nanpercentile(q[k], 0.5)), float(np.nanpercentile(q[k], 99.5))] for k in q if q[k] is not None}
    print('map', fn, 'ranges p0.5-99.5', rng, flush=True)
    return center, float(tt)


RIN_ = C.RIN
HL = C.dumps('hydro_w'); ML = C.dumps('m1')
pairs = [(th, f, C.m1_for(th, ML)) for th, f, _ in HL]
pairs = [p for p in pairs if p[2]]
cfile = IMG + '/prod_center.json'
cen = tuple(json.load(open(cfile))) if os.path.exists(cfile) and not C.FORCE else None
want = []
if pairs:
    last = pairs[-1]
    near = lambda tv: min(pairs, key=lambda p: abs(p[0] - tv))   # noqa: E731
    want = [near(TS + 0.05), near(0.5 * (TS + 0.05 + last[0])), last]
seen, maps = set(), []
for th, hf, mf in want:
    n = C.dnum(hf)
    if n in seen:
        continue
    seen.add(n)
    fn = f'{IMG}/fig_prod_map_{n:05d}.jpg'
    if C.FORCE or not os.path.exists(fn):
        cen, _ = prodslice(hf, mf, fn, cen)
        json.dump(list(cen), open(cfile, 'w'))
    maps.append(dict(file=os.path.basename(fn), nh=n, nm=C.dnum(mf), t=th))
for f in glob.glob(IMG + '/fig_prod_map_*.jpg'):
    if os.path.basename(f) not in [m['file'] for m in maps]:
        os.remove(f)
st['maps'] = maps
st['ndumps'] = len(HL); st['nm1'] = len(ML)
st['last_dump_t'] = HL[-1][0] if HL else None
# theta-wall screen: mass fraction at r 9.5-16.2 in the 3 theta rows next to each wall, latest hydro dump (flag at 1 %)
wfile = IMG + '/prod_wall.json'
wall = json.load(open(wfile)) if os.path.exists(wfile) and not C.FORCE else {}
if HL and str(C.dnum(HL[-1][1])) not in wall:
    h = C.load(HL[-1][1], ['dens'])
    r = h['x1v']; nj = len(h['x2v'])
    sel = (r > 9.5) & (r < 16.2)
    m = (h['dens'] * C.dvol(h))[:, :, sel].sum(axis=(0, 2))
    wall[str(C.dnum(HL[-1][1]))] = dict(t=float(h['Time']), frac=float((m[:3].sum() + m[-3:].sum()) / m.sum()), uniform=6.0 / nj)
    json.dump(wall, open(wfile, 'w'))
if wall:
    st['wall'] = wall[max(wall, key=int)]
json.dump(st, open(IMG + '/prod.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in st.items() if k not in ('chain',)}, indent=1))
print('prod done')
