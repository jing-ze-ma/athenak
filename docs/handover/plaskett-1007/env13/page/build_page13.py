#!/usr/bin/env python3
"""Build the "Plaskett env13 stream run" page (index.html) from figs/numbers_*.json, the
figures made by mkfigs13.py, the run logs and the two inputs (env12, env13).

usage: build_page13.py [--arm q13_s1] [--run <test1007>/run/q13_s1] [--standin]
                       [--state auto|no-data|queued|running|finished|stopped|failed]
                       [--figs DIR] [--out DIR] [--embed]

Same layout and style as test1007/build_page.py (the q12 page, "Plaskett Stream Impact").
Figures are referenced as figs/<name>.png next to index.html (publish them as supporting
files), or inlined as base64 with --embed (one self-contained file).
"""
import argparse
import base64
import datetime
import glob
import html
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
T = '/work/nvme/bivj/jma20/plaskett_1007/test1007'
RUN12 = '/work/nvme/bivj/jma20/plaskett_1007/run'
IN12 = '/work/nvme/bivj/jma20/plaskett_1007/in/plaskett_env12.athinput'
IN13 = f'{T}/in/plaskett_env13.athinput'
TLIM = 0.2291
POLD = [('q12_s1', 'env12', 'PLM + HLLC, nghost 2', '3333330'),
        ('q12f_s1', 'env12f', 'PLM + HLLC + FOFC, nghost 3', '3333331')]

ap = argparse.ArgumentParser()
ap.add_argument('--arm', default='q13_s1')
ap.add_argument('--run', default=f'{T}/run/q13_s1')
ap.add_argument('--standin', action='store_true',
                help='the --run dir is a stand-in for q13_s1 (test of the pipeline)')
ap.add_argument('--state', default='auto')
ap.add_argument('--figs', default=os.path.join(HERE, 'figs'))
ap.add_argument('--out', default=HERE)
ap.add_argument('--embed', action='store_true')
a = ap.parse_args()
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M CDT')
ARM = a.arm


# ------------------------------------------------------------------ helpers
def last_cycle(run):
    try:
        lines = [ln for ln in open(f'{run}/run.log', errors='replace') if 'cycle=' in ln]
        for ln in reversed(lines):
            m = re.search(r'elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)', ln)
            if m:
                el, cyc, t, dt = float(m[1]), int(m[2]), float(m[3]), float(m[4])
                break
        else:
            return None
        m = re.search(r'^time=(\S+) cycle=(\d+)', lines[-1])   # final line at tlim
        if m and float(m[1]) >= t:
            t, cyc = float(m[1]), int(m[2])
        return el, cyc, t, dt
    except Exception:
        return None


def job_log(arm):
    """Newest slurm log of the job that actually ran this arm (skips 'refusing' copies)."""
    best = None
    for fn in glob.glob(f'{T}/logs/{arm}.*.out'):
        txt = open(fn, errors='replace').read()
        if 'refusing' in txt or ' start ' not in txt:
            continue
        if best is None or os.path.getmtime(fn) > os.path.getmtime(best):
            best = fn
    return best


def run_state(run, arm):
    """no-data | running | finished | stopped | failed (from run.log and the job log)."""
    log = f'{run}/run.log'
    if not os.path.isfile(log):
        return 'no-data', 'no run.log yet (job not started)'
    txt = open(log, errors='replace').read()
    if re.search(r'FATAL ERROR|### FATAL', txt):
        return 'failed', 'FATAL ERROR in run.log'
    if 'Terminating on time limit' in txt:
        return 'finished', 'reached tlim'
    if re.search(r'Terminating on (wall|cycle)', txt):
        m = re.search(r'Terminating on [^\n]*', txt)
        return 'stopped', m[0]
    jl = job_log(arm)
    if jl:
        m = re.search(r'^rc (\d+)', open(jl, errors='replace').read(), re.M)
        if m:
            return ('stopped' if m[1] == '0' else 'failed'), f'job ended rc {m[1]} before tlim'
    if 'cpu time used' in txt:
        return 'stopped', 'ended before tlim'
    return 'running', 'run.log growing'


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and x != x):
        return '—'
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return f'{x:d}'
    if isinstance(x, str):
        return html.escape(x)
    if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e4):
        return f'{x:.{nd-1}e}'.replace('e-0', 'e-').replace('e+0', 'e')
    return f'{x:.{nd}g}' if abs(x) < 1 else f'{x:.{nd}f}'.rstrip('0').rstrip('.')


def keys(fn):
    par, blk = {}, None
    try:
        for line in open(fn):
            c = line.split('#', 1)
            ln = c[0].strip()
            m = re.match(r'<(\S+)>', ln)
            if m:
                blk = m[1]
                continue
            if '=' in ln and blk:
                k, v = ln.split('=', 1)
                par[k.strip()] = v.strip()
    except OSError:
        pass
    return par


def load_nums(arm):
    fn = os.path.join(a.figs, f'numbers_{arm}.json')
    try:
        return json.load(open(fn))
    except Exception:
        return {'arm': arm, 'no_data': True, 'reason': 'no numbers file'}


def img(name, alt, missing):
    p = os.path.join(a.figs, name)
    if not os.path.isfile(p):
        return f'<div class="nodata">{html.escape(missing)}</div>'
    if a.embed:
        src = 'data:image/png;base64,' + base64.b64encode(open(p, 'rb').read()).decode()
    else:
        src = f'figs/{name}'
    return f'<img src="{src}" alt="{html.escape(alt)}" loading="lazy">'


# ------------------------------------------------------------------ data
st13, why13 = (a.state, 'set by hand') if a.state != 'auto' else run_state(a.run, ARM)
n13 = load_nums(ARM)
if st13 == 'no-data':
    n13 = {'arm': ARM, 'no_data': True}
nums = {ARM: n13}
for arm, *_ in POLD:
    nums[arm] = load_nums(arm)
ARMS = [ARM] + [p[0] for p in POLD]
k12, k13 = keys(IN12), keys(IN13)
jl13 = job_log(ARM) if not a.standin else None
job13 = re.search(r'\.(\d+)\.out$', jl13)[1] if jl13 else ('stand-in' if a.standin else '—')
q13_label = f'{ARM} <small>(stand-in: {html.escape(os.path.basename(a.run))})</small>' \
    if a.standin else ARM


def card(arm, title, sub, run, state, note, tlim=TLIM):
    lc = last_cycle(run)
    tl = tlim
    if a.standin and arm == ARM and lc:
        tl = float(nums[arm].get('t_final_code') or lc[2]) if state == 'finished' else tlim
    pct = 100*lc[2]/tl if lc else 0
    return f'''
      <div class="arm">
        <div class="arm-head"><span class="arm-name">{title}</span><span class="pill {state}">{html.escape(state)}</span></div>
        <div class="arm-sub">{sub}</div>
        <div class="bar" role="img" aria-label="{pct:.0f} percent of tlim"><span style="width:{min(pct, 100):.1f}%"></span></div>
        <dl class="kv">
          <div><dt>t / t<sub>lim</sub></dt><dd>{fmt(lc[2]) if lc else '—'} / {fmt(tl)} <small>({pct:.0f} %)</small></dd></div>
          <div><dt>cycle</dt><dd>{lc[1] if lc else '—'}</dd></div>
          <div><dt>dt now</dt><dd>{fmt(lc[3]) if lc else '—'}</dd></div>
          <div><dt>wall</dt><dd>{f"{lc[0]/60:.0f} min" if lc else '—'}</dd></div>
        </dl>
        <div class="arm-note">{html.escape(note)}</div>
      </div>'''


status_rows = card(ARM, q13_label,
                   f'env13 · PLM + HLLC, nghost 2 · nx1 {k13.get("nx1", "640")} · job {job13}',
                   a.run, st13, why13)
for arm, inp, scheme, job in POLD:
    stt, why = run_state(f'{RUN12}/{arm}', arm)
    status_rows += card(arm, arm, f'{inp} · {scheme} · nx1 {k12.get("nx1", "500")} · job {job}',
                        f'{RUN12}/{arm}', stt, why)

# ------------------------------------------------------------------ what changed
CHG = [
    ('Envelope top', 'sphere through R<sub>acc</sub> at φ 90° (r<sub>top</sub> 9.18): the tidal caps toward L1/L2 '
     'were missing at t = 0', '<code>env_top_mode = equipotential</code>: ends on ψ = −15 c<sub>ph</sub><sup>2</sup> '
     'at every φ (r<sub>top</sub> 9.16 at φ 90°, 9.69 at φ 0°)',
     'the photosphere follows the Roche equipotential (9.48 at φ 0°, 9.30 at 180°); the q12 IC blew out and '
     'breathed through R<sub>acc</sub>'),
    ('Well-balancing cut', f'<code>hydro/wb_rmax = {k12.get("wb_rmax", "8.52")}</code> (sphere)',
     f'<code>env_wb_depth = {k13.get("env_wb_depth", "0.41")}</code> R<sub>☉</sub> below the photosphere, '
     'on an equipotential', '1.5 penetration depths of the Plaskett impact (d<sub>pen</sub> 0.27); a sphere cut '
     'sits at different depths on the bulged star'),
    ('Measuring radius', f'R<sub>acc</sub> = {fmt(float(k12.get("r_acc", "9.0013")), 5)}',
     f'<code>r_meas = {k13.get("r_meas", "9.615675")}</code> (face above max r<sub>ph</sub> + 10 H<sub>p</sub>)',
     'R<sub>acc</sub> lies inside the star on the binary axis; MR/JR and Menv/Jenv now refer to r<sub>meas</sub>'),
    ('Stream width', 'pgen default c<sub>s</sub>/Ω = 1.53 R<sub>☉</sub>',
     f'<code>stream_width = {k13.get("stream_width", "0.735")}</code> R<sub>☉</sub> at r<sub>out</sub>',
     'Ryu et al. 2025 isothermal L1 solution, carried to r<sub>out</sub> with its transverse pressure-supported width'),
    ('Photospheric density', f'<code>env_rho_ph = {k12.get("env_rho_ph", "0.3")}</code>',
     f'<code>env_rho_ph = {k13.get("env_rho_ph", "0.055")}</code>',
     'ρ<sub>ph</sub> 2.82e−9 g/cc over the stream peak at r<sub>out</sub> 5.11e−8 g/cc (Mdot 1e−4 M<sub>☉</sub>/yr); '
     '0.3 used A = (c/Ω)<sup>2</sup> and 650 km/s'),
    ('Temperatures', f'c<sub>ph</sub> {k12.get("env_cs_ph", "20")} km/s, T<sub>don</sub> {k12.get("t_don", "33000")} K',
     f'c<sub>ph</sub> {k13.get("env_cs_ph", "19.4132")} km/s, T<sub>don</sub> {k13.get("t_don", "26200")} K',
     'Wade et al. 2026 track at the onset of Case A (gainer 28.3 kK, donor 26.2 kK)'),
    ('Ambient / floor', f'env_amb_rho {k12.get("env_amb_rho", "1e-6")}, dfloor {k12.get("dfloor", "1e-7")}',
     f'env_amb_rho {k13.get("env_amb_rho", "1e-7")}, dfloor {k13.get("dfloor", "1e-8")}',
     'scaled with the lower photospheric density'),
    ('Radial grid', f'{k12.get("nx1", "500")} cells, polynomial + plateau stretch',
     f'{k13.get("nx1", "640")} cells, pure plateau: dr 2.48e−3 ≤ H<sub>p</sub>/4 over 8.76–9.64',
     'resolves the whole bulged photosphere; neighbour ratio ≤ 1.095'),
    ('Spin', f'{k12.get("spin", "1.0")}', f'{k13.get("spin", "1.0")} (kept)',
     'user decision 10-07: the 7.2 of the spin-arm plan was an RY Per carry-over, above Plaskett\'s critical 4.72'),
    ('Code', '<code>accretor-1006</code> ea04cd19', '<code>accretor-1007</code> 8cecb89a',
     'new keys are bitwise-off by default (T1: old/new binaries identical on the q12 input)'),
]
chg_rows = ''.join(f'<tr><th scope="row">{k}</th><td class="t">{o}</td><td class="t">{n}</td>'
                   f'<td class="t w">{w}</td></tr>' for k, o, n, w in CHG)

# ------------------------------------------------------------------ numbers table
ROWS = [('t_final_orbits', 'Time reached', 'orbits'),
        ('r_meas', 'Measuring radius', 'R<sub>☉</sub>'),
        ('M_stream_in_cum', 'Stream mass in (cumulative)', 'ρ<sub>s</sub> R<sub>☉</sub><sup>3</sup>'),
        ('M_through_Racc_cum', 'Net mass in through the measuring face', 'ρ<sub>s</sub> R<sub>☉</sub><sup>3</sup>'),
        ('accreted_fraction', 'Accreted fraction (net in / stream in)', ''),
        ('M_through_rmeas_min', 'Lowest net cumulative mass through the face', 'ρ<sub>s</sub> R<sub>☉</sub><sup>3</sup>'),
        ('M_through_rmeas_max', 'Highest net cumulative mass through the face', 'ρ<sub>s</sub> R<sub>☉</sub><sup>3</sup>'),
        ('edge_dr_median_max', 'Envelope edge displacement, median over φ (max in time)', 'R<sub>☉</sub>'),
        ('edge_dr_max_max', 'Envelope edge displacement, max over φ and time', 'R<sub>☉</sub>'),
        ('edge_minus_rph_phi180_initial', 'Edge above the photosphere at φ 180°, first dump', 'R<sub>☉</sub>'),
        ('edge_minus_rph_phi180_final', 'Edge above the photosphere at φ 180°, last dump', 'R<sub>☉</sub>'),
        ('j_through_Racc_cum_over_jK', 'j through the face, cumulative', 'j<sub>K</sub>(r<sub>meas</sub>)'),
        ('j_through_Racc_last0p1orb_over_jK', 'j through the face, last 0.1 orbit', 'j<sub>K</sub>(r<sub>meas</sub>)'),
        ('j_stream_over_jK', 'j brought in by the stream', 'j<sub>K</sub>(r<sub>meas</sub>)'),
        ('j_ballistic_over_jK', 'j of the ballistic L1 orbit at r<sub>meas</sub>', 'j<sub>K</sub>(r<sub>meas</sub>)'),
        ('j_star_rotation_Racc_over_jK', 'j of rigid rotation at r<sub>meas</sub>', 'j<sub>K</sub>(r<sub>meas</sub>)'),
        ('j_Kep_Racc_code', 'j<sub>K</sub>(r<sub>meas</sub>)', 'R<sub>☉</sub> km/s'),
        ('fatal_count', 'FATAL messages', ''),
        ('dt_min_code', 'Smallest dt', 'code')]


def cell(arm, k):
    n = nums[arm]
    if n.get('no_data'):
        return '<td class="na">no data</td>'
    return f'<td>{fmt(n.get(k))}</td>'


table = ''.join(f'<tr><th scope="row">{lab}</th>' + ''.join(cell(arm, k) for arm in ARMS)
                + f'<td class="u">{u}</td></tr>' for k, lab, u in ROWS)
thead = ''.join(f'<th scope="col">{ARM + (" *" if a.standin else "") if arm == ARM else arm}</th>'
                for arm in ARMS)

# ------------------------------------------------------------------ headline (auto text)
n12 = nums['q12_s1']


def headline():
    if n13.get('no_data'):
        return ('<p><b>No q13 data yet.</b> The q12 columns below are the old-setup reference; '
                'the q13 column fills in when the run has written its first history lines.</p>')
    out = []
    f13, f12 = n13.get('accreted_fraction'), n12.get('accreted_fraction')
    if f13 is not None and f12 is not None:
        out.append(f'Net mass in through the measuring face over stream inflow: <b>{fmt(f13)}</b> '
                   f'for {ARM} (r<sub>meas</sub> {fmt(n13.get("r_meas"), 4)}) against '
                   f'<b>{fmt(f12)}</b> for q12_s1 (R<sub>acc</sub> {fmt(n12.get("r_meas"), 4)}).')
    elif f13 is None:
        out.append(f'{ARM} has no stream inflow yet (or runs without a stream), so no accreted fraction.')
    lo13, lo12 = n13.get('M_through_rmeas_min'), n12.get('M_through_rmeas_min')
    if lo13 is not None and lo12 is not None:
        out.append(f'Lowest net cumulative mass through the face (outflow = breathing): {fmt(lo13)} '
                   f'against {fmt(lo12)} ρ<sub>s</sub> R<sub>☉</sub><sup>3</sup> in q12_s1.')
    e13, e12 = n13.get('edge_dr_median_max'), n12.get('edge_dr_median_max')
    if e13 is not None and e12 is not None:
        out.append(f'Envelope edge moves by {fmt(e13)} R<sub>☉</sub> (median over φ, away from the stream) '
                   f'against {fmt(e12)} R<sub>☉</sub> in q12_s1.')
    j13, j12 = n13.get('j_through_Racc_last0p1orb_over_jK'), n12.get('j_through_Racc_last0p1orb_over_jK')
    if j13 is not None and j12 is not None:
        out.append(f'j of the gas through the face over the last 0.1 orbit: {fmt(j13)} j<sub>K</sub> '
                   f'against {fmt(j12)} j<sub>K</sub> (ballistic {fmt(n13.get("j_ballistic_over_jK"))} / '
                   f'{fmt(n12.get("j_ballistic_over_jK"))}).')
    t13 = n13.get('t_final_orbits')
    if t13 is not None and t13 < 0.49:
        out.append(f'<i>{ARM} has reached {fmt(t13)} orbit only; q12 values are at 0.5 orbit, so the '
                   f'cumulative numbers are not yet comparable.</i>')
    return '<ul class="caveats">' + ''.join(f'<li>{s}</li>' for s in out) + '</ul>'


# ------------------------------------------------------------------ banner
if a.standin:
    banner = (f'<b>STAND-IN TEST.</b> The {ARM} column, cards and figures are made from '
              f'<code>{html.escape(a.run)}</code>, not from the production run. '
              f'Built {now} to test the pipeline only.')
    bclass = 'banner warn'
elif st13 == 'finished':
    banner = f'{ARM} reached the end of the half orbit. q12_s1 / q12f_s1 are the old-setup reference runs.'
    bclass = 'banner'
elif st13 == 'no-data':
    banner = (f'Snapshot of {now}. {ARM} has not started yet (no run.log); only the q12 reference '
              f'is shown.')
    bclass = 'banner'
elif st13 == 'running':
    banner = (f'Snapshot of {now}. {ARM} is still running; figures and numbers are partial and '
              f'will be replaced when it ends.')
    bclass = 'banner'
else:
    banner = (f'{ARM} ended without reaching t<sub>lim</sub> ({html.escape(st13)}: '
              f'{html.escape(why13)}). Figures show the data up to where it stopped.')
    bclass = 'banner warn'


def plate(name, alt, cap, missing):
    return (f'<figure>{img(name, alt, missing)}<figcaption>{cap}</figcaption></figure>')


nd13 = f'{ARM}: no data yet'
map13 = plate(f'fig_{ARM}_maps.png', f'theta-mean density maps of {ARM}',
              f'<b>{ARM}.</b> θ-mean density in units of the stream peak. Top: full domain. Bottom: impact zoom '
              '(cyan box). White dashed: ballistic L1 orbit; grey disk: inner boundary; solid blue: '
              'R<sub>acc</sub>; dashed light blue: r<sub>meas</sub>.', nd13)
bud13 = plate(f'fig_{ARM}_budget.png', f'mass and AM budget of {ARM}',
              f'<b>{ARM}.</b> Left: cumulative stream inflow, net mass in through r<sub>meas</sub>, '
              'outflow through r<sub>out</sub>. Right: inertial specific AM of the stream and of the gas '
              'through r<sub>meas</sub>, in j<sub>Kep</sub>(r<sub>meas</sub>).', nd13)

page = f'''<title>Plaskett env13 stream run</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:wght@500;600&family=Public+Sans:wght@400;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
/* layout: single reading column, wide figure plates that break out to 1100px; status strip on top */
:root {{
  --bg: #f4f5f7; --surface: #ffffff; --fg: #1b1f27; --muted: #5b6372; --line: #d9dde4;
  --accent: #c2410c; --accent-soft: #fde7da; --ok: #15803d; --ok-soft: #dcfce7; --run: #1d4ed8; --run-soft: #dbe5fe;
  --bad: #b91c1c; --bad-soft: #fee2e2;
  --f-display: 'Spectral', Georgia, serif; --f-body: 'Public Sans', system-ui, sans-serif;
  --f-mono: 'JetBrains Mono', ui-monospace, Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #12151b; --surface: #1b1f27; --fg: #e6e8ec; --muted: #9aa3b2; --line: #2e3440;
  --accent: #fb923c; --accent-soft: #3b2216; --ok: #4ade80; --ok-soft: #12301d; --run: #93b4ff; --run-soft: #1c2742;
  --bad: #f87171; --bad-soft: #3a1717; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #12151b; --surface: #1b1f27; --fg: #e6e8ec; --muted: #9aa3b2; --line: #2e3440;
  --accent: #fb923c; --accent-soft: #3b2216; --ok: #4ade80; --ok-soft: #12301d; --run: #93b4ff; --run-soft: #1c2742;
  --bad: #f87171; --bad-soft: #3a1717; color-scheme: dark; }}
body {{ background: var(--bg); color: var(--fg); font: 16px/1.6 var(--f-body); margin: 0; }}
.wrap {{ max-width: 1100px; margin: 0 auto; padding-inline: 16px; padding-block: 40px 64px; display: grid; gap: 40px; }}
.prose {{ max-width: 68ch; }}
h1, h2 {{ font-family: var(--f-display); font-weight: 600; text-wrap: balance; line-height: 1.2; margin: 0; }}
h1 {{ font-size: clamp(30px, 5vw, 44px); }}
h2 {{ font-size: 26px; margin-bottom: 12px; }}
.eyebrow {{ font: 500 12px/1 var(--f-mono); letter-spacing: .08em; text-transform: uppercase; color: var(--accent); }}
header {{ display: grid; gap: 12px; }}
.lede {{ color: var(--muted); margin: 0; max-width: 68ch; }}
.banner {{ background: var(--accent-soft); border-left: 3px solid var(--accent); padding: 10px 14px; font-size: 14px; }}
.banner.warn {{ background: var(--bad-soft); border-left-color: var(--bad); }}
.status {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; }}
.arm {{ background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 16px; display: grid; gap: 10px; min-width: 0; }}
.arm-head {{ display: flex; justify-content: space-between; align-items: center; gap: 8px; }}
.arm-name {{ font: 500 18px var(--f-mono); }}
.arm-name small {{ font-size: 12px; color: var(--muted); }}
.arm-sub, .arm-note {{ color: var(--muted); font-size: 13px; }}
.pill {{ font: 500 12px/1 var(--f-mono); padding: 5px 8px; border-radius: 99px; background: var(--run-soft); color: var(--run); white-space: nowrap; }}
.pill.finished {{ background: var(--ok-soft); color: var(--ok); }}
.pill.failed, .pill.stopped {{ background: var(--bad-soft); color: var(--bad); }}
.pill.no-data {{ background: var(--line); color: var(--muted); }}
.bar {{ height: 6px; background: var(--line); border-radius: 3px; overflow: hidden; }}
.bar span {{ display: block; height: 100%; background: var(--accent); }}
.kv {{ display: grid; grid-template-columns: repeat(2, minmax(0,1fr)); gap: 6px 16px; margin: 0; }}
.kv div {{ display: flex; justify-content: space-between; gap: 8px; border-bottom: 1px dotted var(--line); }}
.kv dt {{ color: var(--muted); font-size: 13px; }}
.kv dd {{ margin: 0; font: 13px var(--f-mono); font-variant-numeric: tabular-nums; text-align: right; }}
.setup {{ margin: 0; padding-left: 1.2em; display: grid; gap: 6px; }}
figure {{ margin: 0; display: grid; gap: 8px; }}
figure img {{ width: 100%; height: auto; background: #fff; border: 1px solid var(--line); border-radius: 4px; }}
figcaption {{ color: var(--muted); font-size: 14px; max-width: 80ch; }}
.nodata {{ border: 1px dashed var(--line); border-radius: 4px; padding: 40px 16px; text-align: center; color: var(--muted); font: 14px var(--f-mono); }}
.plates {{ display: grid; gap: 28px; }}
.tablebox {{ overflow-x: auto; }}
table {{ border-collapse: collapse; font-size: 14px; min-width: 560px; width: 100%; }}
th, td {{ padding: 8px 10px; border-bottom: 1px solid var(--line); text-align: right; vertical-align: top; }}
th[scope="row"] {{ text-align: left; font-weight: 400; }}
thead th {{ font: 500 13px var(--f-mono); color: var(--muted); }}
td {{ font-family: var(--f-mono); font-variant-numeric: tabular-nums; }}
td.t {{ font-family: var(--f-body); text-align: left; font-size: 14px; }}
td.w {{ color: var(--muted); font-size: 13px; }}
td.na {{ color: var(--muted); font-family: var(--f-body); font-size: 13px; }}
td.u {{ color: var(--muted); font-family: var(--f-body); text-align: left; font-size: 13px; }}
.caveats {{ margin: 0; padding-left: 1.2em; display: grid; gap: 8px; }}
footer {{ color: var(--muted); font-size: 13px; border-top: 1px solid var(--line); padding-top: 16px; }}
code {{ font-family: var(--f-mono); font-size: .9em; }}
</style>

<div class="wrap">
  <header>
    <div class="eyebrow">ry_per_accretor · env13 · DeltaAI GH200 · half an orbit</div>
    <h1>Plaskett env13 stream run</h1>
    <p class="lede">The L1 stream of Plaskett's star on the resolved envelope of the mass gainer, with the corrected
      env13 setup: an envelope that ends on a Roche equipotential, a Ryu et al. 2025 stream and fluxes measured
      above the bulged photosphere. Compared with the old-setup arms q12_s1 and q12f_s1, whose flux through
      R<sub>acc</sub> was dominated by the breathing of a truncated tidal bulge.</p>
    <div class="{bclass}">{banner}</div>
  </header>

  <section class="status" aria-label="Run status">{status_rows}
  </section>

  <section class="prose">
    <h2>Headline</h2>
    {headline()}
  </section>

  <section>
    <h2>What changed vs q12</h2>
    <div class="tablebox"><table>
      <thead><tr><th></th><th scope="col" style="text-align:left">q12 (env12)</th><th scope="col" style="text-align:left">q13 (env13)</th><th scope="col" style="text-align:left">why</th></tr></thead>
      <tbody>{chg_rows}</tbody>
    </table></div>
  </section>

  <section class="prose">
    <h2>Setup</h2>
    <ul class="setup">
      <li>Plaskett progenitor at the start of Case-A accretion (Wade et al. 2026): M<sub>a</sub> 16 / M<sub>d</sub> 18 M<sub>☉</sub>,
        a = 32.56 R<sub>☉</sub>, P = 3.69 d (P<sub>orb</sub> 0.4583 code), R<sub>acc</sub> = 9.0013 R<sub>☉</sub> (photosphere at φ 90°).</li>
      <li>n = 3 polytropic envelope in the Roche potential (spin 1), c<sub>ph</sub> 19.41 km/s, ρ<sub>ph</sub> = 0.055 ρ<sub>stream</sub>,
        hot ambient 300 km/s; ideal gas, temperatures relaxed on 1e−4 code time (stream effectively isothermal).</li>
      <li>Grid: r = 6.3–13.50 R<sub>☉</sub> (0.85 d<sub>L1</sub>), 640 × 4 × 2048 equatorial wedge, pure plateau radial stretch,
        dr = 2.48 × 10<sup>−3</sup> R<sub>☉</sub> (≤ H<sub>p</sub>/4) over 8.76–9.64.</li>
      <li>Stream: Mdot 1e−4 M<sub>☉</sub>/yr, width 0.735 R<sub>☉</sub> at r<sub>out</sub>, entering at φ ≈ 3.8°. Half an orbit
        (t<sub>lim</sub> 0.2291), bin dumps every P/20.</li>
      <li>Derivation: <code>test1007/setup_derivation/RESULTS.md</code>; input <code>test1007/in/plaskett_env13.athinput</code>.</li>
    </ul>
  </section>

  <section>
    <h2>Envelope breathing</h2>
    <div class="plates">
      {plate('fig_edge_cmp.png', 'envelope edge radius versus phi for each dump, and its displacement in time',
             'Top: radius where the θ-mean density falls to 1 % of ρ<sub>ph</sub>, one curve per bin dump (dark = first, '
             'yellow = last), with the photosphere equipotential (black dashed). The grey band is the stream sector '
             '(φ −35…50°), left out of the numbers. Bottom: displacement of that edge from the first dump, median and '
             'maximum over φ. In q12 the first dump is the spherically truncated IC: the caps blow out and settle about '
             '0.1 R<sub>☉</sub> higher, and the edge keeps moving.', 'no edge data')}
      {plate('fig_envonly.png', 'net mass through the measuring face in the envelope-only tests',
             'Envelope-only tests over 0.05 orbit (no stream; job 3334844): net mass in through the measuring face over '
             'ρ<sub>ph</sub>. t2_sphere = old setup with the spherical top, t2_eq = old setup with the equipotential top, '
             't4_env13 = the env13 envelope. Symlog axis: the sphere top loses ~1 ρ<sub>ph</sub> R<sub>☉</sub><sup>3</sup>; '
             'the equipotential tops stay below 1e−6.', 'no envelope-only test data')}
    </div>
  </section>

  <section>
    <h2>Density maps</h2>
    <div class="plates">
      {map13}
      {plate('fig_q12_s1_maps.png', 'theta-mean density maps of q12_s1',
             '<b>q12_s1</b> (old setup). Same layout; R<sub>acc</sub> is the measuring radius.', 'q12_s1 maps missing')}
    </div>
  </section>

  <section>
    <h2>Mass and angular-momentum budget</h2>
    <div class="plates">
      {bud13}
      {plate('fig_q12_s1_budget.png', 'mass and AM budget of q12_s1',
             '<b>q12_s1.</b> Same layout, measured at R<sub>acc</sub>.', 'q12_s1 budget missing')}
    </div>
  </section>

  <section>
    <h2>Run comparison</h2>
    {plate('fig_plaskett_cmp.png', 'accreted mass, accreted fraction and specific AM through the measuring face for all arms',
           'Left: net mass in through each arm&rsquo;s measuring face (solid) and stream inflow (dashed). Middle: the same over '
           'the final stream inflow. Right: cumulative j of the gas through the face, each in j<sub>Kep</sub> at its own '
           'measuring radius, with the ballistic and rigid-rotation values at each radius.', 'no comparison figure')}
  </section>

  <section>
    <h2>Numbers</h2>
    <div class="tablebox"><table>
      <thead><tr><th></th>{thead}<th scope="col" style="text-align:left">unit</th></tr></thead>
      <tbody>{table}</tbody>
    </table></div>
    {'<p class="lede">* stand-in data, not q13_s1.</p>' if a.standin else ''}
  </section>

  <section class="prose">
    <h2>Caveats</h2>
    <ul class="caveats">
      <li>The measuring face differs: r<sub>meas</sub> 9.616 for q13, R<sub>acc</sub> 9.001 for q12. The q12 face lies inside
        the star on the binary axis, so its flux mixes accretion with the bulge's motion. j values are in j<sub>K</sub> at each
        arm's own radius.</li>
      <li>Masses are in ρ<sub>stream</sub> R<sub>☉</sub><sup>3</sup>, and ρ<sub>stream</sub> is the stream's peak at r<sub>out</sub>.
        The q13 stream is narrower and the density scale was rederived, so compare the accreted <em>fraction</em>, not the
        absolute masses.</li>
      <li>The edge diagnostic uses a 1 % ρ<sub>ph</sub> contour of the θ-mean density, from the bin dumps (P/20 in q13, P/40 in
        q12); breathing faster than the dump cadence is not seen there, but it is in the history fluxes.</li>
      <li>The spikes in the j curves are where the net cumulative mass crosses zero, not features of the flow.</li>
      <li>Half an orbit tests the setup; it is not a steady state.</li>
    </ul>
  </section>

  <footer>Built {now} on DeltaAI with <code>page13/make_all.sh</code> (<code>mkfigs13.py</code>, adapted from viper's
    <code>mkfigs_plaskett.py</code>, and <code>build_page13.py</code>). Runs: <code>{html.escape(a.run)}</code>,
    <code>{RUN12}/q12_s1</code>, <code>q12f_s1</code>.</footer>
</div>
'''
os.makedirs(a.out, exist_ok=True)
fn = os.path.join(a.out, 'index_embedded.html' if a.embed else 'index.html')
open(fn, 'w').write(page)
print('wrote', fn, f'({len(page)/1e6:.2f} MB)', 'state', ARM, st13, '-', why13)
