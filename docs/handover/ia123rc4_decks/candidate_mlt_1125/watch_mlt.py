#!/usr/bin/env python3
"""watch_mlt.py RUNDIR [--once]   T9 MLT full run: one status line per call; writes RUNDIR/ALERT_MLT (reason) and returns 3 when a STOP criterion of the NOTE is met.  Read only otherwise.
Criteria (fixed before the run): (a) CLEAN WARNING, FATAL, inf or nan in the run log / hst; (b) L_top/L_in < 0.8 or > 3 continuously for 2e4 s (t > 2e4 s);
(c) dt < 20 s for 200 consecutive cycles.  Stopping itself stays with the launcher's STOP mechanism."""
import sys, re, os, glob, numpy as np
R = sys.argv[1]; log = sorted(glob.glob(R + '/run.log*') or glob.glob(R + '/log.*'), key=os.path.getmtime); hst = glob.glob(R + '/*.user.hst')
alert = []; cnt = {}
for f in log:
    # the launcher's summary lines ('... CW 0 FATAL 0 inf 0 ...') are counts, not events
    t = ''.join(l for l in open(f, errors='replace') if not re.search(r'CW \d+ FATAL \d+', l))
    for k, p in (('CW', r'CLEAN WARN'), ('FATAL', r'FATAL'), ('DT_RETRY', r'DT_RETRY (?!on)'), ('NaN', r'(?<![A-Za-z])-?nan(?![A-Za-z])'), ('inf', r'(?<![A-Za-z_])-?inf(?![A-Za-z])')):
        cnt[k] = cnt.get(k, 0) + len(re.findall(p, t))
    m = re.findall(r'backward-Euler steps=([0-9.eE+-]+)', t); cnt['BE'] = float(m[-1]) if m else 0
    m = re.findall(r'stage fallbacks=([0-9.eE+-]+)', t); cnt['fallback'] = float(m[-1]) if m else 0
    cnt['clamped'] = cnt.get('clamped', 0) + len(re.findall(r'clamped [1-9]', t))
if cnt.get('CW'): alert.append('CLEAN WARNING %d' % cnt['CW'])
if cnt.get('FATAL'): alert.append('FATAL')
if cnt.get('inf') or cnt.get('NaN'): alert.append('inf or nan in the log')
line = 'no hst'
if hst:
    nm = re.findall(r'\]=(\S+)', [l for l in open(hst[0]) if '[1]' in l][0]); rows = []
    for l in open(hst[0]):
        if l[0] != '#' and len(l.split()) == len(nm): rows.append([float(x) for x in l.split()])
    d = np.array(rows); _, i = np.unique(d[:, 0], return_index=True); d = d[i]; h = {n: d[:, j] for j, n in enumerate(nm)}
    if not np.isfinite(d).all(): alert.append('non-finite hst value')
    t, L, dt = h['time'], h['L_top']/h['L_in'], h['dt']
    bad = ((L < 0.8) | (L > 3.0)) & (t > 2.0e4); k = 0.0; t0 = None
    for a, b in zip(t, bad):
        if b and t0 is None: t0 = a
        if not b: t0 = None
        if b and a - t0 >= 2.0e4: alert.append('L_top/L_in outside 0.8-3 for 2e4 s from t %.4e (now %.2f)' % (t0, L[-1])); break
    run = 0
    for x in (dt[1:] < 20.0):
        run = run + 1 if x else 0
        if run >= 200: alert.append('dt < 20 s for 200 rows'); break
    cyc = len(t) - 1; es = h['P_esrc'][-1]/h['L_in'][-1] if 'P_esrc' in h else float('nan')
    line = 't %.5e rows %d dt %.2f L_top/L_in %.4f L_bot/L_in %.4f KE %.3e w_mlt %.3f L_MLT/L_in %.3f sum_esrc/L_in %.1e' % (t[-1], len(t), dt[-1], L[-1], h['L_bot'][-1]/h['L_in'][-1], h['KE_int'][-1], h.get('w_mlt', [float('nan')])[-1], h.get('L_MLT', [float('nan')])[-1]/h['L_in'][-1], es)
print(line, '|', ' '.join('%s %g' % kv for kv in sorted(cnt.items())))
if alert:
    open(R + '/ALERT_MLT', 'a').write('; '.join(alert) + ' :: ' + line + '\n'); print('STOP CRITERION:', '; '.join(alert)); sys.exit(3)
