#!/usr/bin/env python3
"""picstop-1017: per-arm table from run directories (run.log, rc.txt).  usage: stat.py DIR [DIR ...]
Tail = solves whose ACCEPTED iterate (final PICRES check of the solve: stop=acc/exh or no stop field) has max(rc,rT) > 1e-4 / > 1e-2
(cells above e_floor).  Passes = 'Picard iterations mean' of the run summary (real passes; a residual check that fails becomes a
pass and is counted there).  checks = residual evaluations (shadow passes)."""
import sys, re, os
def num(pat, s, d='nan'):
    m = re.findall(pat, s, re.M); return m[-1] if m else d
print('%-12s %6s %6s %7s %5s %4s %5s %8s %8s %8s %7s %7s %8s %8s %7s %s' % ('arm', 'rc', 'cyc', 'solves', 'pass', 'NC', 'retry',
      'tail>1e-4', 'tail>1e-2', 'max_acc', 'checks', 'fails', 'exh(NC)', 'ce_refus', 'wall_s', 't_end'))
for d in sys.argv[1:]:
    nm = os.path.basename(d.rstrip('/'))
    try: L = open(d + '/run.log', errors='replace').read()
    except OSError: print(nm, 'no run.log'); continue
    rc = open(d + '/rc.txt').read().strip() if os.path.exists(d + '/rc.txt') else 'rc=?'
    wall = num(r'wall=([0-9.]+)', rc); rcv = num(r'rc=(\d+)', rc, '?')
    n = 0; t4 = 0; t2 = 0; mx = 0.0
    for l in L.splitlines():
        if not l.startswith('PICRES'): continue
        if ' stop=cont' in l: continue
        d0 = dict(t.split('=', 1) for t in l.split()[1:] if '=' in t and not t.startswith('n='))
        rcx = float(re.search(r' rc=([^ ]+)', l).group(1)); rTx = float(re.search(r' rT=([^ ]+)', l).group(1))
        w = max(rcx, rTx); n += 1; t4 += w > 1e-4; t2 += w > 1e-2; mx = max(mx, w)
    solves = num(r'implicit transport: solves=([^ ]+)', L); pas = num(r'Picard iterations mean=([^ ]+)', L)
    nc = num(r'NON-CONVERGED=([^ \n]+)', L); retr = len(re.findall(r'^DT_RETRY cycle', L, re.M))
    ps = re.findall(r'picstop: .*', L); ps = ps[-1] if ps else ''
    chk = num(r'checks=([^ ]+)', ps, '-'); fl = num(r'failed\+continued=([^ ]+)', ps, '-'); ex = num(r'NC\)=([^ ]+)', ps, '-')
    ref = num(r'implicit_conv_est_nofloor=\d refused=([^ ]+)', ps, '-'); refk = num(r'implicit_conv_est_guard=[^ ]+ refused=([^ ]+)', ps, '-')
    cyc = num(r'^time=[^ ]+ cycle=(\d+)', L, '?'); te = num(r'^time=([^ ]+) cycle=', L, '?')
    if cyc == '?': cyc = num(r'cycle=(\d+)', L, '?')
    f = lambda a: ('%5.2f%%' % (100.0 * a / n)) if n else '-'
    print('%-12s %6s %6s %7.0f %5.2f %4.0f %5d %8s %8s %8.1e %7s %7s %8s %8s %7s %s' % (nm, rcv, cyc, float(solves), float(pas), float(nc),
          retr, f(t4), f(t2), mx, chk, fl, ex, '%s/%s' % (ref, refk), wall, te))
