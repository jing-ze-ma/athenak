# amlz-1112: L budgets of the small-HG arms (AMH_CONS, AMCS_CONS, AMCS_DIAG) from <arm>/run.log
import re, sys
import numpy as np
fl = r'([-+0-9.eE]+)'
G = ''  # usage: amlz_analyze.py RUNDIR ... (reads RUNDIR/run.log)
for a in sys.argv[1:]:
    t = open(G + a + '/run.log').read()
    h = [list(map(float, m)) for m in re.findall(
        r'sum_dT=\(%s %s %s\).*?resid=\(%s %s %s\) scale=%s \|dep\|=%s Lfluid=\(%s %s %s\)' % ((fl,)*11), t)]
    if h:
        print(a, 'AMH stages %d sum_dT_z %.3e max|resid_z| %.3e Lz %.6e -> %.6e' % (
            len(h), sum(r[2] for r in h), max(abs(r[5]) for r in h), h[0][10], h[-1][10]))
    c = [list(map(float, m)) for m in re.findall(
        r'AMCS_CONS solve \d+ dL_exch %s %s %s dtT %s %s %s resid %s %s %s withheld %s %s %s scale_z %s' % ((fl,)*13), t)]
    if c:
        print(a, 'AMCS solves %d max|resid_z| %.3e max|resid_xy| %.3e mean scale_z %.3e' % (
            len(c), max(abs(r[8]) for r in c), max((r[6]**2 + r[7]**2)**.5 for r in c),
            sum(r[12] for r in c)/len(c)))
    d = re.findall(r'AMCS_DIAG (\d+) %s (.*)' % fl, t)
    if d:
        rows = [list(map(float, x[2].split())) for x in d]
        n = min(len(r) for r in rows)
        s = np.array([r[:n] for r in rows]).sum(0)
        print(a, 'AMCS_DIAG sums |dLr_z| |dLg_z| |dtTpl_z| |dtTrad_z| |del_z| |del_z|/(rhoVR) a_f |Lg_z|:',
              ' '.join('%.3e' % x for x in s))
