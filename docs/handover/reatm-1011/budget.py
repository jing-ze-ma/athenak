"""budget.py run...: energy budget from agcar3d.user.hst: R = Etot - Etot0 - int(L_in - L_top + P_esrc - P_sponge) dt
(trapezoid over hst rows, dcycle 1), relative to int|L_in| dt; also L_top/L_top0, Picard mean, Mdot_top, Min_top."""
import sys, re, numpy as np
for d in sys.argv[1:]:
    f = d + '/agcar3d.user.hst'; x = np.loadtxt(f); hdr = [l for l in open(f) if l.startswith('#  [')][0]
    H = {n: x[:, i] for i, n in enumerate(re.findall(r'\[\d+\]=(\S+)', hdr))}
    t = H['time']; s = H['L_in'] - H['L_top'] + H.get('P_esrc', 0*t) - H.get('P_sponge', 0*t)
    I = np.concatenate([[0], np.cumsum(0.5*(s[1:] + s[:-1])*np.diff(t))])
    R = H['Etot'] - H['Etot'][0] - I; S = np.concatenate([[0], np.cumsum(0.5*np.abs(H['L_in'][1:] + H['L_in'][:-1])*np.diff(t))])
    n = len(t) - 1; idx = [i for i in (10, 30, 100, 200, n) if i <= n]
    print('%-8s n=%d  R/int L_in at rows %s: %s | dEtot/E0 %.2e | L_top/L0 %s | Picard mean %.3f | Mdot_top max %.3g Min_top end %.3g' % (
        d.rstrip('/').split('/')[-1], n, idx, ' '.join('%.2e' % (R[i]/S[i]) for i in idx), (H['Etot'][-1] - H['Etot'][0])/H['Etot'][0],
        ' '.join('%.5f' % (H['L_top'][i]/H['L_top'][0]) for i in idx), H['Picard'][1:].mean(), H['Mdot_top'].max(), H['Min_top'][-1]))
