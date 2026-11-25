"""ana.py: R_ph, T_eff of the restart arms against the production dumps (out/*.json, out/*_w.json); exact times from user.hst by nearest row"""
import json, re, glob, os, sys
import numpy as np
J = '/resnick/groups/carnegie_poc/jingze/'
PD = dict(m35zams='prod_c1_1124', m35mms='prod_c1_1124', bsg='prod_c1_1124', t9l62='prod_c4_1124', t9mlt='prod_mlt_1125')
RST = dict(m35zams=50816, m35mms=22144, bsg=29127, t9l62=6332, t9mlt=None)
def hst(f):
    a = [l.split()[:1] for l in open(f) if not l.startswith('#') and l.strip()]; t = np.array(a, float)[:, 0]; return t
for star in sys.argv[1:]:
    o = json.load(open(J + 'rcont_1124/out/%s.json' % star)); w = {(a, f): s for a, f, s in json.load(open(J + 'rcont_1124/out/%s_w.json' % star))}
    cache = {}
    for x in o:
        s = w[(x['tag'], x['f'])]
        g = lambda p: float(re.search(p, s).group(1))
        x['O23'] = g(r'R23 med \S+ mean (\S+)'); x['Osh'] = g(r'Rsh23 (\S+)'); x['OT'] = g(r'T_eff\(R23mean\) (\S+)'); x['Lhst'] = g(r' L (\S+)')
        d = os.path.dirname(os.path.dirname(x['f'])); h = glob.glob(d + '/*.user.hst')[0]
        if h not in cache: cache[h] = hst(h)
        x['tx'] = float(cache[h][np.argmin(abs(cache[h] - x['t']))]) if x['cyc'] > 0 else x['t']
    P = [x for x in o if x['tag'] == 'P']; Q = [x for x in o if x['tag'] == 'Q']; A = {x['tag']: x for x in o if x['tag'] not in 'PQ'}
    P.sort(key=lambda x: x['t']); Q.sort(key=lambda x: x['t'])
    print('==== %s dtype %s; production dump before: %s c%d t %.6e; R1 c%d, R50 c%d t %.6e; next production dump c%d t %.6e' % (star, P[-1]['dtype'], os.path.basename(P[-1]['f']), P[-1]['cyc'], P[-1]['tx'], A['R1']['cyc'], A['R50']['cyc'], A['R50']['tx'], Q[0]['cyc'], Q[0]['tx']))
    for k, nm in (('RA', 'R(mean tau=2/3)'), ('RA1', 'R(mean tau=1)'), ('RM', 'R col median'), ('TA', 'T(RA flux)'), ('TM', 'T col median'), ('O23', 'rph.py R23_mean'), ('Osh', 'rph.py Rsh23'), ('OT', 'rph.py T_eff')):
        p = np.array([x[k] for x in P], float); tp = np.array([x['tx'] for x in P]); q = [x[k] for x in Q]; tq = [x['tx'] for x in Q]
        dd = np.diff(p) / p[:-1]; n = len(dd); last = dd[-min(20, n):]
        # noise of a 4-point cubic interpolation (leave-one-out on the production dumps incl. 2 after)
        T = np.concatenate([tp, tq]); V = np.concatenate([p, q]); res = []
        for i in range(2, len(T) - 2):
            m = [i - 2, i - 1, i + 1, i + 2]; c = np.polyfit(T[m] - T[i], V[m], 3); res.append((np.polyval(c, 0) - V[i]) / V[i])
        rms = np.sqrt(np.mean(np.array(res) ** 2)) if res else np.nan
        line = '%-17s' % nm + ' P %.8e' % p[-1]
        for a in ('R1', 'R25a', 'R50', 'C50'):
            if a in A: line += ' | %s-P %+.3e' % (a, A[a][k] / p[-1] - 1)
        line += ' || R50-C50 %+.2e' % (A['R50'][k] / A['C50'][k] - 1) if 'C50' in A else ' || C50 missing'
        line += ' || prod dump-to-dump rms %.1e max %.1e (n %d)' % (np.sqrt(np.mean(last ** 2)), np.max(abs(last)), len(last))
        m = [len(tp) - 2, len(tp) - 1]; m4 = [len(tp) - 2, len(tp) - 1] + ([len(tp), len(tp) + 1] if len(q) > 1 else [len(tp)])
        for a in ('R1', 'R50', 'C50'):
            if a in A:
                tt = A[a]['tx']; deg = len(m4) - 1; c = np.polyfit(T[m4] - tt, V[m4], deg); line += ' | %s vs prod cubic %+.2e' % (a, A[a][k] / np.polyval(c, 0) - 1)
        line += ' (loo rms of the cubic at full spacing %.1e)' % rms
        print(line)
