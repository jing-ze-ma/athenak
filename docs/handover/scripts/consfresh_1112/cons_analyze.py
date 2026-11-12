#!/usr/bin/env python3
"""consaudit-1112: per-cycle and accumulated conservation budgets from run.log
(EAUDIT / CAUDIT / CAUDITB / CAUDITH lines of <rad_m1>/dbg_energy_audit + dbg_cons_audit).
usage: analyze.py run.log [ncyc] [ccr]
Quantities (all x volume, rank sums): U = IEN + (c/chat) E, M, Lz_g, Lz_r = (r x F/(c chat)).z,
Pr_g = m_r, Pr_r = F_r/(c chat).  Segments (rad_m1.hpp EAuditMark): 0 outside the M1 chain
(hydro stages, pgen sources, hesdirk2 stage average), 1 closure limits, 2 closure exit -> solve,
3 solve, 4 pgen m1post hook, 5 fluid refresh (C2P floors), 6 discarded attempts (put back).
Boundary: CAUDITH = hydro x1 fluxes, rk2-combined (outward +); EAUDIT B = radiation energy rates
of each solve's stored faces (outward +)."""
import re
import sys
from collections import defaultdict

G = 1.0 - 1.0/2**0.5
B2 = G/(2.0*(1.0 - G))


def nums(s):
    return [float(x) for x in s.split(',')]


def parse(fn):
    ea, ca, cb, ch = [], [], [], {}
    for ln in open(fn, errors='replace'):
        if ln.startswith('EAUDIT '):
            d = dict(re.findall(r'(\w+)=(\S+)', ln))
            ea.append(dict(cyc=int(d['cycle']), t2s=int(d['t2s']), dt=float(d['dt']),
                           s=[sum(nums(d['s%d' % q])) for q in range(6)],
                           B=nums(d['B']), T2=float(d['T2']), W=float(d['W']),
                           P=nums(d['P']), RS=float(d['RSOLVE']), D=float(d['D'])))
        elif ln.startswith('CAUDIT '):
            d = dict(re.findall(r'(\w+)=(\S+)', ln))
            ca.append(dict(cyc=int(d['cycle']), t2s=int(d['t2s']), dt=float(d['dt']),
                           s=[nums(d['s%d' % q]) for q in range(7)]))
        elif ln.startswith('CAUDITB '):
            p = ln.split()
            cb.append((int(p[1].split('=')[1]), [float(x) for x in p[2:]]))
        elif ln.startswith('CAUDITH '):
            d = dict(re.findall(r'(\w+)=(\S+)', ln))
            ch[int(d['cycle'])] = dict(st=int(d['st']), dt=float(d['dt']),
                                       i=nums(d['in']), o=nums(d['out']))
    return ea, ca, cb, ch


def main():
    fn = sys.argv[1]
    ncyc = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
    ccr = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    ea, ca, cb, ch = parse(fn)
    cycs = sorted(set(c['cyc'] for c in ca))
    cycs = cycs[1:1 + ncyc]   # the first cycle holds the audit's own start-up
    by_ca = defaultdict(list)
    for c in ca:
        by_ca[c['cyc']].append(c)
    by_ea = defaultdict(list)
    for e in ea:
        by_ea[e['cyc']].append(e)
    acc = defaultdict(float)
    prevB2 = None
    rows = []
    for cy in cycs:
        lines = by_ca[cy]
        eas = by_ea[cy]
        tot = [sum(sum(l['s'][s][a] for s in range(7)) for l in lines) for a in range(6)]
        seg = [[sum(l['s'][s][a] for l in lines) for a in range(6)] for s in range(7)]
        h = ch.get(cy)
        bh = [(h['i'][a] + h['o'][a]) if h else float('nan') for a in range(4)]  # M E Lz Pr
        dt = h['dt'] if h else eas[-1]['dt']
        # radiation energy through the boundary: (A) per-solve dt_sub sums, (B) ESDIRK weights
        bra = sum(e['dt']*ccr*sum(e['B'][:5]) for e in eas)
        s1 = [e for e in eas if e['t2s'] == 2]
        s2 = [e for e in eas if e['t2s'] == 3]
        if s1 and s2 and prevB2 is not None and len(eas) == 2:
            brb = dt*ccr*(0.5*prevB2 + B2*sum(s1[-1]['B'][:5]) + G*sum(s2[-1]['B'][:5]))
        elif eas and eas[-1]['t2s'] in (0, 1):
            brb = dt*ccr*sum(eas[-1]['B'][:5])
        else:
            brb = float('nan')
        prevB2 = sum(s2[-1]['B'][:5]) if s2 else None
        rs = sum(e['RS'] for e in eas)
        T2 = sum(e['T2'] for e in eas)
        W = sum(e['W'] for e in eas)
        Pun = sum(e['P'][0] for e in eas)
        fixE = sum(seg[s][0] for s in (1, 2, 4, 5))
        RE_a = tot[0] + bh[1] + bra
        RE_b = tot[0] + bh[1] + brb
        # L_z: gas + rad; radiation boundary torque not tallied (in the solve segment)
        dLz = tot[2] + tot[3]
        RLz = dLz + bh[2]
        fixLz = sum(seg[s][2] + seg[s][3] for s in (1, 2, 4, 5))
        solLz = seg[3][2] + seg[3][3]
        hydLz = seg[0][2] + seg[0][3] + bh[2]
        dPr = tot[4] + tot[5]
        fixPr = sum(seg[s][4] + seg[s][5] for s in (1, 2, 4, 5))
        hydPr = seg[0][4] + seg[0][5] + bh[3]     # pressure/curvature + gravity + pgen
        fgas = seg[3][4]                          # radiation force impulse on the gas
        solPr = seg[3][4] + seg[3][5]
        dM = tot[1] + bh[0]
        row = dict(cy=cy, dt=dt, dU=tot[0], BhE=bh[1], BrA=bra, BrB=brb, RE_a=RE_a, RE_b=RE_b,
                   RS=rs, fixE=fixE, s0E=seg[0][0] + bh[1], T2=T2, W=W, Pun=Pun,
                   D=seg[6][0], dLz=dLz, BhLz=bh[2], RLz=RLz, fixLz=fixLz, solLz=solLz,
                   hydLz=hydLz, LzG3=seg[3][2], LzR3=seg[3][3], dPr=dPr, BhPr=bh[3],
                   hydPr=hydPr, fgas=fgas, PrR3=seg[3][5], solPr=solPr, fixPr=fixPr, RM=dM)
        rows.append(row)
        for k, v in row.items():
            if k not in ('cy', 'dt'):
                acc[k] += v
                acc['abs_' + k] += abs(v)
    keys = ['dU', 'BhE', 'BrA', 'BrB', 'RE_a', 'RE_b', 'RS', 'fixE', 's0E', 'T2', 'W', 'Pun',
            'D', 'dLz', 'BhLz', 'RLz', 'solLz', 'LzG3', 'LzR3', 'hydLz', 'fixLz',
            'dPr', 'BhPr', 'hydPr', 'fgas', 'PrR3', 'solPr', 'fixPr', 'RM']
    print('file', fn, 'cycles', cycs[0] if cycs else None, '..', cycs[-1] if cycs else None,
          'n', len(cycs))
    print('cy ' + ' '.join('%11s' % k for k in keys))
    for r in rows[:10] + (rows[-3:] if len(rows) > 13 else []):
        print('%d ' % r['cy'] + ' '.join('%11.3e' % r[k] for k in keys))
    print('SUM ' + ' '.join('%11.3e' % acc[k] for k in keys))
    print('ABS ' + ' '.join('%11.3e' % acc['abs_' + k] for k in keys))
    # relative residuals: per cycle max |R| / scale and accumulated |sum R| / sum scale
    def rel(rk, sk):
        pc = max(abs(r[rk])/max(max(abs(r[s]) for s in sk), 1e-300) for r in rows)
        ac = abs(acc[rk])/max(max(acc['abs_' + s] for s in sk), 1e-300)
        return pc, ac
    out = {}
    out['E_a'] = rel('RE_a', ['dU', 'BhE', 'BrA'])
    out['E_b'] = rel('RE_b', ['dU', 'BhE', 'BrB'])
    out['E_fix'] = rel('fixE', ['dU', 'BhE', 'BrB'])
    out['E_rsolve'] = rel('RS', ['dU', 'BhE', 'BrB'])
    out['Lz_tot'] = rel('RLz', ['solLz', 'hydLz', 'LzG3', 'LzR3', 'BhLz'])
    out['Lz_fix'] = rel('fixLz', ['LzG3', 'LzR3', 'BhLz'])
    out['Lz_solve'] = rel('solLz', ['LzG3', 'LzR3', 'BhLz'])
    out['Pr_fix'] = rel('fixPr', ['BhPr', 'fgas', 'PrR3'])
    out['Pr_solve'] = rel('solPr', ['BhPr', 'fgas', 'PrR3'])
    out['M'] = rel('RM', ['dU'])
    for k, (pc, ac) in out.items():
        print('REL %-9s per-cycle max %.3e  accumulated %.3e' % (k, pc, ac))
    # bins: accumulated over the analysed cycles, solve (g0) and fix (g1) of U, Lz, Pr
    bs = [0.0]*144
    for cy, b in cb:
        if cy in cycs:
            bs = [x + y for x, y in zip(bs, b)]
    for g, gn in ((0, 'solve'), (1, 'fix')):
        for q, qn in ((0, 'U'), (1, 'Lz'), (2, 'Pr')):
            v = bs[g*72 + q*24: g*72 + q*24 + 24]
            print('BIN %-5s %-2s ' % (gn, qn) + ' '.join(
                '%9.2e' % v[rb*3 + lb] for rb in range(8) for lb in range(3)))


if __name__ == '__main__':
    main()
