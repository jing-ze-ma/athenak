"""Tables from res_{lowP,old}.pkl -> tables.md (stdout too)."""
import pickle
import sys
import numpy as np
sys.dont_write_bytecode = True
OUT = '/orion/ptmp/jinma/rsg_wind_1008/nlte/'
RSUN = 6.957e10
STAR = {'golden16': dict(R=669.875*RSUN, vcon=7.46e5), 'm20lgl5.5': dict(R=1311.634*RSUN)}


def fmt_roots(rl, tth=False):
    s = [r for r in rl if r['stable']]
    if not s:
        return 'none'
    out = []
    for r in sorted(s, key=lambda r: -r['T']):
        t = f"{r['T']:.0f}"
        if tth:
            t += f"({r['tth']:.1e})"
        out.append(t)
    return '/'.join(out)


def summary(res, tab, lines):
    grp = {}
    for r in res:
        grp.setdefault((r['star'], r['key'], r['val'], r['x']), []).append(r)
    lines.append(f'\n## table = {tab}\n')
    lines.append('stable roots [K] warm/cool; "range" = min-max over q x A_IR x A_el; '
                 'fc = continuum share of the absorbed stellar power at the warm(first) root\n')
    lines.append('| star | constraint | r/R | LTE | (a) scatt. el. | (b) f=0.5 | (b\') f=0.1 | (b\') f=0.5 '
                 '| bistable: LTE/a/b0.5/bp0.1/bp0.5 | fc (a) |')
    lines.append('|---|---|---|---|---|---|---|---|---|---|')
    for (s, key, val, x), rr in grp.items():
        d = {r['scen']: r for r in rr}
        lte = fmt_roots(d['lte']['roots'])

        cols, bis = [], [str(int(sum(q['stable'] for q in d['lte']['roots']) >= 2))]
        fcs = ''
        for p in ('a', 'b0.5', 'bp0.1', 'bp0.5'):
            sel = [r for r in rr if r['scen'].startswith(p + ' ')]
            ws, cs = [], []
            for r in sel:
                st = sorted([q['T'] for q in r['roots'] if q['stable']])
                if len(st) >= 2:
                    ws.append(st[-1])
                    cs.append(st[0])
                elif st:
                    (ws if st[0] > 1500 else cs).append(st[0])
            c = []
            if ws:
                c.append(f'w {min(ws):.0f}-{max(ws):.0f}')
            if cs:
                c.append(f'c {min(cs):.0f}-{max(cs):.0f}')
            cols.append(', '.join(c) or 'none')
            bis.append(f"{sum(1 for r in sel if sum(q['stable'] for q in r['roots']) >= 2)}"
                       f"/{len(sel)}")
            if p == 'a':
                f = [r['roots'][0]['fcont'] for r in sel if r['roots']]
                fcs = f'{min(f):.2f}-{max(f):.2f}' if f else '-'
        v = f'P={val:.0e} bar' if key == 'P' else f'rho={val:.0e}'
        lines.append(f'| {s} | {v} | {x:.0f} | {lte} | ' + ' | '.join(cols) +
                     f" | {' '.join(bis)} | {fcs} |")


def tth_table(res, tab, lines):
    lines.append(f'\n## t_th at the stable roots, table = {tab}, isobaric (P) only\n')
    lines.append('t_flow = r/v (v = 10 / 30 km/s); t_sh = R/v_con (golden16 v_con 7.46 km/s from '
                 'col_golden16_ft12_clamp.npz; m20lgl5.5 FT v_con from rsglib.ft_vcon)\n')
    lines.append('| star | P bar | r/R | t_flow 10/30 | t_sh | LTE T(t_th s) | a: q=1e-12 AIR100 Ael1e8 '
                 '(slowest) | a: q=1e-10 AIR10 Ael1e6 (fastest) | bp0.5 q=1e-12 AIR100 Ael1e8 |')
    lines.append('|---|---|---|---|---|---|---|---|---|')
    for r0 in res:
        if r0['key'] != 'P' or r0['scen'] != 'lte':
            continue
        s, x, val = r0['star'], r0['x'], r0['val']
        R = STAR[s]['R']
        tf = (x*R/1e6, x*R/3e6)
        tsh = R/STAR[s]['vcon']

        def g(name):
            for r in res:
                if (r['star'], r['x'], r['key'], r['val'], r['scen']) == (s, x, 'P', val, name):
                    return fmt_roots(r['roots'], True)
        lines.append(f'| {s} | {val:.0e} | {x:.0f} | {tf[0]:.1e}/{tf[1]:.1e} | {tsh:.1e} | '
                     f"{fmt_roots(r0['roots'], True)} | {g('a q1e-12 AIR100 Ael1e+08')} | "
                     f"{g('a q1e-10 AIR10 Ael1e+06')} | {g('bp0.5 q1e-12 AIR100 Ael1e+08')} |")


def h2only(res, tab, lines):
    lines.append(f'\n## n_H2-only collider flag (q=1e-10, A_IR=10, A_el=1e6, (a)), table = {tab}\n')
    lines.append('| star | constraint | r/R | stable T (t_th s) | eps_IR at root |')
    lines.append('|---|---|---|---|---|')
    for r in res:
        if r['scen'].startswith('H2only') and r['key'] == 'P':
            e = ','.join(f"{q['eir']:.1e}" for q in r['roots'] if q['stable'])
            lines.append(f"| {r['star']} | P={r['val']:.0e} | {r['x']:.0f} | "
                         f"{fmt_roots(r['roots'], True)} | {e} |")


if __name__ == '__main__':
    sys.path.insert(0, OUT)
    import os
    os.environ['CK_DATA'] = '/orion/ptmp/jinma/rsg_wind_1008/lowP/data/'
    import rsglib as L
    m = L.read_mesa('m20lgl5.5')
    STAR['m20lgl5.5']['vcon'] = float(L.ft_vcon(m)[0])
    lines = [f"m20lgl5.5 FT v_con = {STAR['m20lgl5.5']['vcon']:.3e} cm/s"]
    for tab in ('lowP', 'old'):
        res = pickle.load(open(OUT + f'res_{tab}.pkl', 'rb'))
        summary(res, tab, lines)
        tth_table(res, tab, lines)
        h2only(res, tab, lines)
    open(OUT + 'tables.md', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
