#!/usr/bin/env python3
"""blend-taufree-1022 phase 2: score IDORT runs (arm x metric + robustness/cost counters).
usage: score.py ROOT [deck ...]   ROOT/<deck>/<arm>/rep<r>/ (run_set.sh layout)"""
import glob
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/sp')


def counters(d):
    o = dict(rc='?', cyc=-1, t=np.nan, dt=np.nan, nadm=0, be=0, dtr=0, pic=np.nan, nc=np.nan,
             wall=np.nan, fatal='')
    try:
        o['rc'] = open(d + '/rc.txt').read().split('=')[1].strip()
        o['wall'] = float(open(d + '/time.txt').read().split()[0])
    except Exception:
        pass
    try:
        txt = open(d + '/run.log', errors='ignore').read()
    except Exception:
        return o
    cl = re.findall(r'cycle=(\d+) time=([-+0-9.eE]+) dt=([-+0-9.eE]+)', txt)
    if cl:
        o['cyc'], o['t'], o['dt'] = int(cl[-1][0]), float(cl[-1][1]), float(cl[-1][2])
    m = re.findall(r'time=([-+0-9.eE]+) cycle=(\d+)\n', txt)
    if m:
        o['t'], o['cyc'] = float(m[-1][0]), int(m[-1][1])
    o['nadm'] = txt.count('NOT ADMISSIBLE')
    m = re.search(r'(\d+) cycles took the backward-Euler redo', txt)
    o['be'] = int(m.group(1)) if m else 0
    o['dtr'] = len(re.findall(r'^DT_RETRY cycle', txt, re.M))
    m = re.search(r'Picard iterations mean=([-+0-9.eE]+) max=\S+ NON-CONVERGED=([-+0-9.eE]+)', txt)
    if m:
        o['pic'], o['nc'] = float(m.group(1)), float(m.group(2))
    for ln in txt.splitlines():
        if 'FATAL' in ln or 'DIVERGED' in ln:
            o['fatal'] = ln.strip()[:120]
            break
    return o


def metric(deck, d):
    try:
        if deck.startswith(('cb', 'wide', 'wg', 'ws')):
            import cb as C
            snap = re.findall(r'implicit_beam(\d): .*in-plane ang=([-+0-9.eE]+) deg', open(d + '/run.log').read())
            if not hasattr(C, 'readin0'):
                C.readin0 = C.readin
            def rd(fn, snap=snap):
                p = C.readin0(fn)
                for b, a in snap:
                    p['rad_m1/implicit_beam%s_ang' % b] = a
                return p
            C.readin = rd
            m = C.metrics(d)
            return 'E_L1 %.4f lum_top %.4f peak %s fwhm %s (ex %s)' % (
                m['E_L1'], m['E_lum_top'], '/'.join('%.3f' % x for x in m['E_peak']),
                '/'.join('%.1f' % x for x in m['E_fwhm']), '/'.join('%.1f' % x for x in m['ex_fwhm'])), m['E_L1']
        if deck.startswith('shd'):
            import shd as S
            m = S.metrics(d)
            return 'E_L1 %.4f umbra0.25 %.3e (ex %.3e) umbra0.40 %.3e (ex %.3e)' % (
                m['E_L1'], m['E_umbra_0.25'], m['ex_umbra_0.25'], m['E_umbra_0.40'],
                m['ex_umbra_0.40']), m['E_L1']
        if deck.startswith('hcyl'):
            import hcyl as H
            m, _ = H.metrics(d)
            return 'E_L1_in %.4f E_L1_out %.4f F_L1_out %.4f rF %.4f' % (
                m['E_L1_in'], m['E_L1_out'], m['F_L1_out'], m['rF_code/rF_ex(2-3.5)']), m['E_L1_out']
        if deck.startswith('nltek'):
            from s2 import nlte
            m = nlte(d)
            return 'eps %.0e L1 %.4f max %.4f top %.4f' % (m['eps'], m['L1'], m['max'], m['top']), m['L1']
        if deck.startswith('core'):
            import spl
            from s3_eval import core_exact
            t = spl.tab(sorted(glob.glob(d + '/**/*.m1.*.tab', recursive=True))[-1])
            r = t['x1v']
            J, Hh = core_exact(r)
            rel = t['m1_e']/J - 1.0
            fr = t['m1_f1']*r*r
            c = float(open(d + '/in.athinput').read().split('c_light')[1].split('=')[1].split()[0])
            out = 'L1(r>1.1) %.4f E/J-1: r<0.9 %+.3e 0.9-1.1 %+.3e 1.1-3 %+.3e r>3 %+.3e | Fr2/Hr2(r>3) %.4f' % (
                np.mean(np.abs(rel[r > 1.1])), np.mean(rel[r < 0.9]), np.mean(rel[(r >= .9) & (r < 1.1)]),
                np.mean(rel[(r >= 1.1) & (r < 3)]), np.mean(rel[r >= 3]),
                np.mean(fr[r > 3])/c/(Hh[-1]*r[-1]**2))
            return out, np.mean(np.abs(rel[r > 1.1]))
    except Exception as e:
        return 'ERR %s' % e, np.nan
    return '', np.nan


def main():
    root = sys.argv[1]
    decks = sys.argv[2:] or sorted(os.listdir(root))
    print('| deck | arm | rep | rc | cycles | t | dt | NOT ADMISSIBLE | BE redo cycles | DT_RETRY | Picard mean | NC | wall s | s/cycle | metric |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for dk in decks:
        for d in sorted(glob.glob(os.path.join(root, dk, '*', 'rep*'))):
            arm, rep = d.split('/')[-2], d.split('/')[-1]
            c = counters(d)
            txt, _ = metric(dk, d) if c['rc'] == '0' else (c['fatal'], np.nan)
            spc = c['wall']/c['cyc'] if c['cyc'] > 0 else np.nan
            print('| %s | %s | %s | %s | %d | %.4g | %.3g | %d | %d | %d | %.2f | %g | %.1f | %.3f | %s |' % (
                dk, arm, rep, c['rc'], c['cyc'], c['t'], c['dt'], c['nadm'], c['be'], c['dtr'], c['pic'],
                c['nc'], c['wall'], spc, txt))


if __name__ == '__main__':
    main()
