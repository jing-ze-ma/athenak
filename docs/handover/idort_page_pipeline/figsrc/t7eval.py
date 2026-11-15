"""idort_page Q2: T7 spherical shadow, full-wedge metrics against the exact visible-disk quadrature (s3_eval.shd_exact).
Classes per cell (exact visible fraction f = E_ex/E_free): umbra f < 0.02, lit f > 0.98, cells inside the opaque clump
excluded.  Lit L1 = mean|E/E_ex - 1| at r = 4.05 / 4.99 / 5.77 (all lit cells within 0.36 rad of the clump axis);
umbra <E/E_free> at 4.05; 'front' = the thin cells radially in front of (star side of) an inside-clump cell, exact E/E_free = 1.
Runs: FIN (HFfZ, 120 dirs), PRODd (HFfW), FIN with 432 dirs (Q2 rerun), xs2 FfQ (sptests, 432), and the Q2 one-key
diagnostics (reruns/t7diag/*).   usage: python t7eval.py  -> cache/t7eval.json"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, '/raven/ptmp/jinma/suiteF_1012/sp/ana')
sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python')
os.environ['SPT_RUNS'] = 'runs_fin'
import spl  # noqa: E402
import spt  # noqa: E402
import s3_eval as S  # noqa: E402

RR = '/raven/ptmp/jinma/idort_page/reruns'
# FIN3 page (PG): FIN3 final runs; F2O = the FIN2 control (twin_state mean) of the same binary; old runs in .bak_fin3
import glob
RUNS = {'FIN120': '/raven/ptmp/jinma/idort1020/page/p19runs/sp/shd/HFfZ_32',   # idort-1019 PROD
        'FIN432': '/raven/ptmp/jinma/idort1020/runs/xtra/shd432/p',
        'F2O120': '/raven/ptmp/jinma/idort_page/fin3runs/sp/shd/HFfZ_32'}       # v10 FIN3 (before)
CACHE = '/raven/ptmp/jinma/idort1020/page/cache'


def geometry():
    c = spt.case_cfg('shd')
    pr = c['prob']
    cr, ct, cp, ca = (float(pr['clump_r']), float(pr['clump_th']), float(pr['clump_ph']), float(pr['clump_a']))
    cen = cr*np.array([math.sin(ct)*math.cos(cp), math.sin(ct)*math.sin(cp), math.cos(ct)])
    n1, n2, n3 = 32, 48, 48
    rc = spl.rcent(spl.rfaces(c['r'][0], c['r'][1], n1, c['poly']))
    tc = c['th'][0] + (np.arange(n2) + 0.5)*(c['th'][1] - c['th'][0])/n2
    pc = c['ph'][0] + (np.arange(n3) + 0.5)*(c['ph'][1] - c['ph'][0])/n3
    fn = os.path.join(CACHE, 't7ex.npz')
    if os.path.exists(fn):
        z = np.load(fn)
        return rc, z['EX'], z['FR'], z['ANG'], z['INS']
    EX = np.zeros((n3, n2, n1))
    FR = np.zeros_like(EX)
    ANG = np.zeros((n3, n2))
    INS = np.zeros_like(EX, bool)
    for k in range(n3):
        for j in range(n2):
            u = np.array([math.sin(tc[j])*math.cos(pc[k]), math.sin(tc[j])*math.sin(pc[k]), math.cos(tc[j])])
            ANG[k, j] = math.acos(min(u @ (cen/cr), 1.0))
            for i in range(n1):
                if ANG[k, j] > 0.36:
                    FR[k, j, i] = EX[k, j, i] = 2*(1 - math.sqrt(1 - 1/rc[i]**2))
                    continue
                p = rc[i]*u
                EX[k, j, i], FR[k, j, i] = S.shd_exact(p, cen, ca, ns=30)
                INS[k, j, i] = np.linalg.norm(p - cen) < ca
    np.savez(fn, EX=EX, FR=FR, ANG=ANG, INS=INS)
    return rc, EX, FR, ANG, INS


def main():
    rc, EX, FR, ANG, INS = geometry()
    zx = EX/FR
    lit = (zx > 0.98) & ~INS
    front = np.zeros_like(INS)
    front[:, :, :-1] = INS[:, :, 1:] & ~INS[:, :, :-1]
    out = {}
    for nm, d in RUNS.items():
        if not os.path.isdir(os.path.join(d, 'bin')):
            continue
        fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if '.m1.' in f)
        if not fs:
            continue
        b = spl.binall(os.path.join(d, 'bin', fs[-1]))
        if 'm1_e' not in b:
            continue
        z = b['m1_e']/FR
        o = {'t': float(b['time']), 'lit_L1': {}, 'lit_max': {}}
        for rr in (4.0, 5.0, 5.8):
            i = int(np.argmin(np.abs(rc - rr)))
            m = lit[:, :, i] & (ANG < 0.36)
            o['lit_L1']['%.2f' % rc[i]] = float(np.mean(np.abs(z[:, :, i][m] - 1)))
            o['lit_max']['%.2f' % rc[i]] = float(z[:, :, i][m].max())
        i = int(np.argmin(np.abs(rc - 4.0)))
        um = (zx[:, :, i] < 0.02) & ~INS[:, :, i]
        o['umbra'] = float(z[:, :, i][um].mean())
        o['front_mean'] = float(z[front].mean())
        o['front_max'] = float(z[front].max())
        o['inside_mean'] = float(z[INS].mean())
        out[nm] = o
        print('%-22s lit L1 %s  max %s  umbra %.3f  front %.2f / %.2f  inside %.3f' % (
            nm, ' / '.join('%.3f' % v for v in o['lit_L1'].values()), ' / '.join('%.2f' % v for v in o['lit_max'].values()),
            o['umbra'], o['front_mean'], o['front_max'], o['inside_mean']))
    json.dump(out, open(os.path.join(CACHE, 't7eval.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
