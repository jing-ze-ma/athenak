"""littests-1011 shared readers"""
import glob
import os
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/blend_1022/wt/vis/python')
import bin_convert as bc  # noqa: E402
W = "/raven/ptmp/jinma/order_gpu_1110/beams"


def last(d, var):
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.%s.*.bin' % var)))
    if not fs:
        fs = sorted(glob.glob(os.path.join(d, '*.%s.*.bin' % var)))
    return fs[-1] if fs else None


def load(d):
    """E (nx2, nx1), J_sc, F1, F2, x1v, x2v from the last dumps of run dir d"""
    f = last(d, 'm1')
    v = bc.read_binary_as_athdf(f)
    out = dict(E=v['m1_e'][0], F1=v['m1_f1'][0], F2=v['m1_f2'][0],
               x1=v['x1v'], x2=v['x2v'], t=v.get('Time', None))
    g = last(d, 'm1_vsc')
    if g:
        out['J'] = bc.read_binary_as_athdf(g)['m1_vj'][0]
    return out


def runinfo(d):
    r = {}
    try:
        r['rc'] = int(open(d + '/rc.txt').read().split()[0])
        r['sec'] = float(open(d + '/time.txt').read().split()[-2])
    except Exception:
        pass
    for ln in open(d + '/run.log', errors='ignore'):
        if 'Picard iterations mean' in ln:
            p = ln.split()
            r['picard'] = float(p[p.index('mean=') + 1]) if 'mean=' in p else None
            try:
                r['picard'] = float(ln.split('mean=')[1].split()[0])
                r['nonconv'] = float(ln.split('NON-CONVERGED=')[1].split()[0])
            except Exception:
                pass
    return r
