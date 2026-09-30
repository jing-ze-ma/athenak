"""table.py: RESULTS.md numbers from out/*.npz (ana.py measurements, LSF R = 70000).
RVs are Gaussian-fit centres minus the windless template's own centre (zero point)."""
import numpy as np
from ana import analyse

W = '/viper/ptmp2/jinma/w121prod_0929/synth2_trans/out/'


def rv(name, R=70000.0):
    res, cf, _ = analyse(name, R, ret=True)
    run = name.split('_')[0]
    t, _, _ = analyse(run + '_template', R, ret=True)
    z0 = t[('all', 'both')]
    o = {}
    for k, v in res.items():
        o[k] = dict(fe=v['fe'] - z0['fe'], na=np.mean(v['na']) - np.mean(z0['na']),
                    fe_amp=v['fe_amp'], na_amp=np.mean(v['na_amp']), fe_w=v['fe_w'])
    return o, cf


def ew(name):
    z = np.load(W + name + '.npz')
    mw = np.load(W + name.split('_')[0] + '_template.npz')['mw']
    fe = np.einsum('w,elwu->l', mw, z['Afe'] - z['Afe_c'][..., None])
    na = (z['Ana'] - z['Ana_c'][..., None]).sum((0, 2, 3))
    return fe, na


if __name__ == '__main__':
    import sys
    names = sys.argv[1:]
    for n in names:
        try:
            o, cf = rv(n)
        except FileNotFoundError:
            continue
        print('## %s' % n)
        for k in [('all', 'morning'), ('all', 'evening'), ('all', 'both'),
                  ('first half', 'both'), ('second half', 'both')]:
            v = o[k]
            print('%-12s %-8s Fe %+6.2f (sigma %.1f, amp %.0f ppm)  NaD %+6.2f (amp %.0f ppm)' % (
                k + (v['fe'], v['fe_w'], v['fe_amp']*1e6, v['na'], v['na_amp']*1e6)))
        for k, v in cf.items():
            print('CF %-2s %-8s logp16/50/84 %6.2f %6.2f %6.2f  above x1max %.3f  in old window %.3f'
                  % (k + tuple(v)))
        fe, na = ew(n)
        print('EW-proxy Fe morning/evening %.4g %.4g ; Na %.4g %.4g' % (fe[0], fe[1], na[0], na[1]))
