import sys, glob, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/vis/python')
import types; sys.modules.setdefault('h5py', types.ModuleType('h5py'))
import bin_convert as bc
S = '/resnick/groups/carnegie_poc/jingze/smoke_ia121_1123/m35zams128'
files = sorted(glob.glob(S + '/bin/m35zams.hydro_w.*.bin'))
sel = [files[i] for i in sorted(set([0, 2, 5, 10, 15, len(files) - 1]))]
x1 = np.loadtxt(S + '/m35zams.x1grid.txt', comments='#', max_rows=390)[3:387, 2]
print('x1v', x1.shape, '%.4e %.4e' % (x1[0], x1[-1]))
prof = {}
for f in sel:
    d = bc.read_binary(f)
    names = d['var_names']
    if f == sel[0]: print('vars', names)
    iv3 = [n for n in names if n.endswith('velz') or n.endswith('vel3')][0]
    iv1 = [n for n in names if n.endswith('velx') or n.endswith('vel1')][0]
    rho = np.array(d['mb_data'][[n for n in names if n.endswith('dens')][0]])
    v3 = np.array(d['mb_data'][iv3]); v1 = np.array(d['mb_data'][iv1])
    w = rho.sum(axis=(0, 1, 2))
    vp = np.sqrt((rho * v3**2).sum(axis=(0, 1, 2)) / w)
    vr = np.sqrt((rho * v1**2).sum(axis=(0, 1, 2)) / w)
    vpm = (rho * v3).sum(axis=(0, 1, 2)) / w
    rhom = rho.mean(axis=(0, 1, 2))
    prof[d['time']] = (vp, vr, vpm, rhom)
    print(f.split('/')[-1], 't=%.4e' % d['time'], 'nx1', vp.size)
ts = sorted(prof)
n = prof[ts[0]][0].size
print('i  rhomean(t0)  ' + '  '.join('vphi_rms@%.3g' % t for t in ts) + '  vr_rms@last  <vphi>@last')
for i in list(range(0, n - 60, 40)) + list(range(n - 60, n, 4)):
    print('%3d r=%.4e %.3e ' % (i, x1[i], prof[ts[0]][3][i]) + ' '.join('%.3e' % prof[t][0][i] for t in ts) + ' %.3e %.3e' % (prof[ts[-1]][1][i], prof[ts[-1]][2][i]))
last = prof[ts[-1]][0]; im = int(np.argmax(last)); print('max vphi_rms at last: i', im, '%.3e' % last[im])
