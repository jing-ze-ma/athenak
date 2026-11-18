import sys, types, time; sys.modules.setdefault('h5py', types.ModuleType('h5py'))
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/wt_accsp/vis/python')
import numpy as np, bin_convert as bc
t0=time.time(); h = bc.read_binary(sys.argv[1])
print('read %.1f s' % (time.time()-t0), h['time'], h['var_names'], h['n_mbs'], h['mb_data']['dens'][0].shape, h['mb_data']['dens'][0].dtype)
print('logical', h['mb_logical'].tolist())
for m in range(h['n_mbs']): print(m, h['mb_x3f'][m][0], h['mb_x3f'][m][-1], h['mb_x2f'][m][0], h['mb_x2f'][m][-1])
print('x2v deg', np.degrees(h['mb_x2v'][0]).round(2).tolist())
r=h['mb_x1v'][0]; print('r', r[0], r[-1], len(r), 'n r>8.4', (r>8.4).sum())
