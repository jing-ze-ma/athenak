import sys, glob, os
sys.path.insert(0,os.environ['ATHENA_VIS'])
import bin_convert as bc, numpy as np
fd, hd = sys.argv[1], sys.argv[2]
res = []
for f in sorted(glob.glob(fd + '/bin/*.00001.bin')):
    a = bc.read_binary(f); b = bc.read_binary(f.replace(fd, hd))
    d = max(np.abs(np.array(a['mb_data'][k]) - np.array(b['mb_data'][k])).max() for k in a['var_names'])
    res.append('%s %s' % (os.path.basename(f).split('.')[1], 'same' if d == 0 else 'DIFF %.1e' % d))
print(' '.join(res))
