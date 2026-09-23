"""compare vet_sc plane dumps (.plane.g*) of decomposed runs with a one-block reference.
usage: python3 cmp_planes.py REFDIR RUNDIR [RUNDIR ...]"""
import glob
import sys
import numpy as np


def load(d, call):
    fs = sorted(glob.glob('%s/vet.%06d.txt.plane*' % (d, call)))
    rows = np.vstack([np.loadtxt(f) for f in fs])
    ni = int(rows[:, 0].max()) + 1
    nj = int(rows[:, 1].max()) + 1
    a = np.zeros((nj, ni, rows.shape[1]))
    a[rows[:, 1].astype(int), rows[:, 0].astype(int)] = rows
    return a


ref = sys.argv[1]
calls = sorted(int(f.split('vet.')[1][:6]) for f in glob.glob(ref + '/vet.*.txt.plane*'))
for d in sys.argv[2:]:
    out = []
    for c in calls:
        a, b = load(ref, c), load(d, c)
        kj = np.abs(a[..., 3:9] - b[..., 3:9]).max(axis=-1)
        dd = np.abs(a[..., 16:22] - b[..., 16:22]).max(axis=-1)
        jr = np.abs(a[..., 2] / b[..., 2] - 1.0)
        out.append('call %d: K/J max %.2e rms %.2e | D max %.2e | J rel max %.2e'
                   % (c, kj.max(), np.sqrt((kj**2).mean()), dd.max(), jr.max()))
    print(d)
    for o in out:
        print('   ' + o)
