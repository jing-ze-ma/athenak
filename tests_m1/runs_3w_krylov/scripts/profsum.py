"""launches and host syncs per cycle (rank 0) = (nlim20 - nlim10)/10.
usage: python3 profsum.py runs/<tag>"""
import csv
import glob
import os
import sys
from collections import Counter

SYNC = ('hipStreamSynchronize', 'hipDeviceSynchronize', 'hipEventSynchronize',
        'hipMemcpy', 'hipMemcpyDtoH', 'hipMemcpyHtoD', 'hipStreamWaitEvent')


def count(d):
    kf = glob.glob(d + '/prof/**/rank_0*kernel_trace.csv', recursive=True)
    hf = glob.glob(d + '/prof/**/rank_0*hip_api_trace.csv', recursive=True)
    if not kf or not hf:
        return None
    nk = sum(1 for _ in open(kf[0])) - 1
    c = Counter()
    kt = Counter()
    with open(kf[0]) as f:
        for r in csv.DictReader(f):
            kt[r.get('Kernel_Name', '')[:60]] += 1
    with open(hf[0]) as f:
        for r in csv.DictReader(f):
            c[r.get('Function', r.get('Operation', ''))] += 1
    return nk, c, kt


D = sys.argv[1]
arms = sorted({os.path.basename(p).rsplit('_n', 1)[0] for p in glob.glob(D + '/*_n10')})
print('%-10s %10s %10s  %s' % ('arm', 'launch/cyc', 'sync/cyc', 'sync calls per cycle'))
for a in arms:
    x, y = count(D + '/' + a + '_n10'), count(D + '/' + a + '_n20')
    if not x or not y:
        print(a, 'missing')
        continue
    nk = (y[0] - x[0]) / 10.0
    per = {k: (y[1][k] - x[1][k]) / 10.0 for k in set(y[1]) | set(x[1])}
    sy = {k: v for k, v in per.items() if k in SYNC and v}
    print('%-10s %10.1f %10.1f  %s' % (a, nk, sum(sy.values()),
                                       ' '.join('%s=%.1f' % kv for kv in sorted(sy.items()))))
