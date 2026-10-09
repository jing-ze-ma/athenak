#!/usr/bin/env python3
"""Replay a Kokkos memory-events file (device space only) and print the live allocations
at the peak, ranked by size, plus the static set alive at the first 'cycle' marker-free point:
the end of setup is approximated by the time of the largest gap... we simply report
(a) peak, (b) final-before-teardown (max live time ordering), grouped by label.
usage: ana_mem.py <file.mem_events> [ncells_per_block] [nblocks]"""
import sys
import collections

fn = sys.argv[1]
ncell = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
nmb = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
live = {}
cur = 0
peak = 0
peak_set = None
peak_t = 0.0
lastsnap = None
for line in open(fn):
    if line.startswith('#') or 'Region' in line:
        continue
    p = line.split(None, 5)
    if len(p) < 5:
        continue
    t, ptr, size, space, op = p[:5]
    name = p[5].strip() if len(p) > 5 else ''
    if 'Host' in space and 'Pinned' not in space:
        continue
    size = int(size)
    if op == 'Allocate':
        live[ptr] = (size, name, space)
        cur += size
        if cur > peak:
            peak = cur
            peak_t = float(t)
            peak_set = dict(live)
    else:
        if ptr in live:
            cur -= live[ptr][0]
            del live[ptr]


def table(s, title):
    g = collections.defaultdict(lambda: [0, 0])
    for size, name, space in s.values():
        g[(name, space)][0] += size
        g[(name, space)][1] += 1
    tot = sum(v[0] for v in g.values())
    print('## %s: total %.3f GB (%d allocations)' % (title, tot / 1e9, len(s)))
    print('| label | space | n | GB per rank | GB per block | doubles/cell |')
    print('|---|---|---|---|---|---|')
    for (name, space), (sz, n) in sorted(g.items(), key=lambda x: -x[1][0]):
        if sz < 2e7:
            continue
        dpc = sz / 8.0 / (ncell * nmb) if ncell > 0 else 0
        print('| %s | %s | %d | %.3f | %.3f | %.1f |' % (name, space, n, sz / 1e9,
                                                      sz / 1e9 / nmb, dpc))
    print()


table(peak_set, 'PEAK at t=%.2f s' % peak_t)
