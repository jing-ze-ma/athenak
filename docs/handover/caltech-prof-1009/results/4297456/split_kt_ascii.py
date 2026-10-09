#!/usr/bin/env python3
"""Recover simple-kernel-timer ASCII tables (kokkos-tools master default: printed to stdout,
no .dat) from an out.log where the 2 ranks' tables are interleaved in whole-line chunks.
Pairs each '(Type) t calls ...' line with the name line directly above it, else with the oldest
dangling name. Writes kt_mean.txt (times = mean over ranks, kp_reader-like) for prof_group.py.
usage: split_kt_ascii.py <out.log> <nranks> <outfile>"""
import re, sys
from collections import defaultdict
log, nr, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
num = re.compile(r'^\s*\((ParFor|ParRed|ParScan|Region)\)\s+([0-9.eE+-]+)\s+(\d+)\s')
dang, prev = [], None
ent = defaultdict(list)
tot, ker = [], []
for line in open(log, errors='replace'):
    s = line.rstrip('\n')
    m = re.search(r'Total Execution Time.*?([0-9.]+)\s*seconds', s)
    if m: tot.append(float(m.group(1)))
    m = re.search(r'Total Time in Kokkos kernels.*?([0-9.]+)\s*seconds', s)
    if m: ker.append(float(m.group(1)))
    if s.startswith('- ') and len(s) > 2:
        dang.append(s[2:].strip()); prev = 'name'; continue
    m = num.match(s)
    if m:
        name = dang.pop() if prev == 'name' else (dang.pop(0) if dang else None)
        if name is None: sys.exit('numbers line without a name: ' + s)
        if m.group(1) != 'Region':
            ent[name].append((float(m.group(2)), int(m.group(3))))
        prev = 'num'; continue
    if s.strip(): prev = None
bad = {k: len(v) for k, v in ent.items() if len(v) != nr}
ksum = sum(t for v in ent.values() for t, c in v)
sys.stderr.write('%s: %d kernels, %d with != %d entries %s, ranks total=%s kernels=%s, '
                 'sum of entries %.4f vs sum of kernel totals %.4f, dangling %d\n' % (
                     log, len(ent), len(bad), nr, list(bad.items())[:5], tot, ker, ksum, sum(ker), len(dang)))
with open(out, 'w') as f:
    f.write('Total Execution Time (incl. Kokkos + non-Kokkos): %.5f seconds\n' % (sum(tot) / len(tot)))
    f.write('Total Time in Kokkos kernels: %.5f seconds\n' % (sum(ker) / len(ker)))
    for k, v in ent.items():
        f.write('- %s\n (ParFor) %.6f %d\n' % (k, sum(t for t, c in v) / nr, round(sum(c for t, c in v) / nr)))
