#!/bin/bash
# 10-14: build bsgli_1014/pagerun = prod (2 GPU, stopped) up to the wide4 seed time + wide4/run (4 GPU) after.
# Symlinks for dumps (wide wins on overlapping numbers); hst rows and run.log concatenated. Read-only on both runs.
B=/resnick/groups/carnegie_poc/jingze/bsgli_1014; P=$B/prod; Q=$B/wide4/run; M=$B/pagerun
mkdir -p $M/bin; rm -f $M/bin/*
for f in $P/bin/*.bin; do ln -s $f $M/bin/; done
for f in $Q/bin/*.bin; do ln -sf $f $M/bin/; done
for f in bsg3d.x1grid.txt bsg3d.log; do ln -sf $Q/$f $M/$f; done
python3 - $P $Q $M <<'PY'
import sys
P, Q, M = sys.argv[1:]
for h in ('bsg3d.hydro.hst', 'bsg3d.user.hst'):
    q = open(Q + '/' + h).read().splitlines()
    qd = [l for l in q if l.strip() and not l.startswith('#')]
    t0 = float(qd[0].split()[0])
    p = open(P + '/' + h).read().splitlines()
    out = [l for l in p if l.startswith('#') or (l.strip() and float(l.split()[0]) < t0)] + qd
    open(M + '/' + h, 'w').write('\n'.join(out) + '\n')
PY
cat $P/run.log $Q/run.log > $M/run.log
