#!/usr/bin/env python3
"""sched.py <list> <cores> <one.sh>: run 'one.sh <dir> <nt>' for every '<dir> <nt> [est_s]' line, longest est first (strict order),
each job taking NARM*nt cores (arms concurrently; NARM = number of words in env ARMS, default "ref r"), never more than <cores> in use. Prints each job's output line."""
import sys, subprocess, time
L, C, ONE = sys.argv[1], int(sys.argv[2]), sys.argv[3]
import os
NA = len(os.environ.get("ARMS", "ref r").split())
jobs = []
for ln in open(L):
    w = ln.split()
    if not w or w[0].startswith('#'): continue
    jobs.append((w[0], int(w[1]), float(w[2]) if len(w) > 2 else 1e9))
jobs.sort(key=lambda j: -j[2])
run, used = [], 0
while jobs or run:
    for p, c in run[:]:
        if p.poll() is not None:
            sys.stdout.write(p.stdout.read().decode().strip().split('\n')[-1] + '\n'); sys.stdout.flush()
            run.remove((p, c)); used -= c
    started = False
    if jobs:   # strict order (longest first): no backfill past a job that does not fit
        j = jobs[0]; c = min(NA * j[1], C)
        if used + c <= C:
            run.append((subprocess.Popen(['bash', ONE, j[0], str(j[1])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT), c))
            used += c; jobs.pop(0); started = True
    if not started: time.sleep(2)
