# run "weight<TAB>command" lines with the sum of running weights <= cap (login node)
import sys, subprocess, time
cap = int(sys.argv[2]); jobs = []
for ln in open(sys.argv[1]):
    ln = ln.strip()
    if ln and not ln.startswith('#'):
        w, c = ln.split(None, 1); jobs.append((int(w), c))
run = []
while jobs or run:
    run = [(w, p) for (w, p) in run if p.poll() is None]
    used = sum(w for w, _ in run)
    for i, (w, c) in enumerate(jobs):
        if used + w <= cap:
            run.append((w, subprocess.Popen(c, shell=True))); used += w; jobs.pop(i); break
    else:
        time.sleep(2)
print('SCHED DONE', flush=True)
