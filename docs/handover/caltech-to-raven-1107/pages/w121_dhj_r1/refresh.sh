#!/bin/bash
# dhj_r1 page refresh (pages.sh convention).
#   refresh.sh         quick index.html rebuild on login from cached reductions (no dump reads, ~10 s), then a DETACHED
#                      CPU-only catch-up step in a held allocation (new dumps only; figures; index.html; writes READY).
#   refresh.sh status  one line, instant.
P=/resnick/home/jingze/ATHENAK/pages/dhj_r1
LOCK=$P/logs/catchup.lock
cd $P || exit 1
mkdir -p logs
if [ "$1" = status ]; then
  if ! flock -n $LOCK true 2>/dev/null; then echo "dhj_r1: catch-up running ($(cat logs/catchup.where 2>/dev/null)); last: $(head -1 READY 2>/dev/null || echo none)"
  else echo "dhj_r1: $(head -1 READY 2>/dev/null || echo 'not ready')"; fi
  exit 0
fi
JOB=; NODE=
for c in 5274196:hpc-sm-01-11 5719508:hpc-sm-01-05; do
  j=${c%%:*}; n=${c##*:}
  [ "$(squeue -h -j $j -o %T 2>/dev/null)" = RUNNING ] && { JOB=$j; NODE=$n; break; }
done
# pagefast-1102: the quick html rebuild runs as a 2-core CPU step in the held allocation (login only if none is running)
QH='module load python/3.11.6-gcc-11.3.1-rphh4kv py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1; cd '$P'/tools && timeout 60 python3 ana.py html'
if [ -n "$JOB" ]; then timeout 120 srun --overlap --jobid=$JOB -w $NODE -N1 -n1 -c 2 --mem=16G --gres=none -J dhjhtml nice -n 10 bash -c "$QH" > logs/html.out 2>&1 < /dev/null
else bash -c "$QH" > logs/html.out 2>&1; fi && echo "dhj_r1: index.html rebuilt from cache" || echo "dhj_r1: quick html FAILED (logs/html.out)"
if ! flock -n $LOCK true 2>/dev/null; then echo "dhj_r1: catch-up already running ($(cat logs/catchup.where 2>/dev/null)), not relaunched"; exit 0; fi
[ -z "$JOB" ] && { echo "dhj_r1: no held allocation RUNNING, catch-up skipped"; exit 0; }
cat > logs/step.sh <<EOS
#!/bin/bash
# runs on login under flock; srun = CPU-only step in the held allocation
echo "$JOB $NODE since \$(date +%m-%d\ %H:%M)" > $P/logs/catchup.where
srun --overlap --jobid=$JOB -w $NODE -N1 -n1 -c 4 --mem=40G --gres=none -J dhjpage \\
  bash -c 'module load python/3.11.6-gcc-11.3.1-rphh4kv py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1; cd $P/tools && NPROC=4 python3 ana.py full' \\
  > $P/logs/ana.out 2>&1 < /dev/null
L=\$(grep '^dhj_r1 READY' $P/logs/ana.out | tail -1 | sed 's/^dhj_r1 READY //')
if [ -n "\$L" ]; then
  { echo "READY \$(date +%m-%d\ %H:%M) \$L"; cd $P; echo index.html; ls img/*.png; } > $P/READY.tmp && mv $P/READY.tmp $P/READY
else
  echo "catch-up FAILED \$(date +%m-%d\ %H:%M), see logs/ana.out" >> $P/logs/failures.txt
fi
EOS
chmod +x logs/step.sh
setsid nohup flock -n $LOCK $P/logs/step.sh > /dev/null 2>&1 < /dev/null &
echo "dhj_r1: catch-up launched in $JOB ($NODE); 'refresh.sh status' / READY when done"
# fehook-1102 (USER): fe_hires Fe I update, detached with its own lock (non-blocking; the page builds without it); launched
# here, not from step.sh, so it does not inherit the catch-up lock
if flock -n $P/logs/fe.lock true; then setsid nohup flock -n $P/logs/fe.lock $P/fe_step.sh $JOB $NODE > /dev/null 2>&1 < /dev/null &
  echo "dhj_r1: fe_hires update launched (logs/fe.out)"; else echo "dhj_r1: fe_hires update still running (logs/fe.out)"; fi
