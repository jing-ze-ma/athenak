#!/bin/bash
# idortd6int-1113: the 7 serial-bound 1-D rows (lists/list_omp.txt) with 10 OpenMP threads, arms p (EXE_OMP_P) and f (EXE_OMP_F),
# as Raven jobs/time_omp.sbatch (there: correctness + timing; here correctness only: CPU timings are not results).
# The same 7 decks also run serially in list_mainrk (runs/); cmpout runs_omp vs runs of the same arm is an extra OpenMP check.
# usage: omp_o.sh   env: C W EXE_OMP_P EXE_OMP_F [NLIM]   run dirs: $W/runs_omp/<row>/{p,f}
for r in $(awk '{print $1}' "$C/lists/list_omp.txt"); do for x in p f; do
  case $x in p) dk=$C/decks_d6/$r; B=$EXE_OMP_P;; f) dk=$C/decks_fin/$r; B=$EXE_OMP_F;; esac
  ( d="${W:?}/runs_omp/${r:?}/$x"; mkdir -p "$d"; cd "$d" || exit; [ -f rc.txt ] && [ "$(cat rc.txt)" = 0 ] && exit
    cp "$dk/in.athinput" "$dk/timekeys.txt" .; md5sum "$B" > bin.txt; t0=$(date +%s)
    OMP_NUM_THREADS=10 OMP_PROC_BIND=false "$B" -i in.athinput -d "$d" ${NLIM:+time/nlim=$NLIM} > run.log 2>&1; echo $? > rc.txt
    cp run.log all.log; echo "$(( $(date +%s)-t0 )) s nt=10" > time.txt; rm -rf "${d:?}/rst"
    echo "omp $r $x rc=$(cat rc.txt) fatal=$(grep -c FATAL all.log) $(grep -h -o 'NON-CONVERGED=[^ ]*' all.log | tail -1)" ) &
done; done; wait; echo OMP_DONE
