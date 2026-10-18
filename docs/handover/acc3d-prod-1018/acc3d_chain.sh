#!/bin/bash
# ACC3D production restart chain (one link = one batch job; run it INSIDE the job allocation). Raven original:
# /raven/ptmp/jinma/acc3d_1017/jobs/prod_chain.sbatch (same logic). Required env:
#   BN     = athena binary (accretor-3d-1017 @ 186306cb, PROBLEM=ry_per_accretor, nofma)
#   DK     = the substituted deck (acc3d_prod_mk.athinput)
#   PD     = run dir (static/ and stream/ subdirs are made here)
#   WT     = seconds of simulation wall for this link (job walltime minus ~20 min for the final 2.9 GB rst + margin)
#   LAUNCH = the MPI launch prefix incl. your CUDA_VISIBLE_DEVICES wrapper, e.g. "srun -n 4 --cpu-bind=cores ./gpuwrap.sh"
# Optional: TS (static-phase tlim, default 0.05), XK (extra overrides for both phases).
# Phases: static (problem/stream=false, tlim TS) in $PD/static, then stream (deck as is: stream on, tlim 1e9) in $PD/stream
# from the newest rst; runs to the wall. Keeps the 3 newest stream rst (2.9 GB each). Exit 1 = STOPPED (dependents must not run).
TS=${TS:-0.05}; XK=${XK:-}; T0=$SECONDS
: ${BN:?} ${DK:?} ${PD:?} ${WT:?} ${LAUNCH:?}
md5sum $BN $DK; hostname; echo "PD $PD WT $WT TS $TS XK '$XK'"
stop_check() { for L in $PD/static/run.log $PD/stream/run.log; do [ -f $L ] || continue
    if grep -qE "NON-CONVERGED=[1-9]|FATAL|Picard solve DIVERGED" $L || tail -50 $L | grep -qwi "nan"; then
      echo "STOP: $L"; grep -nE "NON-CONVERGED=[1-9]|FATAL|DIVERGED" $L | tail -3; return 1; fi; done
  grep -q "STOPPED" $PD/chain.txt 2>/dev/null && { echo "STOP: earlier link stopped"; return 1; }; return 0; }
done_() { grep -q "Terminating on time limit" $1/run.log 2>/dev/null; }
newest() { ls -t $1/rst/ryper.*.rst 2>/dev/null | head -1; }
hms() { printf "%02d:%02d:%02d" $(($1/3600)) $(($1%3600/60)) $(($1%60)); }
mkdir -p $PD/static $PD/stream
stop_check || { echo "$(date) STOPPED" >> $PD/chain.txt; exit 1; }
for ph in static stream; do
  D=$PD/$ph; done_ $D && continue
  LEFT=$((WT-(SECONDS-T0))); [ $LEFT -lt 600 ] && { echo "no wall left for $ph"; break; }
  R=$(newest $D); [ $ph = stream ] && [ -z "$R" ] && R=$(newest $PD/static)
  [ -z "$R" ] && IN="-i $DK" || IN="-r $R -i $DK"
  [ $ph = static ] && K="problem/stream=false time/tlim=$TS $XK" || K="$XK"
  cd $D; echo "$(date) $ph: $IN $K -t $(hms $LEFT)" | tee -a $PD/chain.txt
  $LAUNCH $BN $IN $K -d $D -t $(hms $LEFT) < /dev/null >> run.log 2>&1; rc=$?
  echo "$(date) $ph rc $rc wall $((SECONDS-T0)) $(grep 'cycle=' run.log | tail -1)" | tee -a $PD/chain.txt
  [ $ph = stream ] && ls -t $D/rst/ryper.*.rst 2>/dev/null | tail -n +4 | xargs -r rm -f
  { [ $rc -ne 0 ] || ! stop_check; } && { echo "$(date) STOPPED (rc $rc)" >> $PD/chain.txt; exit 1; }
  [ $ph = static ] && ! done_ $D && break
done
exit 0
