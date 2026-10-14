#!/bin/bash
# PU pulse-1014 (DeltaAI): one link of 4 independent single-GPU arms (REDUCED 3-D E1, nx1 536, lateral 32x32 = 1 MeshBlock per arm).
# Env: BIN (fin4-dtretry 5b7e1bc5 build), DECK (resolved deck, md5 16b0bba7 with the Raven paths), S (scratch root), WALL (default 01:52:00),
#      LAUNCH (default "srun --exact -n 1 -c 16 --gpus-per-task=1"). Link n restarts every arm from its newest rst.
export MPICH_GPU_SUPPORT_ENABLED=1 OMP_NUM_THREADS=1
LAUNCH=${LAUNCH:-"srun --exact -n 1 -c 16 --gpus-per-task=1"}
TL=1.25e7
RED="mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32"
OUTK="time/ndiag=200 time/dt_retry=true output3/dt=1e6 output5/dt=1e6 time/tlim=$TL"
ARMS=("A536|" "S0_536|problem/he_seed=0" "S3_536|problem/he_seed=1.0e-3" "F536|problem/mlt_ramp_start=1.0e6 problem/mlt_ramp_time=2.0e6")
for a in "${ARMS[@]}"; do
  IFS='|' read -r NM KY <<< "$a"; D=$S/$NM; mkdir -p $D
  [ -f $D/STOP ] && continue
  R=$(ls $D/rst/*.rst 2>/dev/null | sort | tail -1); L=$(ls $D/run*.log 2>/dev/null | wc -l)
  if [ -z "$R" ]; then ARGS="-i $DECK $RED $OUTK $KY"; else ARGS="-r $R $KY time/tlim=$TL"; fi
  echo "link $L job ${SLURM_JOB_ID:-} $(md5sum $BIN | cut -c1-8) $ARGS" >> $D/cmd.txt
  ( cd $D; $LAUNCH $BIN $ARGS -d $D -t ${WALL:-01:52:00} < /dev/null > run$L.log 2>&1; rc=$?; echo "rc=$rc" > rc$L.txt
    grep -q "FATAL" run$L.log && echo "FATAL link $L" > STOP ) &
done
wait; echo LINK_DONE
