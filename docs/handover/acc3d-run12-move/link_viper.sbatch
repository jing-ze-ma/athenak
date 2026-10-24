#!/bin/bash -l
#SBATCH -J r12viper
#SBATCH -p apu
#SBATCH --constraint=apu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH -o /viper/ptmp2/jinma/run12_viper_1022/logs/link.%j.out
# run12 VIPER twin of the Raven acc3d prod chain (prod_chain.sbatch, job 31291344): same phases, keys, stop rules.
# Static phase (problem/stream=false time/tlim=$TS) in $PD/static from t=0, then stream (deck as is) in $PD/stream
# from the newest rst; links restart from the newest stream rst. Start guard: $W/CANCEL exists -> exit 0.
# Env overrides: PD, WT (s of simulation wall), TS, XK (extra keys, both phases), NR (ranks).
W=/viper/ptmp2/jinma/run12_viper_1022
mkdir -p $W/run12
[ -e $W/CANCEL ] && { echo "$(date) job $SLURM_JOB_ID: $W/CANCEL exists -> exit 0" | tee -a $W/run12/chain.txt; exit 0; }
module purge; module load gcc/16 rocm/7.1 openmpi_gpu/5.0
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1 UCX_RNDV_THRESH=inf UCX_TLS=^rocm_ipc OMP_NUM_THREADS=1
BN=$W/bin/athena_apu_ry3d_ff1f52bc_r71; DK=$W/acc3d_prod_D_B1v6_extrap2_viper.athinput
PD=${PD:-$W/run12}; WT=${WT:-85200}; TS=${TS:-0.05}; XK=${XK:-}; NR=${NR:-$SLURM_NTASKS}; T0=$SECONDS
md5sum $BN $DK; hostname; echo "job $SLURM_JOB_ID PD $PD WT $WT TS $TS XK '$XK' NR $NR nodes $SLURM_JOB_NODELIST"
stop_check() { for L in $PD/static/run.log $PD/stream/run.log; do [ -f $L ] || continue
    if grep -qE "NON-CONVERGED=[1-9]|FATAL|Picard solve DIVERGED" $L || tail -50 $L | grep -qwi "nan"; then
      echo "STOP: $L"; grep -nE "NON-CONVERGED=[1-9]|FATAL|DIVERGED" $L | tail -3; return 1; fi; done
  grep -q "STOPPED" $PD/chain.txt 2>/dev/null && { echo "STOP: earlier link stopped"; return 1; }; return 0; }
done_() { grep -q "Terminating on time limit" $1/run.log 2>/dev/null; }
newest() { ls -t $1/rst/ryper.*.rst 2>/dev/null | head -1; }
hms() { printf "%02d:%02d:%02d" $(($1/3600)) $(($1%3600/60)) $(($1%60)); }
mkdir -p $PD/static $PD/stream
stop_check || { echo "$(date) job $SLURM_JOB_ID STOPPED" >> $PD/chain.txt; exit 1; }
for ph in static stream; do
  D=$PD/$ph; done_ $D && continue
  LEFT=$((WT-(SECONDS-T0))); [ $LEFT -lt 600 ] && { echo "no wall left for $ph"; break; }
  R=$(newest $D); [ $ph = stream ] && [ -z "$R" ] && R=$(newest $PD/static)
  [ -z "$R" ] && IN="-i $DK" || IN="-r $R -i $DK"
  [ $ph = static ] && K="problem/stream=false time/tlim=$TS $XK" || K="$XK"
  cd $D || exit 1; echo "$(date) job $SLURM_JOB_ID $ph: $IN $K -t $(hms $LEFT)" | tee -a $PD/chain.txt
  srun -n $NR --cpu-bind=cores $BN $IN $K -d $D -t $(hms $LEFT) < /dev/null >> run.log 2>&1; rc=$?
  echo "$(date) job $SLURM_JOB_ID $ph rc $rc wall $((SECONDS-T0)) $(grep 'cycle=' run.log | tail -1)" | tee -a $PD/chain.txt
  # Raven pruned stream rst to the newest 3 here; this viper copy keeps all (prune by hand, 3.6 GB each)
  { [ $rc -ne 0 ] || ! stop_check; } && { echo "$(date) job $SLURM_JOB_ID STOPPED (rc $rc)" >> $PD/chain.txt; exit 1; }
  [ $ph = static ] && ! done_ $D && break
done
exit 0
