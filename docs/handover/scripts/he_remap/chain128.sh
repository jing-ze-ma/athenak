#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --constraint=gpu
#SBATCH --gres=gpu:a100:4
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
# He presn wedge 128^2 x 522 (he_mltpp_1002/remap128): continues pp2 at 70 ks from the remapped
# restart rst/hepresn.00033.rst (he_remap_rst.py: radial 470 -> 522 grid, 64^2 -> 128^2).
#   sbatch -t <time> -J he128 -o <RUNDIR>/job.%j.out chain128.sh <RUNDIR>
# Everything the run needs is in the run dir:
#   run.keys  command-line keys (one or more lines), incl. time/tlim; read on EVERY link
#   run.bin   the binary (one line)
# Guards: CANCEL -> exit 0; older RUNNING job of the same name -> exit 0; STOP -> exit 1;
# newest rst t >= tlim -> DONE.  Restart from the newest rst with <outputN>/last_time
# (rst_info.py).  The FIRST link (newest rst = the remapped hepresn.00033.rst) adds
# time/restart_refill_ghosts=true: the remap leaves the radial ghosts approximate.
# FATAL/NaN in run.log or hst -> STOP.  EXTRA env: appended keys (smoke: time/nlim=...).
D=${1:?run dir}
cd $D || exit 1
if [ -f $D/CANCEL ]; then echo "CANCEL present"; exit 0; fi
OLDER=$(squeue -h -u $USER -n "$SLURM_JOB_NAME" -t RUNNING -o %i | awk -v me=$SLURM_JOB_ID '$1<me' | head -1)
if [ -n "$OLDER" ]; then echo "guard: older RUNNING job $OLDER named $SLURM_JOB_NAME: exit"; exit 0; fi
if [ -f $D/STOP ]; then echo "STOP file present: $(cat $D/STOP)"; exit 1; fi
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0
export OMP_NUM_THREADS=1
BIN=$(head -1 $D/run.bin)
KEYS=$(grep -v '^#' $D/run.keys | tr '\n' ' ')
TLIM=$(echo $KEYS | tr ' ' '\n' | grep '^time/tlim=' | tail -1 | cut -d= -f2)
[ -n "$TLIM" ] || { echo "run.keys has no time/tlim"; exit 1; }
RS=$(ls $D/rst/*.rst 2>/dev/null | sort | tail -1)
[ -n "$RS" ] || { echo "no restart in $D/rst"; exit 1; }
INFO=$(python3 /raven/ptmp/jinma/he_raven_1002/rst_info.py $RS) || { echo "rst_info failed"; exit 1; }
T=$(echo "$INFO" | head -1 | awk '{print $1}')
if awk -v t=$T -v l=$TLIM 'BEGIN{exit !(t>=l)}'; then echo "newest rst $RS t=$T >= $TLIM" > $D/DONE; exit 0; fi
LT=$(echo "$INFO" | sed -n 2p)
FIRST=""
[ "$(basename $RS)" = "hepresn.00033.rst" ] && FIRST="time/restart_refill_ghosts=true"
JL=$(squeue -h -j $SLURM_JOB_ID -o %l)
S=$(echo $JL | awk -F'[-:]' '{n=NF; s=$n+60*$(n-1); if(n>=3)s+=3600*$(n-2); if(n>=4)s+=86400*$(n-3); print s}')
S=$((S-900)); [ $S -lt 300 ] && S=$((S+600))
WL=$(printf '%02d:%02d:%02d' $((S/3600)) $((S%3600/60)) $((S%60)))
md5sum $BIN $RS | head -2
echo "job $SLURM_JOB_ID ntasks $SLURM_NTASKS wall $WL restart $RS t=$T"
echo "keys: $LT $KEYS $FIRST $EXTRA"
[ -f run.log ] && mv run.log run.$(date +%m%d%H%M%S).log
srun -n $SLURM_NTASKS $BIN -r $RS -d $D -t $WL $LT $KEYS $FIRST $EXTRA > run.log 2>&1
RC=$?; echo "rc=$RC $(date +%T)"
HST=$(ls $D/*.hst 2>/dev/null | head -1)
if grep -qiwE 'fatal|nan|-nan' run.log || { [ -n "$HST" ] && tail -1 $HST | grep -qiw nan; }; then
  echo "job $SLURM_JOB_ID: FATAL/NaN ($(date +%F_%T))" > $D/STOP; cat $D/STOP; exit 1
fi
exit $RC
