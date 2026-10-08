#!/bin/bash -l
# he_giant_1006 scout128 (10-07): copy of /viper/ptmp2/jinma/bsg_viper_1006/dual_viper.sh (made by mk_dual.py)
# with ONE change: FRESH START from IN (run.cfg) when $D/rst has no rst AND no run log exists (first link only);
# every later link restarts.  Every link:  sbatch -J <name> --nodes=<N> --time=<T> -o $D/dual.%j.out -e $D/dual.%j.err dual_viper.sh $D
# $D/run.cfg (shell KEY=VALUE, read at job start, editable for pending links): BIN, BIN_MD5, TLIM, XKEYS, PK
#   XKEYS = the Raven command-line keys (identical);  PK = problem/ path overrides (viper copies of the Raven
#   IC/opacity tables, md5 verified equal); the restart carries the rest of the input.
# (a) GUARD: a RUNNING job of this user and name that started earlier (same second: lower id) -> exit 0.
# (b) CANCEL -> exit 0; DONE or newest rst time >= TLIM -> exit 0; STOP -> exit 1; md5 != BIN_MD5 -> STOP.
# (c) restart from the newest rst with outputN/last_time reset (rst_info.py); NO fresh start on viper (no rst -> exit 1);
#     before a restart linkcheck.py checks the previous run.log (NaN/FATAL, newer rst, dt collapse).
# (d) code wall limit = job limit - elapsed - 20 min; NaN/FATAL in run.log afterwards -> STOP.
# Args: $1 = run dir (with its run.cfg and rst/), $2 = cycles past the start point (smoke only).
#SBATCH -J hegiant
#SBATCH -p apu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --constraint=apu
#SBATCH --cpus-per-task=24
#SBATCH --time=04:00:00
export SBATCH_EXPORT=NONE
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
INFO=/viper/ptmp2/jinma/bsg_viper_1006/rst_info.py
LCHK=/viper/ptmp2/jinma/bsg_viper_1006/linkcheck.py
D=$1
[ -n "$D" ] && [ -f $D/run.cfg ] || { echo "no run dir / run.cfg ($D): exit 1"; exit 1; }
BIN=""; BIN_MD5=""; TLIM=""; XKEYS=""; PK=""; IN=""
source $D/run.cfg
[ -x "$BIN" ] && [ -n "$BIN_MD5" ] && [ -n "$TLIM" ] || { echo "run.cfg incomplete (BIN=$BIN TLIM=$TLIM): exit 1"; exit 1; }
echo "run.cfg: BIN=$BIN TLIM=$TLIM"
echo "XKEYS: $XKEYS"
echo "PK: $PK"
NCYC=${2:-}
cd $D || exit 1
echo "job $SLURM_JOB_ID ($SLURM_JOB_NAME) $SLURM_NNODES nodes $SLURM_JOB_NODELIST, $SLURM_NTASKS ranks, start $(date +%F_%T)"
[ -f CANCEL ] && { echo "CANCEL exists: exit 0"; exit 0; }
sq() { local i out; for i in 1 2 3; do out=$(squeue -h "$@") && { echo "$out"; return 0; }; sleep 10; done; return 1; }
ME=$(sq -j $SLURM_JOB_ID -o %S) || { echo "GUARD: squeue failed: exit 0"; exit 0; }
RUN=$(sq -u $USER -t R -n $SLURM_JOB_NAME -o "%i,%S") || { echo "GUARD: squeue failed: exit 0"; exit 0; }
[[ "$ME" == 20* ]] || { echo "GUARD: no start time for myself ($ME): exit 0"; exit 0; }
echo "GUARD: me $SLURM_JOB_ID start $ME; running of this name: $(echo $RUN)"
for l in $RUN; do
  id=${l%%,*}; st=${l#*,}
  [ "$id" = "$SLURM_JOB_ID" ] && continue
  if [[ "$st" < "$ME" ]] || { [ "$st" = "$ME" ] && [ "$id" -lt "$SLURM_JOB_ID" ]; }; then
    echo "GUARD: older running job $id (start $st): exit 0"; exit 0
  fi
done
echo "GUARD: I run"
[ -f CANCEL ] && { echo "CANCEL exists: exit 0"; exit 0; }
[ -f DONE ] && { echo "DONE exists: run finished: exit 0"; exit 0; }
[ -f STOP ] && { echo "STOP exists ($(cat STOP)): exit 1"; exit 1; }
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
M=$(md5sum $BIN | cut -d' ' -f1); echo "$M  $BIN"
[ "$M" != "$BIN_MD5" ] && { echo "md5 $M != BIN_MD5 $BIN_MD5: STOP"; echo "md5 mismatch (job $SLURM_JOB_ID)" > STOP; exit 1; }
LIM=$(squeue -h -j $SLURM_JOB_ID -o %l)      # [D-]H:MM:SS
dd=0; [[ $LIM == *-* ]] && { dd=${LIM%%-*}; LIM=${LIM#*-}; }
IFS=: read -r a b c <<< "$LIM"; [ -z "$c" ] && { c=$b; b=$a; a=0; }
LIMS=$(( dd*86400 + 10#$a*3600 + 10#$b*60 + 10#$c ))
MARGIN=1200; [ -n "$NCYC" ] && MARGIN=180
LEFT=$(( LIMS - SECONDS - MARGIN ))
[ $LEFT -lt 600 ] && { echo "less than 10 min of wall left: exit 0"; exit 0; }
TL=$(printf '%02d:%02d:%02d' $((LEFT/3600)) $(((LEFT%3600)/60)) $((LEFT%60)))
R=$(ls $D/rst/*.rst 2>/dev/null | sort | tail -1)
if [ -z "$R" ]; then
  { [ -n "$IN" ] && [ -f "$IN" ] && [ ! -f run.log ] && ! ls run.*.log >/dev/null 2>&1; } || \
    { echo "no rst in $D/rst and no fresh start allowed (IN=$IN, earlier logs?): exit 1"; exit 1; }
  ARGS="-i $IN time/tlim=$TLIM $XKEYS $PK"
  [ -n "$NCYC" ] && ARGS="$ARGS time/nlim=$NCYC time/ndiag=1"
  echo fresh > .link_start
  echo "FRESH start from $IN ($(md5sum $IN | cut -d' ' -f1))"
  echo "srun -n $SLURM_NTASKS $BIN $ARGS -d $D -t $TL   $(date +%T)"
  srun -n $SLURM_NTASKS --cpu-bind=cores $BIN $ARGS -d $D -t $TL > run.log 2>&1
  rc=$?
  echo "rc=$rc $(date +%T), newest rst $(ls $D/rst/*.rst 2>/dev/null | sort | tail -1)"
  if grep -q "FATAL" run.log || tail -n 60 run.log | grep -qiw "nan"; then
    echo "NaN/FATAL in run.log: STOP"; echo "NaN/FATAL in run.log of job $SLURM_JOB_ID (rc $rc)" > STOP; exit 1
  fi
  [ $rc -ne 0 ] && { echo "srun rc=$rc: STOP"; echo "srun rc=$rc in job $SLURM_JOB_ID" > STOP; exit 1; }
  exit 0
fi
INF=$(/usr/bin/python3 $INFO $R) || { echo "rst_info failed on $R: exit 1"; exit 1; }
read -r T NC TL0 <<< "$(echo "$INF" | head -1)"
awk -v t=$T -v l=$TLIM 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "rst t=$T >= tlim: DONE"; touch DONE; exit 0; }
if [ -f run.log ]; then
  CK=$(/usr/bin/python3 $LCHK $D 2>&1) || { echo "linkcheck refused: $CK"; echo "linkcheck ($SLURM_JOB_ID): $CK" > STOP; exit 1; }
  echo "$CK"
  mv run.log run.$(date +%m%d_%H%M%S).log
fi
ARGS="-r $R $(echo "$INF" | tail -1) time/tlim=$TLIM $XKEYS $PK"
basename $R > .link_start
echo "restart from $R (t = $T, cycle $NC)"
[ -n "$NCYC" ] && ARGS="$ARGS time/nlim=$(( NC + NCYC )) time/ndiag=1"
echo "srun -n $SLURM_NTASKS $BIN $ARGS -d $D -t $TL   $(date +%T)"
( while true; do echo "$(date +%T) $(rocm-smi --showmemuse --showuse --csv 2>/dev/null | tr '\n' ' ')" >> gpumem.$SLURM_JOB_ID.log; sleep 60; done ) &
MON=$!
srun -n $SLURM_NTASKS --cpu-bind=cores $BIN $ARGS -d $D -t $TL > run.log 2>&1
rc=$?
kill $MON 2>/dev/null
R2=$(ls $D/rst/*.rst 2>/dev/null | sort | tail -1)
echo "rc=$rc $(date +%T), newest rst $R2"
if grep -q "FATAL" run.log || tail -n 60 run.log | grep -qiw "nan"; then
  echo "NaN/FATAL in run.log: STOP"; echo "NaN/FATAL in run.log of job $SLURM_JOB_ID (rc $rc)" > STOP; exit 1
fi
[ $rc -ne 0 ] && { echo "srun rc=$rc: STOP"; echo "srun rc=$rc in job $SLURM_JOB_ID" > STOP; exit 1; }
if [ -n "$R2" ]; then
  T2=$(/usr/bin/python3 $INFO $R2 | head -1 | cut -d' ' -f1)
  awk -v t=$T2 -v l=$TLIM 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "rst t=$T2 >= tlim: DONE"; touch DONE; }
fi
exit 0
