#!/bin/bash
# AG Car production link for DeltaAI (1 node x 4 GH200, 4 ranks, 1 MeshBlock 480x64x64 per GPU).
# Same logic as viper rc_1009/jobs/rc_raven.sh + bsg truerepro2 link (chain, rst_info restart, guards).
# usage: sbatch -A bivj-dtai-gh -p <ghx4|ghx4-interactive> --time=<max> -J agcB_dai \
#          -o $D/link.%j.out agcar_dai_link.sh $D [NCYC]
#   $D/run.cfg: BIN BIN_MD5 IN TLIM [XKEYS].  No rst in $D/rst -> FRESH start from IN (t = 0); else restart from
#   the newest rst with outputN/last_time reset ($D/rst_info.py, or RST_INFO=<path>).
#   NCYC (smoke only): stop NCYC cycles after the start point (time/nlim), ndiag=1.
#   $D/CANCEL or $D/DONE -> exit 0; $D/STOP -> exit 1; md5 mismatch, NaN/FATAL, rc != 0 -> write STOP, exit 1.
#   GUARD: an older RUNNING job with the same name -> exit 0.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1   # job-env rule (ROCm-only, no-op on CUDA)
D=$1
INFO=${RST_INFO:-$D/rst_info.py}   # sbatch runs a spool copy of this script: keep rst_info.py IN the run dir
NCYC=${2:-}
[ -n "$D" ] && [ -f "$D/run.cfg" ] || { echo "no run dir / run.cfg ($D): exit 1"; exit 1; }
BIN=""; BIN_MD5=""; IN=""; TLIM=""; XKEYS=""
source "$D/run.cfg"
cd "$D" || exit 1
echo "job $SLURM_JOB_ID ($SLURM_JOB_NAME) $SLURM_JOB_PARTITION $SLURM_JOB_NODELIST, $SLURM_NTASKS ranks, start $(date +%F_%T)"
[ -f CANCEL ] && { echo "CANCEL: exit 0"; exit 0; }
ME=$(squeue -h -j "$SLURM_JOB_ID" -o %S)
for l in $(squeue -h -u "$USER" -t R -n "$SLURM_JOB_NAME" -o "%i,%S"); do
  id=${l%%,*}; st=${l#*,}
  [ "$id" = "$SLURM_JOB_ID" ] && continue
  if [[ "$st" < "$ME" ]] || { [ "$st" = "$ME" ] && [ "$id" -lt "$SLURM_JOB_ID" ]; }; then
    echo "GUARD: older running job $id: exit 0"; exit 0
  fi
done
[ -f DONE ] && { echo "DONE: exit 0"; exit 0; }
[ -f STOP ] && { echo "STOP ($(cat STOP)): exit 1"; exit 1; }
[ -f "$INFO" ] || { echo "no rst_info.py ($INFO): exit 1"; exit 1; }
[ -x "$BIN" ] && [ -f "$IN" ] && [ -n "$TLIM" ] || { echo "run.cfg incomplete: exit 1"; exit 1; }
M=$(md5sum "$BIN" | cut -d' ' -f1)
[ "$M" != "$BIN_MD5" ] && { echo "md5 $M != $BIN_MD5: STOP"; echo "md5 mismatch" > STOP; exit 1; }
cp -p run.cfg "run.cfg.$SLURM_JOB_ID"
md5sum "$BIN" "$IN" $(grep -E '^ *(he_ic_file|he_opac_table|he_planck_table) ' "$IN" | awk '{print $3}')
LIM=$(squeue -h -j "$SLURM_JOB_ID" -o %l)
dd=0; [[ $LIM == *-* ]] && { dd=${LIM%%-*}; LIM=${LIM#*-}; }
IFS=: read -r a b c <<< "$LIM"; [ -z "$c" ] && { c=$b; b=$a; a=0; }
LIMS=$(( dd*86400 + 10#$a*3600 + 10#$b*60 + 10#$c ))
MARGIN=900; [ $LIMS -le 3600 ] && MARGIN=240
LEFT=$(( LIMS - SECONDS - MARGIN ))
TL=$(printf '%02d:%02d:%02d' $((LEFT/3600)) $(((LEFT%3600)/60)) $((LEFT%60)))
mkdir -p rst
R=$(ls rst/*.rst 2>/dev/null | sort | tail -1)
if [ -z "$R" ]; then
  ARGS="-i $IN time/tlim=$TLIM $XKEYS"
  NC=0
  echo "FRESH start from $IN"
else
  INF=$(python3 "$INFO" "$R") || { echo "rst_info failed on $R: exit 1"; exit 1; }
  read -r T NC TL0 <<< "$(echo "$INF" | head -1)"
  awk -v t="$T" -v l="$TLIM" 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "rst t=$T >= tlim: DONE"; touch DONE; exit 0; }
  ARGS="-r $R $(echo "$INF" | tail -1) time/tlim=$TLIM $XKEYS"
  echo "restart from $R (t = $T, cycle $NC)"
fi
[ -n "$NCYC" ] && ARGS="$ARGS time/nlim=$(( NC + NCYC )) time/ndiag=1"
[ -f run.log ] && mv run.log "run.$(date +%m%d_%H%M%S).log"
WRAP=$D/gpu_wrap.sh
printf '#!/bin/bash\nexport CUDA_VISIBLE_DEVICES=$SLURM_LOCALID\nexec "$@"\n' > "$WRAP"; chmod +x "$WRAP"
GDT=60; [ -n "$NCYC" ] && GDT=10
nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -l $GDT > "gpumem.$SLURM_JOB_ID.log" 2>&1 &
SP=$!
echo "srun $BIN $ARGS -d $D -t $TL  $(date +%T)"
srun -n "$SLURM_NTASKS" -c 72 --cpu-bind=cores "$WRAP" "$BIN" $ARGS -d "$D" -t "$TL" > run.log 2>&1
rc=$?
kill $SP 2>/dev/null
cp -p run.log "run.log.$SLURM_JOB_ID"
echo "rc=$rc fatal=$(grep -c FATAL run.log) nonconv_last=$(grep -o 'NON-CONVERGED=[0-9.e+-]*' run.log | tail -1)" \
     "gpumem_max_MiB=$(awk -F, '{if($2>m)m=$2} END{print m+0}' gpumem.$SLURM_JOB_ID.log) $(date +%F_%T)"
if grep -q "FATAL" run.log || tail -n 60 run.log | grep -qiw "nan"; then
  echo "NaN/FATAL in run.log: STOP"; echo "NaN/FATAL in run.log of job $SLURM_JOB_ID (rc $rc)" > STOP; exit 1
fi
[ $rc -ne 0 ] && { echo "srun rc=$rc: STOP"; echo "srun rc=$rc in job $SLURM_JOB_ID" > STOP; exit 1; }
R2=$(ls rst/*.rst 2>/dev/null | sort | tail -1)
if [ -n "$R2" ] && [ -z "$NCYC" ]; then
  T2=$(python3 "$INFO" "$R2" | head -1 | cut -d' ' -f1)
  awk -v t="$T2" -v l="$TLIM" 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "rst t=$T2 >= tlim: DONE"; touch DONE; }
fi
exit 0
