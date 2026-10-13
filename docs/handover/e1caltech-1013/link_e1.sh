#!/bin/bash -l
#SBATCH -J agcE1c
# e1caltech-1013 link_e1.sh RUN_DIR: one AG Car E1 link on 1 node x 2 H200 (2 MPI ranks x 2 MeshBlocks, 1 rank per GPU).
# = Raven rc_1009/jobs/rc_raven.sh (the stored agcmE1 31106649 script) + e1chain_hook.sh, ported. Submit with YOUR BSG link's
# sbatch flags (partition/account/gres for 2 H200, --ntasks=2, -t 12:00:00), e.g.
#   sbatch <your flags> -t 12:00:00 link_e1.sh <RUN_DIR>
# Everything swappable is read from RUN_DIR/run.cfg AT JOB START:
#   E1C_DIR BIN BIN_MD5 IN KEYS TLIM STOP_NC [E1C_NCYC (tests)] [E1C_GATE=1 E1C_VERIFIED=<smoke dir>/VERIFIED (production)]
#   [LAUNCH="<your srun line up to the binary>"]  (default below; use your BSG production launch if it differs)
# RUN_DIR/CANCEL -> exit 0. Restart = newest VALID rst (hook), never a fresh start over existing rst files.
export OMP_NUM_THREADS=1
D=$1; [ -n "$D" ] || { echo "no run dir"; exit 1; }
cd "$D" || exit 1
[ -f CANCEL ] && { echo "CANCEL exists: exit 0"; exit 0; }
[ -f run.cfg ] || { echo "no $D/run.cfg"; exit 1; }
BIN=; BIN_MD5=; IN=; RST=; KEYS=; LAUNCH=
source ./run.cfg
# E1C_DIR = the directory holding these scripts (set in run.cfg; under sbatch $0 is a spooled copy, so it cannot be derived)
[ -f "${E1C_DIR:-/nonexistent}/e1chain_hook.sh" ] || { echo "run.cfg must set E1C_DIR (dir of e1chain_hook.sh)"; exit 1; }
source "$E1C_DIR/e1chain_hook.sh"          # DONE / STOP / NOGO / GUARD exits; sets RST, KEYS; EXIT trap = STOP rules + prune
[ -x "$BIN" ] && [ -n "$BIN_MD5" ] && [ -f "$IN" ] || { echo "run.cfg incomplete"; exit 1; }
M=$(md5sum $BIN | cut -d' ' -f1); [ "$M" == "$BIN_MD5" ] || { echo "md5 $M != $BIN_MD5"; exit 1; }
cp -p run.cfg run.cfg.$SLURM_JOB_ID
TABS=$(grep -E '^ *(he_ic_file|he_opac_table|he_planck_table) ' $IN | awk '{print $3}')
md5sum $BIN $IN $TABS $RST
LIM=$(squeue -h -j $SLURM_JOB_ID -o %l)      # [D-]H:MM:SS
dd=0; [[ $LIM == *-* ]] && { dd=${LIM%%-*}; LIM=${LIM#*-}; }
IFS=: read -r a b c <<< "$LIM"; [ -z "$c" ] && { c=$b; b=$a; a=0; }
LIMS=$(( dd*86400 + 10#$a*3600 + 10#$b*60 + 10#$c ))
MARGIN=600; [ $LIMS -le 1200 ] && MARGIN=120
LEFT=$(( LIMS - SECONDS - MARGIN ))
TL=$(printf '%02d:%02d:%02d' $((LEFT/3600)) $(((LEFT%3600)/60)) $((LEFT%60)))
if [ -n "$RST" ]; then SRC="-r $RST"; else SRC="-i $IN"; fi
: "${LAUNCH:=srun -n 2 --cpu-bind=cores}"
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST rundir $D $SRC keys: $KEYS  -t $TL  start $(date +%F_%T)"
ln -sf run.log.$SLURM_JOB_ID run.log
$LAUNCH $BIN $SRC -d "$D" -t $TL $KEYS > run.log.$SLURM_JOB_ID 2>&1
rc=$?
echo "rc=$rc fatal=$(grep -c FATAL run.log.$SLURM_JOB_ID) $(date +%F_%T)"
