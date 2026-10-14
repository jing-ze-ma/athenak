#!/bin/bash -l
#SBATCH -J bsgLIc
# bsgli-caltech-1014 link_bsgli.sh RUN_DIR: one link of the Linhao-IC BSG production on 1 node x 2 H200 (2 MPI ranks x 8 MeshBlocks
# 256x64x64, 1 rank per GPU). = e1caltech-1013/link_e1.sh (your E1 link) with the hook bsgli_hook.sh and the deck's 4 file paths
# in the md5 line. Submit with your BSG link's sbatch flags (1 node x 2 H200, --ntasks=2, -t 12:00:00):
#   sbatch <your flags> --ntasks=2 -t 12:00:00 -J bsgLIc link_bsgli.sh <RUN_DIR>
# Everything swappable is read from RUN_DIR/run.cfg AT JOB START:
#   BLC_DIR BIN BIN_MD5 IN KEYS TLIM STOP_NC [BLC_NCYC (tests)] [BLC_GATE=1 BLC_VERIFIED=<smoke dir>/VERIFIED (production)]
#   [LAUNCH=<your srun line up to the binary>]  (use your BSG/E1 1-node launch, e.g. e1_1013/e1launch.sh)
# RUN_DIR/CANCEL -> exit 0. Restart = newest VALID rst (hook), never a fresh start over existing rst files.
export OMP_NUM_THREADS=1
D=$1; [ -n "$D" ] || { echo "no run dir"; exit 1; }
cd "$D" || exit 1
[ -f CANCEL ] && { echo "CANCEL exists: exit 0"; exit 0; }
[ -f run.cfg ] || { echo "no $D/run.cfg"; exit 1; }
BIN=; BIN_MD5=; IN=; RST=; KEYS=; LAUNCH=
source ./run.cfg
# BLC_DIR = the directory holding these scripts (set in run.cfg; under sbatch $0 is a spooled copy, so it cannot be derived)
[ -f "${BLC_DIR:-/nonexistent}/bsgli_hook.sh" ] || { echo "run.cfg must set BLC_DIR (dir of bsgli_hook.sh)"; exit 1; }
source "$BLC_DIR/bsgli_hook.sh"          # DONE / STOP / NOGO / GUARD exits; sets RST, KEYS; EXIT trap = STOP rules + prune
[ -x "$BIN" ] && [ -n "$BIN_MD5" ] && [ -f "$IN" ] || { echo "run.cfg incomplete"; exit 1; }
M=$(md5sum $BIN | cut -d' ' -f1); [ "$M" == "$BIN_MD5" ] || { echo "md5 $M != $BIN_MD5"; exit 1; }
cp -p run.cfg run.cfg.$SLURM_JOB_ID
TABS=$(grep -E '^ *(he_ic_file|he_opac_table|he_planck_table|he_ic_linhao_file) ' $IN | awk '{print $3}')
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
