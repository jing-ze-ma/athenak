#!/bin/bash -l
# bsgli-val-1014 val.sh (raven -> DeltaAI, TASK-2026-10-14-deltaai-bsgli-val): VALIDATION ONLY of fin4-linhao-bc on 1 node x 4 GH200.
# Full 256^3 BSG, LI deck + problem/he_bc_linhao = true, 4 ranks x 4 MeshBlocks (256x64x64), 1 rank per GPU.
#   C    : fresh, 60 cycles
#   L1/L2: fresh 30 cycles (ends with the link-end rst) -> NEW srun restarted from rst_pick.py's pick (+ its last_time keys) to 60
#   Coff : C with rad_m1/vet_gd_halo_ipc=false (speed-only key; must be BITWISE vs C)
# Gate: chaincmp.py C L1 L2 (every hst row + every scheduled dump, data arrays) and cmplast.py C L2 / C Coff -> BITWISE.
# Adapt the #SBATCH lines and LAUNCH to your usual BSG setup (account, partition, gpuwrap.sh). Usage:
#   sbatch --export=ALL,W=<dir of these scripts>,BIN=<binary>,DECK=<resolved deck>,S=<scratch dir> val.sh
#SBATCH -J bsgli_val
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:45:00
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1
W=${W:?W}; BIN=${BIN:?BIN}; DECK=${DECK:?DECK}; S=${S:?S}/j.${SLURM_JOB_ID:-dry}; PY=${PY:-python3}
: "${LAUNCH:=srun -n 4 --ntasks-per-node=4 --cpu-bind=cores}"   # use your BSG launch (e.g. ... -c 72 --cpu-bind=cores <gpuwrap.sh>)
mkdir -p $S; echo "job ${SLURM_JOB_ID:-} ${SLURM_JOB_NODELIST:-} start $(date)"
md5sum $BIN $DECK $(grep -E '^ *(he_ic_file|he_opac_table|he_planck_table|he_ic_linhao_file) ' $DECK | awk '{print $3}')
OC="time/ndiag=1 output1/dt=200 output2/dt=1000 output3/dt=1000 output4/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30"
run() { local n=$1; shift; local D=$S/$n; mkdir -p $D; cd $D; echo "$*" > cmd.txt; local t0=$SECONDS
  $LAUNCH $BIN "$@" -d $D < /dev/null > run.log 2>&1; local rc=$?
  echo "ARM $n rc=$rc wall=$((SECONDS-t0))s fatal=$(grep -c FATAL run.log) nan=$(tail -n 60 run.log | grep -ciw nan) diverged=$(grep -c DIVERGED run.log)" \
       "| $(grep -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' run.log | tail -1)"
  grep -m1 -E 'vet_gd_halo_ipc:' run.log; grep -m2 -E 'he_bc_linhao' run.log; grep -m2 -E 'he_ic_linhao' run.log
  grep -E 'cycle=' run.log | awk '{for(i=1;i<=NF;i++) if($i ~ /^dt=/) print $i}' | sort | uniq -c | head -5
  nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader 2>/dev/null; }
run C   -i $DECK $OC time/nlim=60
run L1  -i $DECK $OC time/nlim=30
PK=$($PY $W/rst_pick.py $S/L1/rst); echo "PICK $PK"
run L2  -r $(echo "$PK" | head -1 | cut -d' ' -f1) $(echo "$PK" | tail -1) time/nlim=60
run Coff -i $DECK $OC time/nlim=60 rad_m1/vet_gd_halo_ipc=false
export BINCONV_DIR=$W
echo "== chaincmp C L1 L2";  $PY $W/chaincmp.py $S/C $S/L1 $S/L2 | tee $S/gate_chaincmp.txt | tail -5
echo "== cmplast C L2";      $PY $W/cmplast.py $S/C $S/L2   | tee $S/gate_cmplast.txt  | tail -1
echo "== cmplast C Coff";    $PY $W/cmplast.py $S/C $S/Coff | tee $S/gate_coff.txt     | tail -1
for h in $(cd $S/C; ls *.hst); do cmp -s <(grep -v '^#' $S/C/$h) <(grep -v '^#' $S/Coff/$h) && echo "COFF HST $h BITWISE" || echo "COFF HST $h DIFF"; done
echo "BSGLI_VAL_DONE $S $(date)"
