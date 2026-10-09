#!/bin/bash
#SBATCH -J radbeam
#SBATCH -A bivj-delta-gpu
#SBATCH -p gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH -o /work/nvme/bivj/jma20/delta_1008/radbeam/j.%j.out
# rad-beam-1008 AG Car A 3-D arms (TASK-2026-10-08-delta-beam-batchN). 1 node x 4 A100, 4 ranks.
# usage: sbatch radbeam_delta.sh <binary> <armsfile> <arm name> [wall hh:mm:ss]
#   runs ONE arm (the line of <armsfile> whose first word is <arm name>) fresh from t = 0 for
#   <wall> (default 00:13:00), bins (hydro_w, m1, m1_fs) every 100 cycles and at the end;
#   exits non-zero if the run fails (rc or FATAL).  arm line: <name> <keys...>
# env (tests only): W (run root), IN0 (production A input), SRUN (launcher)
W=${W:-/work/nvme/bivj/jma20/delta_1008/radbeam}
IN0=${IN0:-/work/nvme/bivj/jma20/delta_1008/agcar_files/agcar_shakeA_ge.athinput}
SRUN=${SRUN:-srun -n 4 --cpu-bind=cores}
BIN=$1; ARMS=$2; ARM=$3; WALL=${4:-00:13:00}
[ -x "$BIN" ] || { echo "no binary $BIN"; exit 1; }
[ -f "$IN0" ] || { echo "no input $IN0"; exit 1; }
module unload cudatoolkit 2>/dev/null
module load cuda/12.9 2>/dev/null
export MPICH_GPU_SUPPORT_ENABLED=1
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
export OMP_NUM_THREADS=1
export HSA_XNACK=1; export HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
JID=${SLURM_JOB_ID:-test$$}
awk -v a="$ARM" '$1 == a {f = 1} END {exit !f}' "$ARMS" || { echo "arm $ARM not in $ARMS"; exit 1; }
keys=$(awk -v a="$ARM" '$1 == a {$1 = ""; print; exit}' "$ARMS")
echo "job $JID node ${SLURMD_NODENAME:-local} arm $ARM keys:$keys wall $WALL"
# per-job COPY of the input holding every key the command line overrides (a fresh run
# refuses absent keys and blocks): dcycle in output1-3, the rad-beam keys in <rad_m1>
# (defaults = production arithmetic), and output6 = m1_fs bin
mkdir -p "$W/in"
IN=$W/in/agcarA.$ARM.$JID.athinput
awk '{print} /^<output1>/{print "dcycle = 1"} /^<output2>/||/^<output3>/{print "dcycle = 100"}
  /^<rad_m1>/{print "vet_source_noesrc = false"; print "implicit_recon = dc";
              print "implicit_recon_dgpass = false"}' "$IN0" > "$IN"
{ echo "<output6>"; echo "file_type = bin"; echo "variable = m1_fs"; echo "dcycle = 100"; } >> "$IN"
md5sum "$BIN" "$IN0" "$IN"
BASE="time/ndiag=1 output4/dt=1e30 output5/dt=1e30"
D=$W/A.$ARM.$JID
mkdir -p "$D/bin"
cd "$D" || exit 1
echo "start $(date +%T)"
$SRUN "$BIN" -i "$IN" -d "$D" -t "$WALL" $BASE $keys < /dev/null > run.log 2>&1
rc=$?
nf=$(grep -c FATAL run.log)
echo "rc=$rc fatal=$nf $(date +%T)"
ls "$D/bin" | tail -4
if [ $rc -ne 0 ] || [ "$nf" -ne 0 ]; then exit 2; fi
exit 0
