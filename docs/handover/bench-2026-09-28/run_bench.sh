#!/bin/bash
# bench-2026-09-28 cross-cluster timing: ONE fresh-start WASP-121b run of 2000 cycles.
#
#   run_bench.sh <athena binary> <nranks> <1x|10x> <run dir>
#
# Environment (set by the machine header below or by the caller):
#   MACHINE   viper | caltech | deltaai   (selects the header; default: none, then LAUNCH needed)
#   CKDATA    directory that replaces ATHENAK_CK_DATA   (1x tables: repo data/exo_fms_ck)
#   CKDATA10  directory that replaces ATHENAK_CK_DATA10 (10x tables: ckdata10, with cia/ ray/
#             sw_flux/ as in the Caltech handover)
#   LAUNCH    MPI launcher prefix; "%N" is replaced by <nranks>.  The header sets a default.
#   EXTRA     extra block/key=value overrides (default none; do not change physics for timing)
# Writes <run dir>/bench.log (athena stdout+stderr, preceded by one "BENCH ..." header line)
# and <run dir>/dhj.hydro.hst.  Analyse with: python3 ana_bench.py <run dir>/bench.log ...
# The run dir must not exist or must be empty of an old bench.log (the script refuses to reuse).
set -u
BIN=$1; NR=$2; ARM=$3; D=$4
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MACHINE=${MACHINE:-}

# ---------------------------------------------------------------- machine headers
type module > /dev/null 2>&1 || source /etc/profile > /dev/null 2>&1   # module in a non-login shell
case $MACHINE in
  viper)
    # MPCDF viper MI300A (apu / apudev): 2 GPUs per node, 1 rank per GPU, 24 cores per rank.
    # Job header: -p apudev (<15 min) or -p apu, --constraint=apu --gres=gpu:2
    #   --ntasks-per-node=2 --cpus-per-task=24; 4 GPUs = -p apu -N 2.  SBATCH_EXPORT=NONE on
    #   viper: the job loads its own modules.
    module purge; module load gcc/14 rocm/6.3 openmpi_gpu/5.0
    export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1          # validated ROCm settings (1.32x)
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-24}
    CKDATA=${CKDATA:-/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck}
    CKDATA10=${CKDATA10:-/viper/ptmp2/jinma/wasp121_0925/ckdata10}
    # 1 rank inside a 2-GPU allocation must take one GPU only (--exact, --gres=gpu:1)
    if [ "$NR" = 1 ]; then
      LAUNCH=${LAUNCH:-"srun -n 1 -N 1 -c ${SLURM_CPUS_PER_TASK:-24} --gres=gpu:1 --exact"}
    else
      LAUNCH=${LAUNCH:-"srun -n %N -c ${SLURM_CPUS_PER_TASK:-24}"}
    fi ;;
  caltech)
    # Caltech resnick gpu partition: 4 H200 per node.  TIME ON H200 ONLY
    #   (--gres=gpu:nvidia_h200:<n>), exclude the bad nodes:
    #   #SBATCH -A carnegie_poc -p gpu --exclude=hpc-sm-01-09,hpc-sm-02-16
    #   #SBATCH -N 1 --ntasks-per-node=<n> --gres=gpu:nvidia_h200:<n> --cpus-per-task=8 --mem=100G
    #   An H100 row may be added only if clearly labelled H100.
    module purge
    module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
    CKDATA=${CKDATA:-/resnick/home/jingze/ATHENAK/athenak_data/exo_fms_ck}
    CKDATA10=${CKDATA10:-/resnick/home/jingze/ATHENAK/athenak_data/ckdata10}
    LAUNCH=${LAUNCH:-"srun --mpi=pmix -n %N -c ${SLURM_CPUS_PER_TASK:-8} --gpus-per-task=1 --exact"} ;;
  deltaai)
    # NCSA DeltaAI GH200 (4 per node, ARM Grace host).  TO BE FILLED by the DeltaAI session
    # (modules, account, partition, launcher, data dirs) per TASK-2026-09-28-deltaai-bringup.md.
    : "${CKDATA:?set CKDATA}" "${CKDATA10:?set CKDATA10}" "${LAUNCH:?set LAUNCH}" ;;
  *)
    : "${CKDATA:?set CKDATA or MACHINE}" "${CKDATA10:?set CKDATA10 or MACHINE}" \
      "${LAUNCH:?set LAUNCH or MACHINE}" ;;
esac
# --------------------------------------------------------------------------------

case $ARM in
  1x)  INP=$HERE/inputs/w121_bench_1x.athinput;  IC=$HERE/inputs/ic_w121_1x.txt ;;
  10x) INP=$HERE/inputs/w121_bench_10x.athinput; IC=$HERE/inputs/ic_w121.txt ;;
  *) echo "arm must be 1x or 10x"; exit 2 ;;
esac
mkdir -p "$D"
if [ -e "$D/bench.log" ]; then echo "$D/bench.log exists: use a new run dir"; exit 2; fi
sed -e "s#ATHENAK_CK_DATA10#$CKDATA10#g" -e "s#ATHENAK_CK_DATA#$CKDATA#g" "$INP" > "$D/run.athinput"
cp "$IC" "$D/"
L=${LAUNCH//%N/$NR}
MD5=$(md5sum "$BIN" | cut -c1-32)
{
  echo "BENCH arm=$ARM nranks=$NR machine=${MACHINE:-custom} host=$(hostname) job=${SLURM_JOB_ID:-none}" \
       "bin=$BIN md5=$MD5 date=$(date -Iseconds)"
  echo "BENCH launch: $L $BIN -i run.athinput ${EXTRA:-}"
} > "$D/bench.log"
( cd "$D" && $L "$BIN" -i run.athinput ${EXTRA:-} >> bench.log 2>&1 )
RC=$?
echo "BENCH rc=$RC end=$(date -Iseconds)" >> "$D/bench.log"
echo "run_bench $ARM n=$NR rc=$RC $D"
exit $RC
