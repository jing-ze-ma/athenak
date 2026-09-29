#!/bin/bash
# bench-2026-09-29-hebox cross-cluster timing: ONE fresh-start 3-D He-box run (implicit M1 +
# vet_sc, hesdirk2, cfl 0.3) of time/nlim cycles (input: 1200).
#
#   run_bench_hebox.sh <athena binary (PROBLEM=box_convection)> <nranks> <run dir>
#
# Environment (set by the machine header below or by the caller):
#   MACHINE      viper | caltech | deltaai   (selects the header; default none, then LAUNCH needed)
#   HE_BOX_DATA  directory with the five he_box/ files of athenak_data_2026-09-26.tar.gz
#                (ic_m1_V3edd_pgen.txt m1_rad_ic_V3edd.txt rosseland_he_x0.0_z0.02.txt
#                 planck_he_x0.0_z0.02_ferg+tops.txt arad_V3edd.txt)
#   LAUNCH       MPI launcher prefix; "%N" is replaced by <nranks>.  The header sets a default.
#   EXTRA        extra block/key=value overrides (smoke only; do not change physics for timing)
# Writes <run dir>/bench.log (athena stdout+stderr after two "BENCH" header lines) and the hst.
# Analyse with: python3 ana_bench_hebox.py <run dir>/bench.log ...
# nranks must be 1, 2 or 4 (the box has 4 MeshBlocks of 84 x 52 x 52).
set -u
BIN=$1; NR=$2; D=$3
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MACHINE=${MACHINE:-}
case $NR in 1|2|4) ;; *) echo "nranks must be 1, 2 or 4 (4 MeshBlocks)"; exit 2 ;; esac

# ---------------------------------------------------------------- machine headers
type module > /dev/null 2>&1 || source /etc/profile > /dev/null 2>&1   # module in a non-login shell
case $MACHINE in
  viper)
    # MPCDF viper MI300A (apu / apudev): 2 GPUs per node, 1 rank per GPU, 24 cores per rank.
    # Job header: -p apudev (<15 min) or -p apu, --constraint=apu --gres=gpu:2
    #   --ntasks-per-node=2 --cpus-per-task=24; 4 GPUs = -p apu -N 2.  SBATCH_EXPORT=NONE.
    # Stack: ROCm 7.2 (builds/build_inc_viper.sh target box_gpu72, adopted 09-29).
    module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
    export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1          # validated ROCm settings (1.32x)
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-24}
    HE_BOX_DATA=${HE_BOX_DATA:-/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/he_box}
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
    module purge
    module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
    HE_BOX_DATA=${HE_BOX_DATA:-/resnick/home/jingze/ATHENAK/athenak_data/he_box}
    LAUNCH=${LAUNCH:-"srun --mpi=pmix -n %N -c ${SLURM_CPUS_PER_TASK:-8} --gpus-per-task=1 --exact"} ;;
  deltaai)
    # NCSA DeltaAI GH200 (4 per node, 72-core Grace each; GPU r <-> NUMA r <-> cores 72r..72r+71).
    # Job header: -A bivj-dtai-gh -p ghx4 (or ghx4-interactive, x2 charge) -N 1 --ntasks-per-node=<n>
    #   --gpus-per-node=<n> --cpus-per-task=16.  Default modules (PrgEnv-gnu, gcc-native/14,
    #   cudatoolkit/25.5_12.9, cray-mpich/9.0.1, craype-accel-nvidia90).  Same header as bench-2026-09-28.
    # Launch = the AthenaK docs' DeltaAI recipe: every rank sees all GPUs, Kokkos maps device by local rank,
    #   cores bound per rank.  NOT --gpus-per-task=1/--gpu-bind=closest (cgroups block CUDA IPC).
    module reset > /dev/null 2>&1
    export MPICH_GPU_SUPPORT_ENABLED=1 OMP_NUM_THREADS=1
    export SLURM_CPU_BIND=cores KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
    HE_BOX_DATA=${HE_BOX_DATA:-/u/jma20/ATHENAK/data/athenak_data/he_box}
    LAUNCH=${LAUNCH:-"srun -n %N -c ${SLURM_CPUS_PER_TASK:-16} --cpu-bind=cores"} ;;
  *)
    : "${HE_BOX_DATA:?set HE_BOX_DATA or MACHINE}" "${LAUNCH:?set LAUNCH or MACHINE}" ;;
esac
# --------------------------------------------------------------------------------

for f in ic_m1_V3edd_pgen.txt m1_rad_ic_V3edd.txt rosseland_he_x0.0_z0.02.txt \
         planck_he_x0.0_z0.02_ferg+tops.txt arad_V3edd.txt; do
  [ -r "$HE_BOX_DATA/$f" ] || { echo "missing $HE_BOX_DATA/$f"; exit 2; }
done
mkdir -p "$D"
if [ -e "$D/bench.log" ]; then echo "$D/bench.log exists: use a new run dir"; exit 2; fi
sed -e "s#HE_BOX_DATA#$HE_BOX_DATA#g" "$HERE/inputs/hebox_bench.athinput" > "$D/run.athinput"
L=${LAUNCH//%N/$NR}
MD5=$(md5sum "$BIN" | cut -c1-32)
{
  echo "BENCH arm=c03 nranks=$NR machine=${MACHINE:-custom} host=$(hostname) job=${SLURM_JOB_ID:-none}" \
       "bin=$BIN md5=$MD5 date=$(date -Iseconds)"
  echo "BENCH launch: $L $BIN -i run.athinput ${EXTRA:-}"
} > "$D/bench.log"
( cd "$D" && $L "$BIN" -i run.athinput ${EXTRA:-} >> bench.log 2>&1 )
RC=$?
echo "BENCH rc=$RC end=$(date -Iseconds)" >> "$D/bench.log"
echo "run_bench_hebox n=$NR rc=$RC $D"
exit $RC
