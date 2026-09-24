#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/h2fast_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/h2fast_0924/gpu/log.err.%j
#SBATCH -J h2fast
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-h2fast timing: 3-D He box 84x104x104 (4 blocks), 1 GPU, 120 cycles, ms/cycle 20-120.
# Env: ARMS = "name|binary-tag|input|overrides;..."  REP = repeats (even repeats reversed)
# PROF=1: rocprofv3 kernel trace, 40 cycles, one pass over the arms
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/h2fast_0924
P="problem/vpert=1.0e-2"
IFS=';' read -ra AL <<< "$ARMS"
n=${#AL[@]}
for r in $(seq 1 ${REP:-2}); do
  ORD=$(seq 0 $((n-1))); [ $((r % 2)) -eq 0 ] && ORD=$(seq $((n-1)) -1 0)
  for x in $ORD; do
    IFS='|' read -r name tag inp ov <<< "${AL[$x]}"
    E=$W/bin/athena_${tag}_gpu
    d=$W/gpu/${name}_$r; rm -rf $d; mkdir -p $d; cd $d
    echo "#### ${name}_$r $(md5sum $E | cut -c1-8) $(date +%T)"
    if [ "${PROF:-0}" == 1 ]; then
      rocprofv3 --kernel-trace --stats --output-format csv -d $d/prof -o rank_0 -- \
        $E -i $W/inp/$inp.athinput -d . -t 00:03:00 time/nlim=40 rad_m1/implicit_picard_log=12 $P $ov > run.log 2> run.err
      echo "rc=$?"; python3 $W/tools/split.py $d/prof 0 40 5 > split.txt 2>&1
    else
      srun -n 1 $E -i $W/inp/$inp.athinput -d . -t 00:01:40 time/nlim=120 $P $ov > run.log 2> run.err
      echo "rc=$? $(date +%T)"
    fi
  done
done
python3 $W/gsum.py gpu
