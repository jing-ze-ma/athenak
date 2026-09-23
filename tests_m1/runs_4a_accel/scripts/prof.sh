#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/accel_0923/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/accel_0923/log.err.%j
#SBATCH -J acprof
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# rocprofv3 kernel trace, 3-D He box, 40 cycles.  Env: BIN, ARMS ("name|input|overrides;...")
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/accel_0923
md5sum $BIN
IFS=';' read -ra AL <<< "$ARMS"
for a in "${AL[@]}"; do
  IFS='|' read -r name inp ov <<< "$a"
  d=$W/runs/$name; rm -rf $d; mkdir -p $d; cd $d
  echo "#### $name $(date +%T)"
  rocprofv3 --kernel-trace --stats --output-format csv -d $d/prof -o rank_0 -- \
    $BIN -i $W/inp/$inp.athinput -d . -t 00:04:00 time/nlim=40 rad_m1/implicit_picard_log=12 $ov > run.log 2> run.err
  echo "rc=$? $(date +%T)"; grep -E "Picard iterations|inner iterations" run.log
  python3 $W/tools/split.py $d/prof 0 40 5 > split.txt; cat split.txt
done
