#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/h2div_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/h2div_0924/gpu/log.err.%j
#SBATCH -J h2div
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-h2div: the diffuse-wall shadow (vet_sc, n=64, t=20) on the GPU, base vs f2, hesdirk2 and be,
# interleaved, 2 repeats; main-loop seconds from "cpu time used"
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/h2div_0924
md5sum $W/bin/athena_base_none_gpu $W/bin/athena_f2_none_gpu
for rep in 1 2; do
  for arm in f2:h2 base:h2 f2:be base:be; do
    b=${arm%%:*}; s=${arm#*:}; d=$W/gpu/runs/${b}_${s}_r$rep; rm -rf $d; mkdir -p $d; cd $d
    srun -n 1 $W/bin/athena_${b}_none_gpu -i $W/inp/sh_$s.athinput -d . -t 00:02:30 > log.txt 2>&1
    echo "#### $b $s rep $rep rc=$? $(grep -h 'cpu time used' log.txt) $(grep -h 'stage fallbacks' log.txt | sed 's/.*\(stage fallbacks=[^ ]*\).*/\1/') $(grep -h 'NON-CONVERGED=' log.txt | sed 's/.*\(NON-CONVERGED=[^ ]*\).*/\1/')"
  done
done
echo GPU DONE
