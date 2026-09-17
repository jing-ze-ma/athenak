#!/bin/bash -l
#SBATCH -o gpu_csctl2.out.%j
#SBATCH -e gpu_csctl2.err.%j
#SBATCH -J adi_csctl2
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# Is the 4x4-MeshBlocks-per-panel degradation a function of the BLOCK COUNT or of the
# BLOCK SIZE?  n = 128 per panel edge, 1x1 / 2x2 / 4x4 blocks per panel (128/64/32 cells),
# explicit and both implicit solvers.  The production grid is 4x4 blocks of 80 cells.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
B=$W/build_gpu_cs_test/src/athena
I=$W/inputs/tests/cubed_sphere_raddiff.athinput
run(){ d=$1; shift
  rm -rf $d; mkdir -p $d
  (cd $d && srun -n 2 $B -i $I hydro/rad_ang_verbose=true "$@" > log.txt 2>&1)
  echo -n "$d  "; python3 measure_cs.py $d 2 2>/dev/null | sed "s/^[^ ]* *//"
  grep -o "|sum V de|/sum V|de| = [0-9.e+-]*" $d/log.txt | tail -1
  rm -rf $d/bin
}
for mb in 128 64 32; do
  run Y128_exp_mb$mb mesh/nx2=128 mesh/nx3=128 meshblock/nx2=$mb meshblock/nx3=$mb \
      time/cfl_number=0.0075 time/nlim=80
  for s in sts adi; do
    run Y128_${s}_mb$mb mesh/nx2=128 mesh/nx3=128 meshblock/nx2=$mb meshblock/nx3=$mb \
        time/cfl_number=0.0075 time/nlim=80 \
        hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
  done
done
echo GPU_CSCTL2_DONE
