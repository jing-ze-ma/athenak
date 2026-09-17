#!/bin/bash -l
#SBATCH -o gpu_csctl.out.%j
#SBATCH -e gpu_csctl.err.%j
#SBATCH -J adi_csctl
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# CONTROL: is the degradation at 4x4 MeshBlocks per panel a property of the IMPLICIT
# transverse operator or of the cubed-sphere mesh machinery itself?  The explicit
# operator (rad_implicit_ang = false) is run on the same layouts, and the implicit ones
# are run with the metric cross term off (rad_cs_exact = false), which is the only term
# that reads the x2x3 DIAGONAL ghost.  nghost is varied too.
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
  (cd $d && srun -n 1 $B -i $I hydro/rad_ang_verbose=true "$@" > log.txt 2>&1)
  echo -n "$d  "; python3 measure_cs.py $d 2 2>/dev/null | sed "s/^[^ ]* *//"
  grep -o "|sum V de|/sum V|de| = [0-9.e+-]*" $d/log.txt | tail -1
  rm -rf $d/bin
}
for ng in 2 3; do
 for mb in 64 32 16; do
  run X64_exp_mb${mb}_ng$ng mesh/nghost=$ng mesh/nx2=64 mesh/nx3=64 \
      meshblock/nx2=$mb meshblock/nx3=$mb time/cfl_number=0.0075 time/nlim=80
  for s in sts adi; do
   run X64_${s}_mb${mb}_ng$ng mesh/nghost=$ng mesh/nx2=64 mesh/nx3=64 \
       meshblock/nx2=$mb meshblock/nx3=$mb time/cfl_number=0.0075 time/nlim=80 \
       hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
   run X64_${s}nox_mb${mb}_ng$ng mesh/nghost=$ng mesh/nx2=64 mesh/nx3=64 \
       meshblock/nx2=$mb meshblock/nx3=$mb time/cfl_number=0.0075 time/nlim=80 \
       hydro/rad_cs_exact=false hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
  done
 done
done
echo GPU_CSCTL_DONE
