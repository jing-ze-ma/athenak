#!/bin/bash -l
#SBATCH -o gpu_cs.out.%j
#SBATCH -e gpu_cs.err.%j
#SBATCH -J adi_cs_gpu
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# The cubed-sphere arm of the implicit transverse operator, on 2 MI300A APUs (so the
# panel seams are also RANK boundaries and the module's own exchange is exercised over
# MPI).  Same runs as the CPU ladder, at n = 16 and 32.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
B=$W/build_gpu_cs_test/src/athena
I=$W/inputs/tests/cubed_sphere_raddiff.athinput
run(){ d=$1; shift; rm -rf $d; mkdir -p $d; (cd $d && srun -n 2 $B -i $I "$@" > log.txt 2>&1); }
for p in "16 0.03 20" "32 0.0075 80"; do
  set -- $p; n=$1; c=$2; l=$3
  for s in exp sts adi; do
    extra=""
    [ $s = sts ] && extra="hydro/rad_implicit_ang=true hydro/rad_ang_solver=sts"
    [ $s = adi ] && extra="hydro/rad_implicit_ang=true hydro/rad_ang_solver=adi"
    run G_${s}_$n mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$n meshblock/nx3=$n \
        time/cfl_number=$c time/nlim=$l hydro/rad_ang_verbose=true $extra
  done
done
# the stability arm: transverse diffusion number ~ 2e3 and ~7e5
for kf in 3.0e-3 1.0e-5; do
  for s in sts adi; do
    mx=200; [ $s = sts ] && mx=2000
    run GS_${s}_$kf mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32 \
        time/nlim=40 hydro/rad_kappa_fac=$kf hydro/rad_implicit_x1=true \
        hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s hydro/rad_ang_maxit=$mx \
        hydro/rad_ang_verbose=true
  done
done
echo GPU_CS_DONE
