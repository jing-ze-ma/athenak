#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r10/relax.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r10/relax.err
#SBATCH -J thickrelax
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
ROOT=/viper/u2/jinma/ATHENAK/bench/wt_he4
IN=$ROOT/tests_r10/two_stream_sph_thick_relax.athinput
EXE=$ROOT/build_gpu_rg/src/athena
cd $ROOT/tests_r10
mapfile -t CASES < <(grep -v '^#' cases.txt | grep -v '^[[:space:]]*$')
for row in "${CASES[@]}"; do
  read tag n tau ratio x1max kappa lstar finner <<< "$row"
  echo "=================== CASE $tag  n=$n tau=$tau rout/rin=$ratio"
  cd $ROOT/tests_r10/$tag
  rm -rf rst
  srun -n 1 $EXE -i $IN \
    mesh/x1max=$x1max problem/kappa_const=$kappa problem/lstar=$lstar \
    hydro/rad_flux_inner=$finner problem/rt_col3_skip_sweep=true \
    problem/mlt_dump=mltfaces_t0.txt time/nlim=200 \
    < /dev/null > run_stage1.log 2>&1
  echo "stage1 exit $?  $(tail -3 run_stage1.log)"
  RST=$(ls -t rst/*.rst | head -1)
  srun -n 1 $EXE -r $RST \
    problem/mlt_dump=mltfaces_relax.txt time/nlim=201 \
    < /dev/null > run_stage2.log 2>&1
  echo "stage2 exit $?  $(tail -3 run_stage2.log)"
  rm -rf rst *.bin
  cd $ROOT/tests_r10
done
