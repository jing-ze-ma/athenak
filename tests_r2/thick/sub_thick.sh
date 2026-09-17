#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r2/thick/thick.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r2/thick/thick.err
#SBATCH -J rgthick
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
IN=$ROOT/inputs/tests/two_stream_sph_thick.athinput
EXE=$ROOT/build_gpu_rg/src/athena
cd $ROOT/tests_r2/thick
# NOTE: the case loop must NOT be fed by a pipe -- srun inherits stdin and swallows the
# rest of the list, so only the first case ever ran.  Read the file into an array first,
# and give srun its own /dev/null.
mapfile -t CASES < <(grep -v '^#' cases.txt | grep -v '^[[:space:]]*$')
for row in "${CASES[@]}"; do
  read tag n tau ratio x1max kappa lstar finner <<< "$row" 
  echo "=================== CASE $tag  n=$n tau=$tau rout/rin=$ratio"
  cd $ROOT/tests_r2/thick/$tag
  srun -n 1 $EXE -i $IN \
    mesh/x1max=$x1max \
    problem/kappa_const=$kappa \
    problem/lstar=$lstar \
    hydro/rad_flux_inner=$finner \
    < /dev/null > run.log 2>&1
  echo "exit $?  $(ls -la mltfaces.txt 2>/dev/null)"
  tail -5 run.log
  cd $ROOT/tests_r2/thick
done
