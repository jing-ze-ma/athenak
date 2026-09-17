#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r4/g3b_n0.log
#SBATCH -J g3b
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
EXE=$ROOT/tests_r4/athena_v2
IN3=$ROOT/inputs/tests/two_stream_sph_thick.athinput
mapfile -t CASES < <(grep -v '^#' $ROOT/tests_r2/thick/cases.txt | grep -v '^[[:space:]]*$')
for row in "${CASES[@]}"; do
  read tag n tau ratio x1max kappa lstar finner <<< "$row"
  D=$ROOT/tests_r4/g3_n0/$tag
  mkdir -p $D && cd $D
  ln -sf $ROOT/tests_r2/thick/$tag/ic_thick.txt ic_thick.txt
  srun -n 1 $EXE -i $IN3 mesh/x1max=$x1max problem/kappa_const=$kappa \
    problem/lstar=$lstar hydro/rad_flux_inner=$finner \
    problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false \
    < /dev/null > run.log 2>&1
  echo "=== G3 n0 $tag exit $?"
done
echo ALLDONE
