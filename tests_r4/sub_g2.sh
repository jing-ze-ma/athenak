#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r4/g2.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r4/g2.err
#SBATCH -J g2g3
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
EXE=$ROOT/tests_r4/athena_m3
IN=$ROOT/tests_m3/thin_ex.athinput
# ---- G2: the TRANSPARENT gate, L = 4 pi r^2 F must be constant ----
for arm in "m0:0:thomas:false" "m3s:3:thomas:true" "m3p:3:pcr:true"; do
  IFS=: read tag mode solv skip <<< "$arm"
  mkdir -p $ROOT/tests_r4/g2/$tag && cd $ROOT/tests_r4/g2/$tag
  srun -n 1 $EXE -i $IN problem/rt_implicit_column=$mode \
     problem/rt_impl_solver=$solv problem/rt_col3_skip_sweep=$skip \
     problem/rt_impl_mixed=0 problem/rt_impl_warm=0 \
     < /dev/null > run.log 2>&1
  echo "=== G2 $tag exit $?"
done
# ---- G3: the THICK gate ----
IN3=$ROOT/inputs/tests/two_stream_sph_thick.athinput
mapfile -t CASES < <(grep -v '^#' $ROOT/tests_r2/thick/cases.txt | grep -v '^[[:space:]]*$')
for md in 3; do
for row in "${CASES[@]}"; do
  read tag n tau ratio x1max kappa lstar finner <<< "$row"
  D=$ROOT/tests_r4/g3_m$md/$tag
  mkdir -p $D && cd $D
  srun -n 1 $EXE -i $IN3 mesh/x1max=$x1max problem/kappa_const=$kappa \
    problem/lstar=$lstar hydro/rad_flux_inner=$finner \
    problem/rt_implicit_column=$md \
    problem/rt_col3_skip_sweep=$( [ $md = 3 ] && echo true || echo false ) \
    < /dev/null > run.log 2>&1
  echo "=== G3 m$md $tag exit $? $(ls -la mltfaces.txt 2>/dev/null | awk '{print $5}')"
done
done
echo ALLDONE
